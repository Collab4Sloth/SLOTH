/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief 2D coalescence bubbles solved by Cahn-Hilliard equations
 * @version 0.1
 * @date 2025-07-04
 *
 * Copyright CEA (c) 2025
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

// We need this class for test case sources
struct TestParameters {
  int order = 1;
  int nx = 32;
  int ny = 32;
  int nz = 32;
  int refinement = 3;
  int verbosity = -1;
  int post_processing = 1; // default value : enabled
  double duration = 5.0;
  double dt = 5.e-2;
  int tcase = 0;
};


void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
  args.AddOption(&p.order, "-o", "--order", "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement", "refinement level of the mesh, default = 3");
  args.AddOption(&p.nx, "-nx", "--nx", "number of elements in dimension X, default = 32");
  args.AddOption(&p.ny, "-ny", "--ny", "number of elements in dimension Y, default = 32");
  args.AddOption(&p.nz, "-nz", "--nz", "number of elements in dimension Z, default = 32");
  args.AddOption(&p.dt, "-dt", "--delta-t", "timestep incriment, default = 0.05");
  args.AddOption(&p.duration, "-d", "--duration", "timestep incriment, default = 5s");
  args.AddOption(&p.verbosity, "-v", "--verbosity", "verbosity level, default is -1");
  args.AddOption(&p.post_processing, "-p", "--post-processing", "run post processing step");
  args.AddOption(&p.tcase, "-tc", "--test-case", "0: {GMRES, ILU}, 1: {PCG, GAMG}");

  args.Parse();

  if (!args.Good()) {
    if (mfem::Mpi::WorldRank() == 0) {
      args.PrintUsage(mfem::out);
      std::exit(EXIT_FAILURE);
    }
  }
  if(mfem::Mpi::WorldRank() == 0) args.PrintOptions(mfem::out);
}

int64_t sum(int64_t in) {
  int64_t res = 0;
  MPI_Reduce(&in, &res, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
};

template <typename Arg>
void Message(Arg a_msg) {
  if (mfem::Mpi::WorldRank() == 0) {
    mfem::out << a_msg << std::endl;
  }
}

template <typename Arg, typename... Args>
void Message(Arg a_msg, Args... a_msgs) {
  if (mfem::Mpi::WorldRank() == 0) {
    mfem::out << a_msg << " ";
    Message(a_msgs...);
  }
}

template <typename Mesh, typename FES>
void print_mesh_information(Mesh& mesh, FES& fespace) {

  // get the number of vertices
  int64_t numbers_of_vertices_local = mesh.GetNV();
  int64_t numbers_of_vertices = sum(numbers_of_vertices_local);

  // get the number of elements
  int64_t numbers_of_elements_local = mesh.GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);

  // get n dofs
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);

  Message("INFO: number of vertices -> ", numbers_of_vertices);
  Message("INFO: number of elements -> ", numbers_of_elements);
  Message("INFO: number of dofs     -> ", unknowns);
}

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Verbose);

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 3;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using NLFI = CahnHilliardNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
  ThermodynamicsPotentials::F, Mobility::Constant>;

  using LHS_NLFI = TimeCHNLFormIntegrator<VARS>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, LHS_NLFI>;
  using PB = Problem<OPE, VARS, PST>;
  using PB1 = MPI_Problem<VARS, PST>;

  // ################ //
  // ################ //
  //   Read options   //
  // ################ //
  // ################ //
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string mesh_type =
    "InlineSquareWithHexaedres";  // type of mesh // "InlineSquareWithTriangles"
  const int order_fe =  p.order;           // finite element order
  const int refinement_level = p.refinement;   // number of levels of uniform refinement
  const int nx = p.nx;
  const int ny = p.ny;
  const int nz = p.nz;
  const double lx = 2. * M_PI;
  const double ly = 2. * M_PI;
  const double lz = 2. * M_PI;
  const std::tuple<int, int, int, double, double, double>& tuple_of_dimensions = std::make_tuple(
      nx, ny, nz, lx, ly, ly);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  print_mesh_information(*(spatial.get_mesh()), *(spatial.get_finite_element_space()));

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("bottom", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.),
    Boundary("top", 2, "Neumann", 0.),    Boundary("left", 3, "Neumann", 0.),
    Boundary("front", 4, "Neumann", 0.),  Boundary("rear", 5, "Neumann", 0.)};
  auto bcs_phi = BCS(&spatial, boundaries);
  auto boundaries_mu = {Boundary("bottom", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.),
    Boundary("top", 2, "Neumann", 0.),    Boundary("left", 3, "Neumann", 0.),
    Boundary("front", 4, "Neumann", 0.),  Boundary("rear", 5, "Neumann", 0.)};
  auto bcs_mu = BCS(&spatial, boundaries_mu);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Interface thickness
  const double epsilon(0.02);
  // Interfacial energy
  const double sigma(1.);
  // Two-phase mobility
  const double mob(1.);
  const double lambda = (epsilon * epsilon);
  const double omega = 1.;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
      Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################

  auto user_func_solution =
    std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
  const double xx = x[0];
  const double yy = x[1];
  const double zz = x[2];
  const double r1 = (xx - M_PI + 1) * (xx - M_PI + 1) + (yy - M_PI) * (yy - M_PI) +
  (zz - M_PI) * (zz - M_PI);
  const double r2 = (xx - M_PI - 1) * (xx - M_PI - 1) + (yy - M_PI) * (yy - M_PI) +
  (zz - M_PI) * (zz - M_PI);
  const double r3 = (xx - M_PI - 1) * (xx - M_PI - 1) + (yy - M_PI) * (yy - M_PI) +
  (zz - M_PI) * (zz - M_PI);
  double sol = 0.;
  if (r1 < 1 || r2 < 1 || r3 < 1) {
  sol = 1.;
  } else {
  sol = -1.;
  }
  return sol;
  });

  auto mu_user_func_solution =
    std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
  const double xx = x[0];
  const double yy = x[1];
  const double zz = x[2];
  const double r1 = (xx - M_PI + 1) * (xx - M_PI + 1) + (yy - M_PI) * (yy - M_PI) +
  (zz - M_PI) * (zz - M_PI);
  const double r2 = (xx - M_PI - 1) * (xx - M_PI - 1) + (yy - M_PI) * (yy - M_PI) +
  (zz - M_PI) * (zz - M_PI);
  const double r3 = (xx - M_PI - 1) * (xx - M_PI - 1) + (yy - M_PI) * (yy - M_PI) +
  (zz - M_PI) * (zz - M_PI);
  double sol = 0.;
  if (r1 < 1 || r2 < 1 || r3 < 1) {
  sol = 0;
  } else {
  sol = 0;
  }
  return sol;
  });

  auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto mu_initial_condition = AnalyticalFunctions<DIM>(mu_user_func_solution);
  const std::string& var_name_1 = "phi";
  const std::string& var_name_2 = "mu";
  auto v1 = VAR(&spatial, bcs_phi, var_name_1, 2, phi_initial_condition);
  auto v2 = VAR(&spatial, bcs_mu, var_name_2, 2, mu_initial_condition);
  auto vars = VARS(v1, v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  std::string calculation_path = "CahnHilliard";
  const double threshold = 10.;
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
    {var_name_1, {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  Parameters p_pst;
  const int frequency = p.post_processing ? 10 : 10000;
  p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
      Parameter("calculation_path", calculation_path), 
      Parameter("frequency", frequency),
      Parameter("level_of_detail", level_of_detail),
      Parameter("integral_to_compute", map_threshold_integral),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));

  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial, &spatial};
  OPE oper(spatials, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", p.verbosity),
  Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));

  switch (p.tcase)
  {
    case 0: oper.overload_solver(HypreSolverType::HYPRE_GMRES, Parameters(Parameter("tol", 1.e-12), Parameter("kdim", 100.0), Parameter("print_level", -1.0), Parameter("iter_max", 5000)));
      oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU);
      Message("Use a HypreGMRES solver with a HypreILU preconditioner");
      break;
    case 1: oper.overload_solver(HypreSolverType::HYPRE_PCG, Parameters(Parameter("tol", 1.e-12), Parameter("kdim", 100.0), Parameter("print_level", -1.0), Parameter("iter_max", 5000)));
      oper.overload_preconditioner(HyprePreconditionerType::HYPRE_BOOMER_AMG);
      Message("Use a HyprePCG solver with a BoomerAMG preconditioner");
      break;
  }

  auto pst = PST(&spatial, p_pst);
  PB problem1(oper, vars, pst);

  // Coupling 1
  auto cc = Coupling("CahnHilliard Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = p.duration;
  const double dt = p.dt;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
      Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  // time.get_tree();
  time.solve();
  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
