/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a square
 * @version 0.1
 * @date 2024-05-23
 *
 * @copyright Copyright (c) 2024
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
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Debug);

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using NLFI = CahnHilliardNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                            ThermodynamicsPotentials::W, Mobility::Constant>;

  using LHS_NLFI = TimeCHNLFormIntegrator<VARS>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI, LHS_NLFI>;
  using PB = Problem<OPE, VARS, PST>;
  using PB1 = MPI_Problem<VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string mesh_type =
      "InlineSquareWithQuadrangles";  // type of mesh // "InlineSquareWithTriangles"
  const int order_fe = 1;             // finite element order
  const int refinement_level = 0;     // number of levels of uniform refinement
  const int nx = 256;
  const int ny = 256;
  const double lx = 2. * M_PI;
  const double ly = 2. * M_PI;
  const std::tuple<int, int, double, double>& tuple_of_dimensions =
      std::make_tuple(nx, ny, lx, ly);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.),
                     Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Interface thickness
  const double epsilon(0.1);
  // Interfacial energy
  const double sigma(1.);
  // Two-phase mobility
  const double mob(1.);
  const double lambda = epsilon * epsilon;
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
        const double r1 = (xx - M_PI + 1) * (xx - M_PI + 1) + (yy - M_PI) * (yy - M_PI);
        const double r2 = (xx - M_PI - 1) * (xx - M_PI - 1) + (yy - M_PI) * (yy - M_PI);
        double sol = 0.;
        if (r1 < 1 || r2 < 1) {
          sol = 1.;
        } else {
          sol = 0.;
        }
        return sol;
      });

  auto mu_user_func_solution =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        const double xx = x[0];
        const double yy = x[1];
        const double r1 = (xx - M_PI + 1) * (xx - M_PI + 1) + (yy - M_PI) * (yy - M_PI);
        const double r2 = (xx - M_PI - 1) * (xx - M_PI - 1) + (yy - M_PI) * (yy - M_PI);
        double sol = 0.;
        if (r1 < 1 || r2 < 1) {
          sol = 0;
        } else {
          sol = 0;
        }
        return sol;
      });

  auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto mu_initial_condition = AnalyticalFunctions<DIM>(mu_user_func_solution);
  auto v1 = VAR(&spatial, bcs, "phi", 2, phi_initial_condition);
  auto v2 = VAR(&spatial, bcs, "mu", 2, mu_initial_condition);
  auto vars = VARS(v1, v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const int frequency = 50;
  std::string calculation_path = "CahnHilliard";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  std::vector<SPA*> spatials{&spatial, &spatial};
  OPE oper(spatials, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_ILU;
  oper.overload_solver(solver);
  oper.overload_preconditioner(precond);

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst);
  PB problem1(oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("CahnHilliard Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 50.;
  const double dt = 1.e-3;
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
