/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief 2D spinodal decomposition solved by Cahn-Hilliard equations
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
                                            ThermodynamicsPotentials::WW, Mobility::Constant>;
  using LHS_NLFI = TimeCHNLFormIntegrator<VARS>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, LHS_NLFI>;
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
  const int nx = 100;
  const int ny = 100;
  const double lx = 200;
  const double ly = 200;
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
  // Interfacial energy
  const double sigma(1.);
  // Two-phase mobility
  const double mob(5.);
  const double lambda = 2.;
  const double omega = 5.;
  auto params =
      Parameters(Parameter("sigma", sigma), Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################

  auto user_func_solution =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double co = 0.5;
        double epsilon = 0.01;
        double xx = x[0];
        double yy = x[1];

        double sol =
            co + epsilon * (std::cos(0.105 * xx) * std::cos(0.11 * yy) +
                            (std::cos(0.13 * xx) * std::cos(0.087 * yy)) *
                                (std::cos(0.13 * xx) * std::cos(0.087 * yy)) +
                            (std::cos(0.025 * xx - 0.15 * yy) * std::cos(0.07 * xx - 0.02 * yy)));

        return sol;
      });

  auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto mu_initial_condition = 0.0;
  const std::string& var_name_1 = "phi";
  const std::string& var_name_2 = "mu";
  auto v1 = VAR(&spatial, bcs, var_name_1, 2, phi_initial_condition);
  auto v2 = VAR(&spatial, bcs, var_name_2, 2, mu_initial_condition);
  auto vars = VARS(v1, v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const int frequency = 100;
  std::string calculation_path = "CahnHilliard";
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
      {var_name_1, {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
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
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", "Newton solver "),
                                     Parameter("print_level", -1), Parameter("rel_tol", 1.e-12),
                                     Parameter("abs_tol", 1.e-12), Parameter("iter_max", 1000)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_ILU;
  oper.overload_solver(solver);
  oper.overload_preconditioner(precond);
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
  const double t_final = 10.;  // 5.e4;
  const double dt = 1.;
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
