/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a square with the splitting version of the equations
 * @version 0.1
 * @date 2025-07-07
 *
 * @copyright Copyright (c) 2025
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
  using LHS_NLFI = TimeCHNLFormIntegrator<VARS>;

  using NLFI = BlockAllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                              ThermodynamicsPotentials::W, Mobility::Constant>;
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
  const std::tuple<int, int, double, double>& tuple_of_dimensions = std::make_tuple(
      30, 30, 1.e-3, 1.e-3);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 1.),
                     Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Dirichlet", 0.)};
  auto boundaries_eta = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 0.),
                         Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Dirichlet", 0.)};
  auto bcs = BCS(&spatial, boundaries);
  auto bcs_eta = BCS(&spatial, boundaries_eta);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Interface thickness
  const double epsilon(5.e-4);
  // Interfacial energy
  const double sigma(6.e-2);
  // Two-phase mobility
  const double mob(1.e-5);
  const double lambda = 3. * sigma * epsilon / 2.;
  const double omega = 12. * sigma / epsilon;
  auto params =
      Parameters(Parameter("epsilon", epsilon), Parameter("epsilon", epsilon),
                 Parameter("sigma", sigma), Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################
  const double center_x = 0.;
  const double center_y = 0.;
  const double a_x = 1.;
  const double a_y = 0.;
  const double thickness = 5.e-5;
  const double radius = 5.e-4;

  auto initial_condition = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, thickness, radius);
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, epsilon, radius);
  auto v2 = VAR(&spatial, bcs_eta, "eta", 2, 0.);
  auto v1 = VAR(&spatial, bcs, "phi", 2, initial_condition, analytical_solution);
  auto vars = VARS(v1, v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const int frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial, &spatial};
  OPE oper(spatials, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", "Newton solver "),
                                     Parameter("print_level", 1), Parameter("abs_tol", 1.e-12)));
  oper.overload_solver(HypreSolverType::HYPRE_GMRES);
  oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU);
  auto pst = PST(&spatial, p_pst);
  PB problem1(oper, vars, pst);

  // Coupling 1
  auto cc = Coupling("AllenCahn-Splitting Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 1.;
  const double dt = 0.25;
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
