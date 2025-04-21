/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a cube
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

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
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

  using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string& mesh_type = "InlineSquareWithTetraedres";  // type of mesh
  const int order_fe = 1;                                       // finite element order
  const int refinement_level = 0;  // number of levels of uniform refinement
  const std::tuple<int, int, int, double, double, double>& tuple_of_dimensions = std::make_tuple(
      30, 30, 30, 1.e-3, 1.e-3, 1.e-3);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("rear", 0, "Neumann", 0.),    Boundary("lower", 1, "Neumann", 0.),
                     Boundary("right", 2, "Dirichlet", 1.), Boundary("upper", 3, "Neumann", 0.),
                     Boundary("left", 4, "Dirichlet", 0.),  Boundary("front", 5, "Neumann", 0.)};
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
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                           Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################
  const auto& center_x = 0.;
  const auto& center_y = 0.;
  const auto& center_z = 0.;
  const auto& a_x = 1.;
  const auto& a_y = 0.;
  const auto& a_z = 0.;
  const auto& thickness = 5.e-5;
  const auto& radius = 5.e-4;
  auto initial_condition =
      AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y,
                               center_z, a_x, a_y, a_z, thickness, radius);
  auto analytical_solution =
      AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y,
                               center_z, a_x, a_y, a_z, epsilon, radius);

  auto vars = VARS(VAR(&spatial, bcs, "phi", 2, initial_condition, analytical_solution));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst);
  PB problem1("AllenCahn", oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("Default Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.25;
  const auto& dt = 0.25;
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
