/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a 2D pellet fragment with a constant enthalpy melting
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
#include <vector>

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
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using NLFI = AllenCahnConstantMeltingNLFormIntegrator<
      VARS, ThermodynamicsPotentialDiscretization::Implicit, ThermodynamicsPotentials::W,
      Mobility::Constant, ThermodynamicsPotentials::H>;

  using LHS_NLFI = TimeNLFormIntegrator<VARS>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, LHS_NLFI>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  SPA spatial("GMSH", 1, refinement_level, "camembert2D.", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Neumann", 0.),
                     Boundary("upper", 1, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //    Melting factor
  const auto& alpha(7.e3);
  // Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("epsilon", epsilon),
                           Parameter("sigma", sigma), Parameter("lambda", lambda),
                           Parameter("omega", omega), Parameter("melting_factor", alpha));
  // ####################
  //     variables     //
  // ####################
  const auto& pellet_radius = 0.00465;
  const auto& center_x = 0.;
  const auto& center_y = 0.;
  const auto& a_x = 1.;
  const auto& a_y = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = 1.e-2 * pellet_radius;

  auto initial_condition = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, thickness, radius);

  auto vars = VARS(VAR(&spatial, bcs, "phi", 2, initial_condition));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst1 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst = PST(&spatial, p_pst1);

  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));

  PB problem1("AllenCahn", oper, vars, pst);

  // Coupling 1
  auto cc = Coupling("Default Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 1.;
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
