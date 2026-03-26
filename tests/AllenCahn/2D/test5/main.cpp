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
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

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
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using OPE = SteadyOperator<FECollection, DIM>;

  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string& mesh_type = "InlineSquareWithQuadrangles";  // type of mesh
  const int order_fe = 1;                                        // finite element order
  const int refinement_level = 0;  // number of levels of uniform refinement
  const std::tuple<int, int, double, double>& tuple_of_dimensions = std::make_tuple(
      200, 200, 1.e-3, 1.e-3);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 1.),
                     Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
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
  const auto& epsilon(2.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-3);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficients coef_ac(double_well_imp, capillary, mobility, grad_energy);
  // ####################
  //     variables     //
  // ####################
  const auto& center_x = 0.;
  const auto& center_y = 0.;
  const auto& a_x = 1.;
  const auto& a_y = 0.;
  const auto& thickness = 1.e-6;
  const auto& radius = 5.e-4;

  auto initial_condition = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, thickness, radius);
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, epsilon, radius);

  auto vars = VARS(
      VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition, analytical_solution));

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
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"AllenCahn"});
  auto pst = PST(&spatial, p_pst);
  PB problem1("Steady AllenCahn", oper, vars, {coef_ac}, pst);

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
