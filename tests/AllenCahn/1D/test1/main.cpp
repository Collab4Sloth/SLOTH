/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 1D AllenCahn problem along a radius
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
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 1;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////
  using OPE = TransientOperator<FECollection, DIM>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string& mesh_type = "InlineLineWithSegments";  // type of mesh
  const int order_fe = 1;                                   // finite element order
  const int refinement_level = 0;  // number of levels of uniform refinement
  const std::tuple<int, double>& tuple_of_dimensions =
      std::make_tuple(30, 1.e-3);  // Number of elements and maximum length

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("left", 0, "Neumann"), Boundary("right", 1, "Neumann")};
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
  // ####################
  //     variables     //
  // ####################
  const auto& center_x = 0.;
  const auto& a_x = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = 5.e-4;

  auto user_func = std::function<double(const mfem::Vector&, double)>(
      [center_x, a_x, radius, thickness](const mfem::Vector& x, [[maybe_unused]] double time) {
        const auto xx = a_x * (x[0] - center_x);
        const auto r = xx;
        const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::from("HyperbolicTangent"), center_x, a_x, epsilon, radius);
  auto vars = VARS(
      VAR(&spatial, bcs, "phi1", Glossary::PhaseField, 2, initial_condition, analytical_solution));
  auto vars2 = VARS(
      VAR(&spatial, bcs, "phi2", Glossary::PhaseField, 2, initial_condition, analytical_solution));
  auto vars3 = VARS(
      VAR(&spatial, bcs, "phi3", Glossary::PhaseField, 2, initial_condition, analytical_solution));
  auto vars4 = VARS(
      VAR(&spatial, bcs, "phi4", Glossary::PhaseField, 2, initial_condition, analytical_solution));
  auto vars5 = VARS(
      VAR(&spatial, bcs, "phi5", Glossary::PhaseField, 2, initial_condition, analytical_solution));

  // ####################
  //     coefficients  //
  // ####################

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient double_well_exp(Glossary::FreeEnergy, Scheme::Explicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  // ####################
  //     operators     //
  // ####################
  std::string calculation_path = "Problem1";
  auto p_pst1 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));

  // Problem 1:
  OPE oper({&spatial}, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  Coefficients coef_pb1(double_well_imp, capillary, mobility, grad_energy);
  auto pst = PST(&spatial, p_pst1);
  PB problem1(oper, vars, {coef_pb1}, pst);

  // Problem 2:
  calculation_path = "Problem2";
  auto p_pst2 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  OPE oper2({&spatial}, {"AllenCahn"}, TimeScheme::EulerExplicit, "TimeDerivative");
  Coefficients coef_pb2(double_well_exp, capillary, mobility, grad_energy);
  auto pst2 = PST(&spatial, p_pst2);
  PB problem2(oper2, vars2, {coef_pb2}, pst2);

  // Problem 3:
  calculation_path = "Problem3";
  auto p_pst3 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  OPE oper3({&spatial}, {"AllenCahn"}, TimeScheme::SDIRK33, "TimeDerivative");
  auto pst3 = PST(&spatial, p_pst3);
  PB problem3(oper3, vars3, {coef_pb1}, pst3);

  // Problem 4:
  calculation_path = "Problem4";
  auto p_pst4 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  OPE oper4({&spatial}, {"AllenCahn"}, TimeScheme::SDIRK23, "TimeDerivative");
  auto pst4 = PST(&spatial, p_pst4);
  PB problem4(oper4, vars4, {coef_pb1}, pst4);

  // Problem 5:
  calculation_path = "Problem5";
  auto p_pst5 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  OPE oper5({&spatial}, {"AllenCahn"}, TimeScheme::RungeKutta4, "TimeDerivative");
  auto pst5 = PST(&spatial, p_pst5);
  PB problem5(oper5, vars5, {coef_pb1}, pst5);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", problem1, problem2, problem3, problem4, problem5);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 50.0;
  const auto& dt = 0.01;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

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
