/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief 1D AllenCahn problem with adaptive time-stepping: BDF2 scheme
 * @version 0.1
 * @date 2026-10-30
 *
 * @copyright Copyright (c) 2026
 *
 * Same case as test5 with the BDF2 scheme. The Newton solver is limited to 3 iterations: the
 * first time-steps fail and are rejected (variables, history and previous time-step of BDF2
 * restored). The growth factor (2.5) is larger than the ratio allowed between two successive
 * time-steps with BDF2 (time_step_max_ratio, 2 by default): the increase is limited to 2.
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
  setVerbosity(Verbosity::Debug);

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  using namespace Sloth1D;
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
  const auto& mob(1.e-3);
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
      VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition, analytical_solution));

  // ####################
  //     coefficients  //
  // ####################

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficient explicit_time_A(Glossary::ExplicitTime_A, 1.0);
  Coefficient explicit_time_B(Glossary::ExplicitTime_B, 1.0);
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
  const std::string& calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));

  TransientOPE oper({&spatial}, {"AllenCahn"}, TimeScheme::BDF2, "TimeDerivative");
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", "Newton solver "),
                                     Parameter("print_level", 1), Parameter("rel_tol", 1.e-10),
                                     Parameter("abs_tol", 1.e-12), Parameter("iter_max", 3)));
  Coefficients coefs(double_well, capillary, mobility, grad_energy, explicit_time_A,
                     explicit_time_B);
  auto pst = PST(&spatial, p_pst);
  TransientPB problem1(oper, vars, {coefs}, pst);

  auto cc = Coupling("coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 100.;
  // Initial time-step: half of the time interval
  const auto& dt = t_final * 0.5;
  auto time_params =
      Parameters(Parameter("initial_time", t_initial), Parameter("final_time", t_final),
                 Parameter("time_step", dt), Parameter("adaptive_time_step", true),
                 Parameter("time_step_min", 1.e-6), Parameter("time_step_cutback_factor", 0.5),
                 Parameter("time_step_growth_factor", 2.5), Parameter("newton_iterations_low", 2),
                 Parameter("newton_iterations_high", 4), Parameter("max_time_step_rejections", 30));
  auto time = TimeDiscretization(time_params, cc);

  time.get_tree();
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
