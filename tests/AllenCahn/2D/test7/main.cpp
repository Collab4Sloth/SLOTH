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

  using OPE = PhaseFieldOperator<FECollection, DIM>;
  using OPE2 = DiffusionOperator<FECollection, DIM>;
  using PB = Problem<OPE, VARS, PST>;
  using PB2 = Problem<OPE2, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  int order = 1;
  int NN = 20;
  double L = 1.;
  SPA spatial("InlineSquareWithQuadrangles", order, refinement_level,
              std::make_tuple(NN, NN, L, L));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.),
                     Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Dirichlet", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  // Two-phase mobility
  const auto& mob(1.e-2);
  const auto& lambda = 1.;
  const auto& omega = 0.;
  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficients coef_ac(double_well_imp, capillary, mobility, grad_energy);
  // ####################
  //     variables     //
  // ####################
  auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& x, [[maybe_unused]] double time) {
        if (x[0] > 0.5) {
          return 1.0;
        } else {
          return 0.0;
        }
      });
  auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto ex_func_solution =
      std::function<double(const mfem::Vector&, double)>([&](const mfem::Vector& x, double time) {
        double Lc = std::sqrt(4. * mob * time);
        double exact = 0.5 + 0.5 * std::erf((x[0] - L * 0.5) / Lc);

        return exact;
      });
  auto exact_condition = AnalyticalFunctions<DIM>(ex_func_solution);

  auto vars =
      VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 1, initial_condition, exact_condition));

  auto vars2 =
      VARS(VAR(&spatial, bcs, "c", Glossary::PhaseField, 1, initial_condition, exact_condition));
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
  std::string calculation_path2 = "Problem2";
  auto p_pst2 = Parameters(
      Parameter("main_folder_path", main_folder_path),
      Parameter("calculation_path", calculation_path2), Parameter("frequency", frequency),
      Parameter("level_of_detail", level_of_detail), Parameter("enable_compute_energies", false));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"AllenCahn"}, TimeScheme::EulerExplicit, "TimeDerivative");

  auto pst = PST(&spatial, p_pst);
  PB problem1(oper, vars, {coef_ac}, pst);

  // Problem 2:

  Coefficient D(Glossary::Diffusivity, mob);
  Coefficients CoeffDiffusion(D);
  OPE2 oper2(spatials, {"Fick"}, TimeScheme::EulerExplicit, "TimeDerivative");
  auto pst2 = PST(&spatial, p_pst2);
  PB2 problem2(oper2, vars2, {CoeffDiffusion}, pst2);

  // Coupling 1
  auto cc = Coupling("AllenCahn + Diffusion", problem1, problem2);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 1.;
  double c = 0.5;
  const auto& dt = c * (1 / (static_cast<double>(NN * NN))) / (4. * mob * order * order);
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
