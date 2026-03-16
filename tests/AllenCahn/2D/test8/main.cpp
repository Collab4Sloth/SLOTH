/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a square : assessment of hypre solvers
 * @version 0.1
 * @date 2024-12-18
 *
 * Copyright CEA (c) 2024
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
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
  auto refinement_level = 0;
  SPA spatial("InlineSquareWithQuadrangles", 1, refinement_level,
              std::make_tuple(30, 30, 1.e-3, 1.e-3));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 1.),
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
  //  Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-4);
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
  const auto& thickness = 5.e-5;
  const auto& radius = 5.e-4;

  auto initial_condition = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, thickness, radius);
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, epsilon, radius);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;

  std::vector<HypreSolverType> vect_solver{
      HypreSolverType::HYPRE_FGMRES, HypreSolverType::HYPRE_GMRES, HypreSolverType::HYPRE_PCG};
  std::vector<HyprePreconditionerType> vect_precond{HyprePreconditionerType::HYPRE_ILU,
                                                    HyprePreconditionerType::HYPRE_BOOMER_AMG,
                                                    HyprePreconditionerType::HYPRE_DIAG_SCALE};

  int is = 0;
  for (const auto& solver : vect_solver) {
    int ip = 0;
    for (const auto& precond : vect_precond) {
      auto vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition,
                           analytical_solution));
      std::string calculation_path = "Problem_" + std::to_string(is) + "_" + std::to_string(ip);
      auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                              Parameter("calculation_path", calculation_path),
                              Parameter("frequency", frequency),
                              Parameter("level_of_detail", level_of_detail));
      // ####################
      //     operators     //
      // ####################

      // Problem 1:
      std::vector<SPA*> spatials{&spatial};
      OPE oper(spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
      oper.overload_solver(solver);
      oper.overload_preconditioner(precond);

      auto pst = PST(&spatial, p_pst);
      PB problem1(oper, vars, {coef_ac}, pst);

      // Coupling
      auto cc = Coupling("AllenCahn with ", problem1);

      // ###########################################
      // ###########################################
      //            Time-integration              //
      // ###########################################
      // ###########################################
      const auto& t_initial = 0.0;
      const auto& t_final = 30.;
      const auto& dt = 0.25;
      auto time_params = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);

      time.solve();

      ip++;
    }
    is++;
  }
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
