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
  constexpr int DIM = Test<2>::dim;
  using FECollection = Test<2>::FECollection;
  using VARS = Test<2>::VARS;
  using VAR = Test<2>::VAR;
  using PSTCollection = Test<2>::PSTCollection;
  using PST = Test<2>::PST;
  using SPA = Test<2>::SPA;
  using BCS = Test<2>::BCS;
  /////////////////////////

  using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
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
  auto params =
      Parameters(Parameter("epsilon", epsilon), Parameter("epsilon", epsilon),
                 Parameter("sigma", sigma), Parameter("lambda", lambda), Parameter("omega", omega));
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
      auto vars = VARS(VAR(&spatial, bcs, "phi", 2, initial_condition, analytical_solution));
      std::string calculation_path = "Problem_" + std::to_string(is) + "_" + std::to_string(ip);
      auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                              Parameter("calculation_path", calculation_path),
                              Parameter("frequency", frequency),
                              Parameter("level_of_detail", level_of_detail));
      // ####################
      //     operators     //
      // ####################

      // Problem 1:
      const auto crit_cvg_1 = 1.e-12;
      OPE oper(&spatial, params, TimeScheme::EulerImplicit);
      oper.overload_mobility(Parameters(Parameter("mob", mob)));
      oper.overload_solver(solver);
      oper.overload_preconditioner(precond);

      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      auto pst = PST(&spatial, p_pst);
      PB problem1(oper, vars, pst, convergence);

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
