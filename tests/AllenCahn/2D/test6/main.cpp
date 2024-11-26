/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief MMS benchmark from PFHUB : https://pages.nist.gov/pfhub/benchmarks/benchmark7.ipynb/
 * @version 0.1
 * @date 2024-10-25
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
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
  const auto DIM = 2;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
  using PB = Problem<OPE, VAR, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################

  std::vector<int> vect_order{1, 2};
  std::vector<int> vect_NN{4, 8, 16};
  for (const auto& order : vect_order) {
    for (const auto& NN : vect_NN) {
      // ##############################
      //           Meshing           //
      // ##############################
      auto refinement_level = 0;
      auto L = 1.;
      mfem::Vector x_translation({L, 0.0});
      std::vector<mfem::Vector> translations = {x_translation};
      SpatialDiscretization<FECollection, DIM> spatial(
          "InlineSquareWithQuadrangles", order, refinement_level,
          std::make_tuple(2 * NN, NN, L, L / 2), translations);
      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {
          Boundary("lower", 0, "Dirichlet", 1.), Boundary("right", 1, "Periodic", 0.),
          Boundary("upper", 2, "Dirichlet", 0.), Boundary("left", 3, "Periodic", 0.)};
      auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

      // ###########################################
      // ###########################################
      //            Physical models               //
      // ###########################################
      // ###########################################
      // ####################
      //     parameters    //
      // ####################
      const auto& epsilon(1.);
      const auto& sigma(1.);
      const auto& mob(1.);
      const auto& lambda(0.0004);
      const auto& omega(1.);
      auto params = Parameters(Parameter("epsilon", epsilon), Parameter("mobility", mob),
                               Parameter("sigma", sigma), Parameter("lambda", lambda),
                               Parameter("omega", omega));
      // ####################
      //     Parameters     //
      // ####################
      const auto A_1 = 0.0075;
      const auto A_2 = 0.03;
      const auto B_1 = 8. * M_PI;
      const auto B_2 = 22. * M_PI;
      const auto C_2 = 0.0625 * M_PI;
      const auto kappa = lambda;
      // ####################
      //   Analytical funct  //
      // ####################

      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [A_1, A_2, B_1, B_2, C_2, kappa](const mfem::Vector& v, double time) {
            // Exact solution cf. PFHUB website
            const double x = v[0];
            const double y = v[1];
            const double t = time;
            const auto func = 0.5 - 0.5 * std::tanh((1.0 / 2.0) * M_SQRT2 *
                                                    (-A_1 * t * std::sin(B_1 * x) -
                                                     A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                                    std::sqrt(kappa));
            return func;
          });
      auto user_func_source_term = std::function<double(const mfem::Vector&, double)>(
          [A_1, A_2, B_1, B_2, C_2, kappa](const mfem::Vector& v, double time) {
            // Source term for MMS cf. PFHUB website
            const double x = v[0];
            const double y = v[1];
            const double t = time;
            const auto func =
                -kappa *
                    (0.5 *
                         (1 - std::pow(std::tanh((1.0 / 2.0) * M_SQRT2 *
                                                 (-A_1 * t * std::sin(B_1 * x) -
                                                  A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                                 std::sqrt(kappa)),
                                       2)) *
                         std::pow(-A_1 * B_1 * t * std::cos(B_1 * x) -
                                      A_2 * B_2 * std::cos(B_2 * x + C_2 * t),
                                  2) *
                         std::tanh((1.0 / 2.0) * M_SQRT2 *
                                   (-A_1 * t * std::sin(B_1 * x) -
                                    A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                   std::sqrt(kappa)) /
                         kappa +
                     0.5 *
                         (1 - std::pow(std::tanh((1.0 / 2.0) * M_SQRT2 *
                                                 (-A_1 * t * std::sin(B_1 * x) -
                                                  A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                                 std::sqrt(kappa)),
                                       2)) *
                         std::tanh((1.0 / 2.0) * M_SQRT2 *
                                   (-A_1 * t * std::sin(B_1 * x) -
                                    A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                   std::sqrt(kappa)) /
                         kappa -
                     0.25 * M_SQRT2 *
                         (1 - std::pow(std::tanh((1.0 / 2.0) * M_SQRT2 *
                                                 (-A_1 * t * std::sin(B_1 * x) -
                                                  A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                                 std::sqrt(kappa)),
                                       2)) *
                         (A_1 * std::pow(B_1, 2) * t * std::sin(B_1 * x) +
                          A_2 * std::pow(B_2, 2) * std::sin(B_2 * x + C_2 * t)) /
                         std::sqrt(kappa)) -
                0.5 *
                    (2.0 - 2.0 * std::tanh((1.0 / 2.0) * M_SQRT2 *
                                           (-A_1 * t * std::sin(B_1 * x) -
                                            A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                           std::sqrt(kappa))) *
                    (-0.5 * std::tanh((1.0 / 2.0) * M_SQRT2 *
                                      (-A_1 * t * std::sin(B_1 * x) -
                                       A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                      std::sqrt(kappa)) -
                     0.5) *
                    std::tanh((1.0 / 2.0) * M_SQRT2 *
                              (-A_1 * t * std::sin(B_1 * x) - A_2 * std::sin(B_2 * x + C_2 * t) +
                               y - 0.25) /
                              std::sqrt(kappa)) -
                0.25 * M_SQRT2 *
                    (1 - std::pow(std::tanh((1.0 / 2.0) * M_SQRT2 *
                                            (-A_1 * t * std::sin(B_1 * x) -
                                             A_2 * std::sin(B_2 * x + C_2 * t) + y - 0.25) /
                                            std::sqrt(kappa)),
                                  2)) *
                    (-A_1 * std::sin(B_1 * x) - A_2 * C_2 * std::cos(B_2 * x + C_2 * t)) /
                    std::sqrt(kappa);
            return func;
          });

      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);

      auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, initial_condition,
                                                  analytical_solution));

      // ###########################################
      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // ###########################################
      const auto& t_initial = 0.0;
      const auto& t_final = 1e-1;
      const auto dt = 0.01;

      const std::string& main_folder_path = "Saves_order_" + std::to_string(order) + "_Nx" +
                                            std::to_string(NN) + "_dt" + std::to_string(dt);
      const auto& level_of_detail = 1;
      int frequency = 1;
      std::string calculation_path = "Problem1";
      auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                              Parameter("calculation_path", calculation_path),
                              Parameter("frequency", frequency),
                              Parameter("level_of_detail", level_of_detail));

      // ####################
      //     operators     //
      // ####################

      // Problem 1:
      const auto crit_cvg_1 = 1.e-12;
      auto src_term = AnalyticalFunctions<DIM>(user_func_source_term);

      OPE oper(&spatial, params, TimeScheme::EulerImplicit, src_term);

      auto nl_params = Parameters(Parameter("description", "Newton Algorithm"),
                                  Parameter("iterative_mode", false));
      auto s_params =
          Parameters(Parameter("description", "MINRES solver "), Parameter("print_level", 1));
      auto p_params =
          Parameters(Parameter("description", "Jacobi preconditionner"), Parameter("type", 1));

      oper.overload_nl_solver(NLSolverType::NEWTON, nl_params);
      oper.overload_solver(HypreSolverType::HYPRE_GMRES, s_params,
                           HyprePreconditionerType::HYPRE_ILU, p_params);
      oper.overload_mass_solver(HypreSolverType::HYPRE_GMRES, s_params,
                                HyprePreconditionerType::HYPRE_ILU, p_params);

      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

      auto pst = PST(&spatial, p_pst);
      PB problem1(oper, vars, pst, convergence);

      // Coupling 1
      auto cc = Coupling("Default Coupling", problem1);

      // ###########################################
      // ###########################################
      //            Time-integration              //
      // ###########################################
      // ###########################################

      auto time_params = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);

      // time.get_tree();
      time.solve();
      //---------------------------------------
      // Profiling stop
      //---------------------------------------
      Profiling::getInstance().print();
    }
  }
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
