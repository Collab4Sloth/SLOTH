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
#include <string>
#include <vector>

#include "./Coefficient.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Verbose);

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  using namespace Sloth2D;

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
      SPA spatial("InlineSquareWithQuadrangles", order, refinement_level,
                  std::make_tuple(2 * NN, NN, L, L / 2), translations);
      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Dirichlet", 1.), Boundary("right", 1, "Periodic"),
                         Boundary("upper", 2, "Dirichlet", 0.), Boundary("left", 3, "Periodic")};
      auto bcs = BCS(&spatial, boundaries);

      // ###########################################
      // ###########################################
      //            Physical models               //
      // ###########################################
      // ###########################################
      // ####################
      //     parameters    //
      // ####################
      const auto& mob(1.);
      const auto& lambda(0.0004);
      const auto& omega(1.);
      Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
      Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
      Coefficient capillary(Glossary::Capillary, lambda);
      Coefficient mobility(Glossary::Mobility, mob);
      Coefficient source(Glossary::SourceTerm, Scheme::Explicit, SourceCoefficient());
      Coefficients coef_ac(double_well_imp, capillary, mobility, grad_energy, source);
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
      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);

      auto vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition,
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

      std::vector<SPA*> spatials{&spatial};
      TransientOPE oper(spatials, {"AllenCahn", "SourceTerm"}, TimeScheme::EulerImplicit,
                        "TimeDerivative");

      auto pst = PST(&spatial, p_pst);
      TransientPB problem1(oper, vars, {coef_ac}, pst);

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
