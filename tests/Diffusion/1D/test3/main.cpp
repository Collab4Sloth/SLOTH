/**
 * @file main.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Comparaison analytical and numerical solution
 * @version 0.1
 * @date 2024-11-28
 *
 * Copyright CEA (c) 2024
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
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
  const int DIM = 1;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using NLFI = MassDiffusionFluxNLFormIntegrator<VARS>;

  using LHS_NLFI = TimeNLFormIntegrator<VARS>;
  using OPE = DiffusionOperator<FECollection, DIM, NLFI, Density::Constant, LHS_NLFI>;
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
  double L = 1e-3;
  std::vector<int> vect_order{1, 2};
  std::vector<int> vect_NN{20, 40, 80, 160, 320};
  for (const auto& order : vect_order) {
    for (const auto& NN : vect_NN) {
      SPA spatial("InlineLineWithSegments", order, refinement_level, std::make_tuple(NN, L));
      // ##############################
      //     Boundary conditions     //
      // // ##############################
      auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};
      auto bcs = BCS(&spatial, boundaries);

      // ###########################################
      // ###########################################
      //            Physical models               //
      // ###########################################
      // ###########################################
      // ####################
      //     parameters    //
      // ####################
      const auto& stabCoeff(1e-8);
      const auto& diffusionCoeff(0.);
      //  ####################
      //      variables     //
      //  ####################

      auto user_func = std::function<double(const mfem::Vector&, double)>(
          [L](const mfem::Vector& x, double time) {
            const auto xx = x[0];
            const auto epsilon = 1e-4;
            auto func = (0.5 + 0.3 * std::tanh((xx - L / 2) / epsilon));

            const auto noiseLevel = 0.;
            static std::random_device rd;
            static std::mt19937 gen(rd());
            std::normal_distribution<> d(0, noiseLevel);

            return func + d(gen);
          });

      auto user_func_analytical = std::function<double(const mfem::Vector&, double)>(
          [L, diffusionCoeff](const mfem::Vector& x, double time) {
            const auto xx = x[0];
            const auto epsilon = 1e-3;
            auto func = 0.5 * (1 + std::erf((xx - L / 2) / std::sqrt(4 * diffusionCoeff * time)));

            return func;
          });

      auto initial_condition = AnalyticalFunctions<DIM>(user_func);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_analytical);

      auto vars = VARS(VAR(&spatial, bcs, "c", 2, initial_condition, initial_condition));

      // ###########################################
      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // ###########################################
      const std::string& main_folder_path =
          "Saves_order_" + std::to_string(order) + "_Nx" + std::to_string(NN);
      const auto& level_of_detail = 1;
      const auto& frequency = 1;
      std::string calculation_path = "Problem1";
      auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                              Parameter("calculation_path", calculation_path),
                              Parameter("frequency", frequency),
                              Parameter("level_of_detail", level_of_detail));
      // ####################
      //     operators     //
      // ####################
      // Thermal diffusion Parameters
      auto td_parameters = Parameters(Parameter("last_component", "U"));
      auto fictitious_mobO = VAR(&spatial, bcs, "Mo", 2, diffusionCoeff);
      // Fictitious mobilities
      fictitious_mobO.set_additional_information("C1_MO2", "O", "mob");
      auto fictitious_mobU = VAR(&spatial, bcs, "Mu", 2, diffusionCoeff);
      fictitious_mobU.set_additional_information("C1_MO2", "U", "mob");
      auto mobo_var = VARS(fictitious_mobO);
      auto mobu_var = VARS(fictitious_mobU);
      // Fictitious chemical potentials and mobilities
      auto fictitious_muo = VAR(&spatial, bcs, "muO", 2, 0.);
      fictitious_muo.set_additional_information("O", "mu");
      auto mu_var = VARS(fictitious_muo);
      auto fictitious_muu = VAR(&spatial, bcs, "muU", 2, 0.);
      fictitious_muu.set_additional_information("U", "mu");
      auto muu_var = VARS(fictitious_muu);

      auto fictitious_Mo = VAR(&spatial, bcs, "Mo", 2, 1.);
      fictitious_Mo.set_additional_information("O", "inter_mob");
      auto fictitious_Mu = VAR(&spatial, bcs, "Mu", 2, 1.);
      fictitious_Mu.set_additional_information("U", "inter_mob");

      auto fictitious_mob = VARS(fictitious_Mo, fictitious_Mu);
      // Problem 1:
      const auto crit_cvg_1 = 1.e-12;
      std::vector<SPA*> spatials{&spatial};
      OPE oper(spatials, td_parameters, TimeScheme::EulerImplicit);
      oper.overload_diffusion(Parameters(Parameter("D", stabCoeff)));

      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      auto pst = PST(&spatial, p_pst);

      PB problem1(oper, vars, pst, convergence, mu_var, muu_var, mobo_var, mobu_var,
                  fictitious_mob);

      // Coupling 1
      auto cc = Coupling("coupling 1 ", problem1);

      // ###########################################
      // ###########################################
      //            Time-integration              //
      // ###########################################
      // ###########################################
      const auto& t_initial = 0.0;
      const auto& t_final = 0.2;
      const auto& dt = 1e-4;
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
