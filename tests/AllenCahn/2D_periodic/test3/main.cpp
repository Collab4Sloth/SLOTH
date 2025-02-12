/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief Allen-Cahn problem solved in a square using MMS
 * @version 0.1
 * @date 2024-07-30
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
  // ---------------------------------------
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

  using OPE = SteadyAllenCahnOperator<FECollection, DIM, NLFI>;

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
  auto L = M_PI;
  std::vector<int> vect_order{1, 2};
  std::vector<int> vect_NN{8, 16, 32};

  for (const auto& order : vect_order) {
    for (const auto& NN : vect_NN) {
      //---------------------------------------
      // Profiling start
      Profiling::getInstance().enable();
      //---------------------------------------
      mfem::Vector x_translation({L, 0.0});
      mfem::Vector y_translation({0.0, L});
      std::vector<mfem::Vector> translations = {x_translation, y_translation};
      SPA spatial("InlineSquareWithQuadrangles", order, refinement_level,
                  std::make_tuple(NN, NN, L, L), translations);
      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                         Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
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
      const auto& epsilon(1.);
      // Interfacial energy
      const auto& sigma(1.);
      // Two-phase mobility
      const auto& mob(1.);
      const auto& lambda = 1.;
      const auto& omega = 1.;
      auto params = Parameters(Parameter("epsilon", epsilon), Parameter("mobility", mob),
                               Parameter("sigma", sigma), Parameter("lambda", lambda),
                               Parameter("omega", omega));
      // ####################,
      //     variables     //
      // ####################

      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double x = v[0];
            const double y = v[1];
            const auto func = std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2);
            return func;
          });
      auto user_func_source_term = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double x = v[0];
            const double y = v[1];
            const double H = 1.;
            const double a = 1.;
            const auto func =
                2 * H * (-2 * std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) + 1) *
                    (-std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) + 1) *
                    std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) -
                a * (18 * std::pow(std::sin(2 * x), 2) * std::pow(std::sin(3 * y), 2) -
                     26 * std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) +
                     8 * std::pow(std::cos(2 * x), 2) * std::pow(std::cos(3 * y), 2));
            return func;
          });

      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);

      auto vars = VARS(VAR(&spatial, bcs, "phi", 2, initial_condition, analytical_solution));

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
                              Parameter("level_of_detail", level_of_detail),
                              Parameter("force_clean_output_dir", false));
      // ####################
      //     operators     //
      // ####################

      // Problem 1:
      const auto crit_cvg_1 = 1.e-12;
      auto src_term = AnalyticalFunctions<DIM>(user_func_source_term);
      OPE oper(&spatial, params, src_term);
      oper.overload_mobility(Parameters(Parameter("mob", mob)));
      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      auto pst = PST(&spatial, p_pst);
      PB problem1("PSteady AllenCahn", oper, vars, pst, convergence);

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
    }
  }

  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
