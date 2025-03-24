/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a square : linear manufactured solution
 * @version 0.1
 * @date 2024-07-19
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
  /////////////////////////
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////
  using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;

  using OPE = SteadyAllenCahnOperator<FECollection, DIM, NLFI>;

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
  auto L = 1.;
  std::vector<int> vect_order{1};
  std::vector<int> vect_NN{4};
  for (const auto& order : vect_order) {
    for (const auto& NN : vect_NN) {
      //---------------------------------------
      // Profiling start
      Profiling::getInstance().enable();
      //---------------------------------------
      SPA spatial("InlineSquareWithQuadrangles", order, refinement_level,
                  std::make_tuple(NN, 5, L, L));

      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 2.),
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
      const auto& epsilon(1.);
      const auto& sigma(1.);
      const auto& mob(1.);
      const auto& lambda(1.);
      const auto& omega(1.);
      auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                               Parameter("lambda", lambda), Parameter("omega", omega));
      // ####################
      //     variables     //
      // ####################
      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double x = v[0];
            const double y = v[1];
            const auto func = 2. * x;
            return func;
          });
      auto user_func_source_term = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double x = v[0];
            const double y = v[1];
            const auto func = 4 * x * (1 - 4 * x) * (1 - 2 * x);
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
      auto p_pst1 = Parameters(Parameter("main_folder_path", main_folder_path),
                               Parameter("calculation_path", calculation_path),
                               Parameter("frequency", frequency),
                               Parameter("level_of_detail", level_of_detail));
      auto pst = PST(&spatial, p_pst1);
      // ####################
      //     operators     //
      // ####################

      // Problem 1:
      const auto crit_cvg_1 = 1.e-12;
      auto src_term = AnalyticalFunctions<DIM>(user_func_source_term);
      OPE oper(&spatial, params, src_term);
      oper.overload_mobility(Parameters(Parameter("mob", mob)));
      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      PB problem1("Steady AllenCahn", oper, vars, pst, convergence);

      auto vars1 = VARS(VAR(&spatial, bcs, "MPI rank", 2, 0.));

      calculation_path = "ProblemMPI_";
      auto p_pst2 = Parameters(Parameter("main_folder_path", main_folder_path),
                               Parameter("calculation_path", calculation_path),
                               Parameter("frequency", frequency),
                               Parameter("level_of_detail", level_of_detail));
      auto pst2 = PST(&spatial, p_pst2);
      PB1 problem2(vars1, pst2, convergence);
      // Coupling 1
      auto cc = Coupling("Steady AllenCahn-MPI Coupling", problem2, problem1);

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
