/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a 2D periodic domain
 * @version 0.1
 * @date 2024-05-23
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <cmath>
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
  // Profiling
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
                                         ThermodynamicsPotentials::F, Mobility::Constant>;
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
  // SpatialDiscretization<FECollection, DIM> spatial("GMSH", 1, 1, "periodic.msh", true);

  std::vector<int> vect_NN{16};  // 16, 32, 64};
  std::vector<std::string> vect_TimeScheme{"EulerImplicit", "EulerExplicit"};

  auto refinement_level = 0;
  for (const auto& time_scheme : vect_TimeScheme) {
    for (const auto& NN : vect_NN) {
      auto L = 2. * M_PI;
      // Create translation vectors defining the periodicity
      mfem::Vector x_translation({L, 0.0});
      mfem::Vector y_translation({0.0, L});
      std::vector<mfem::Vector> translations = {x_translation, y_translation};
      SPA spatial("InlineSquareWithQuadrangles", 1, refinement_level, std::make_tuple(NN, NN, L, L),
                  translations);

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
      const auto& epsilon(0.3);
      const auto& mob(1.);
      const auto& lambda = 1.;
      const auto& omega = 1. / (epsilon * epsilon);
      auto params = Parameters(Parameter("lambda", lambda), Parameter("omega", omega));
      // ####################
      //     variables     //
      // ####################
      auto analytical_solution = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Sinusoide, 1.);

      auto vars = VARS(VAR(&spatial, bcs, "phi", 2, analytical_solution, analytical_solution));

      // ###########################################
      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // ###########################################
      const std::string& main_folder_path = "Saves";
      const auto& level_of_detail = 1;
      const auto& frequency = 1;
      std::string calculation_path = "Problem1_" + time_scheme;
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
      auto source_terme = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Sinusoide2, omega);
      OPE oper(&spatial, params, TimeScheme::from(time_scheme), source_terme);
      oper.overload_mobility(Parameters(Parameter("mob", mob)));

      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      auto pst = PST(&spatial, p_pst);
      PB problem1("AllenCahn", oper, vars, pst, convergence);

      auto user_func = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, double time) { return 0.; });

      auto initial_rank = AnalyticalFunctions<DIM>(user_func);
      auto vars1 = VARS(VAR(&spatial, bcs, "MPI rank", 2, initial_rank));
      calculation_path = "ProblemMPI_";
      auto p_pst2 = Parameters(Parameter("main_folder_path", main_folder_path),
                               Parameter("calculation_path", calculation_path),
                               Parameter("frequency", frequency),
                               Parameter("level_of_detail", level_of_detail));
      auto pst2 = PST(&spatial, p_pst2);
      PB1 problem2(vars1, pst2, convergence);
      // Coupling 1
      auto cc = Coupling("AllenCahn-MPI Coupling", problem2, problem1);

      // ###########################################
      // ###########################################
      //            Time-integration              //
      // ###########################################
      // ###########################################
      const auto& t_initial = 0.0;
      const auto& t_final = 0.5;
      const auto& dt = 1. / std::pow(static_cast<double>(NN), 2.);
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
