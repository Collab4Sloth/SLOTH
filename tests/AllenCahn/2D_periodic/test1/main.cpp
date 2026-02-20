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
  setVerbosity(Verbosity::Debug);
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
      Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
      Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, Fw(omega));
      Coefficient capillary(Glossary::Capillary, lambda);
      Coefficient mobility(Glossary::Mobility, mob);
      Coefficients coef_ac(double_well_imp, capillary, mobility, grad_energy);
      // ####################
      //     variables     //
      // ####################
      auto analytical_solution = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Sinusoide, 1.);

      auto vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, analytical_solution,
                           analytical_solution));

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
      std::vector<SPA*> spatials{&spatial};
      std::vector<AnalyticalFunctions<DIM> > src_term;
      src_term.emplace_back(AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Sinusoide2, omega));
      OPE oper(spatials, {"AllenCahn"}, TimeScheme::from(time_scheme), "TimeDerivative", src_term);

      auto pst = PST(&spatial, p_pst);
      PB problem1("AllenCahn", oper, vars, {coef_ac}, pst);

      // Coupling 1
      auto cc = Coupling("AllenCahn Coupling", problem1);

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
