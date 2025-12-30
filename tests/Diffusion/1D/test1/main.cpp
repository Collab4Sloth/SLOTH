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
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using OPE = DiffusionOperator<FECollection, DIM>;
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
  std::vector<int> vect_NN{10};

  for (const auto& NN : vect_NN) {
    SPA spatial("InlineLineWithSegments", 1, refinement_level, std::make_tuple(NN, L));
    // ##############################
    //     Boundary conditions     //
    // // ##############################
    auto boundaries = {Boundary("left", 0, "Dirichlet", 1.), Boundary("right", 1, "Neumann", 0.)};
    auto bcs = BCS(&spatial, boundaries);

    // ###########################################
    // ###########################################
    //            Physical models               //
    // ###########################################
    // ###########################################
    // ####################
    //     parameters    //
    // ####################
    const auto& diffusionCoeff(1e-8);
    // ####################
    //     variables     //
    // ####################

    auto user_func_analytical = std::function<double(const mfem::Vector&, double)>(
        [diffusionCoeff](const mfem::Vector& x, double time) {
          const auto xx = x[0];
          const auto L_c = std::sqrt(4 * diffusionCoeff * time);
          const auto func = (1 - std::erf((xx) / L_c));
          return func;
        });

    auto initial_condition = 0.;
    auto analytical_solution = AnalyticalFunctions<DIM>(user_func_analytical);

    auto vars = VARS(
        VAR(&spatial, bcs, "c", Glossary::MoleFraction, 2, initial_condition, analytical_solution));

    // ###########################################
    // ###########################################
    //      Post-processing                     //
    // ###########################################
    // ###########################################
    const std::string& main_folder_path = "Saves_Nx_" + std::to_string(NN);
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
    Coefficient D(Glossary::Diffusivity, diffusionCoeff);
    Coefficient FW(Glossary::FreeEnergy, Scheme::Implicit, Log());
    Coefficients CoeffDiffusion(D, FW);

    // Problem 1:
    std::vector<SPA*> spatials{&spatial};
    OPE oper(spatials, {"Fick"}, TimeScheme::EulerImplicit, "TimeDerivative");

    auto pst = PST(&spatial, p_pst);

    PB problem1("Problem 1", oper, vars, {CoeffDiffusion}, pst);

    // Coupling 1
    auto cc = Coupling("coupling 1 ", problem1);

    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const auto& t_initial = 0.0;
    const auto& t_final = 5.;  // 0.5;
    const auto& dt = 0.1;      // 0.5 * ((L/NN)*(L/NN))/ (2 * diffusionCoeff);
    auto time_params = Parameters(Parameter("initial_time", t_initial),
                                  Parameter("final_time", t_final), Parameter("time_step", dt));
    auto time = TimeDiscretization(time_params, cc);

    time.solve();
    //---------------------------------------
    // Profiling stop
    //---------------------------------------
    Profiling::getInstance().print();
  }
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
