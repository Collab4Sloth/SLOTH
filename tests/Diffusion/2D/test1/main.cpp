/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Diffusion problem solved in a square (similar to the test 16 in mfem.org page)
 * @version 0.1
 * @date 2024-06-06
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
  //
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
  auto refinement_level = 2;
  SPA spatial("GMSH", 2, refinement_level, "star2D.msh", false);
  // ##############################
  //     Boundary conditions     //
  // // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann")};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     Coefficients  //
  // ####################

  std::string function = "0.5 + 0.01*c";
  std::string gx = "0.01";
  std::vector<std::string> grad_functions{gx};
  std::vector<Coefficients> coeffs;
  Coefficient D(Glossary::Diffusivity, Scheme::Explicit, grad_functions, function, "c");
  Coefficient FW(Glossary::FreeEnergy, Scheme::Implicit, Log());
  Coefficients CoeffDiffusion(D, FW);
  coeffs.emplace_back(CoeffDiffusion);
  // ####################
  //     variables     //
  // ####################
  auto user_func = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& x, [[maybe_unused]] double time) {
        if (x.Norml2() < 0.5) {
          return 1.0;
        } else {
          return 0.0;
        }
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);

  auto vars = VARS(VAR(&spatial, bcs, "c", Glossary::MoleFraction, 2, initial_condition));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"Fick"}, TimeScheme::EulerImplicit, "TimeDerivative");

  auto pst = PST(&spatial, p_pst);

  PB problem1("Problem 1", oper, vars, coeffs, pst);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.5;
  const auto& dt = 0.01;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  time.solve();
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
