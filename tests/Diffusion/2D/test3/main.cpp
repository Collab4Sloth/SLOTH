/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Non linear Diffusion problem solved in a square (non linear version of the test 16 in
 * mfem.org page)
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

#include "kernel/Coefficients/Coefficient.hpp"
#include "kernel/Coefficients/Coefficients.hpp"
#include "kernel/Coefficients/EnergyCoefficients.hpp"
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
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  //
  using OPE = SteadyDiffusionOperator<FECollection, DIM>;
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
  SPA spatial("InlineSquareWithQuadrangles", 1, refinement_level,
              std::make_tuple(200, 200, 1., 1.));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Dirichlet", 0.), Boundary("right", 1, "Dirichlet", 0.),
                     Boundary("upper", 2, "Dirichlet", 0.), Boundary("left", 3, "Dirichlet", 0.)};

  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     Coefficients  //
  // ####################
  std::vector<Coefficients> coeffs;
  Coefficient D(Glossary::Diffusivity, Scheme::Implicit, "-1");
  Coefficient FW(Glossary::FreeEnergy, Scheme::Implicit, energylog());
  Coefficients CoeffDiffusion(D, FW);
  coeffs.emplace_back(CoeffDiffusion);
  // ####################
  //     variables     //
  // ####################
  auto initial_condition = 0.0;

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

  auto user_func_source_term =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& v, double time) {
        const auto func = 1.;
        return func;
      });

  std::vector<AnalyticalFunctions<DIM> > src_term;
  src_term.emplace_back(AnalyticalFunctions<DIM>(user_func_source_term));

  // Problem 1:
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"Fick"}, src_term);

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
  const auto& t_final = 0.1;
  const auto& dt = 0.1;
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
