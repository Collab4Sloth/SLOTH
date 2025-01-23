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
  constexpr int DIM = Test<2>::dim;
  using FECollection = Test<2>::FECollection;
  using VARS = Test<2>::VARS;
  using VAR = Test<2>::VAR;
  using PSTCollection = Test<2>::PSTCollection;
  using PST = Test<2>::PST;
  using SPA = Test<2>::SPA;
  using BCS = Test<2>::BCS;
  /////////////////////////
  using NLFI =
      DiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Implicit, Diffusion::Linear>;

  using OPE = DiffusionOperator<FECollection, DIM, NLFI, Density::Constant>;
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
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  const auto& alpha(1.e-2);
  const auto& kappa(0.5);
  // ####################
  //     variables     //
  // ####################
  auto user_func =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        if (x.Norml2() < 0.5) {
          return 1.0;
        } else {
          return 0.0;
        }
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);

  auto vars = VARS(VAR(&spatial, bcs, "c", 2, initial_condition));

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
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, TimeScheme::EulerImplicit);
  oper.overload_diffusion(Parameters(Parameter("D_0", kappa), Parameter("D_1", alpha)));

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst);

  PB problem1("Problem 1", oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.1;  // 0.5;
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
