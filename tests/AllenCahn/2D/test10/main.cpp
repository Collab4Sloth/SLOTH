/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a 2D pellet fragment with a constant enthalpy melting
 * @version 0.1
 * @date 2024-05-23
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
  auto refinement_level = 0;
  SPA spatial("GMSH", 1, refinement_level, "pellet2Dinclusion.msh", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("external", 2, "Neumann"),
                     Boundary("upper", 1, "Neumann")};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //    Melting factor
  const double alpha(-7.e3);
  // Interface thickness
  const double epsilon(5.e-4);
  // Interfacial energy
  const double sigma(6.e-2);
  // Two-phase mobility
  const double mob(1.e-5);
  const double lambda = 3. * sigma * epsilon / 2.;
  const double omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("melting_factor", alpha));
  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient interpolation(Glossary::InterpolationFunction, Scheme::Implicit, H());
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficients coef_ac(double_well_imp, capillary, mobility, interpolation, grad_energy);
  // ####################
  //     variables     //
  // ####################
  auto vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, 1., {"cluster"}));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string main_folder_path = "Saves";
  const int level_of_detail = 1;
  const int frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst1 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst = PST(&spatial, p_pst1);

  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"AllenCahn", "MeltingConstant"}, params, TimeScheme::EulerImplicit,
           "TimeDerivative");

  auto nl_params = Parameters(Parameter("description", "Newton Algorithm"),
                              Parameter("abs_tol", 1.e-20), Parameter("rel_tol", 1.e-20));

  oper.overload_nl_solver(NLSolverType::NEWTON, nl_params);

  PB problem1("AllenCahn", oper, vars, {coef_ac}, pst);

  // Coupling 1
  auto cc = Coupling("Default Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 10.;
  const double dt = 0.25;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  // time.get_tree();
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
