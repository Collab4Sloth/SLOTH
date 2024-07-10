/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a square with two periodic surface
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

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Couplings/Coupling.hpp"
#include "Integrators/AllenCahnNLFormIntegrator.hpp"
#include "Operators/PhaseFieldOperator.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Profiling/Profiling.hpp"
#include "Spatial/Spatial.hpp"
#include "Time/Time.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI
  //---------------------------------------
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  const auto DIM = 2;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI>;
  using PB = Problem<OPE, VAR, PST>;
  //###########################################
  //###########################################
  //        Spatial Discretization           //
  //###########################################
  //###########################################
  //##############################
  //          Meshing           //
  //##############################

  auto NN = 50;
  auto L = 2.e-3;
  // Create translation vectors defining the periodicity
  mfem::Vector x_translation({L, 0.0});
  std::vector<mfem::Vector> translations = {x_translation};
  auto refinement_level = 0;
  SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithQuadrangles", 1,
                                                   refinement_level, std::make_tuple(NN, NN, L, L),
                                                   translations);

  //##############################
  //    Boundary conditions     //
  //##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Periodic", 0.),
                     Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Periodic", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  //###########################################
  //###########################################
  //           Physical models               //
  //###########################################
  //###########################################
  //####################
  //    parameters    //
  //####################
  const auto& epsilon(3.e-4);
  // Two-phase mobility
  const auto& mob(5.e-5);
  const auto& sigma = 6.e-2;
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("mobility", mob), Parameter("lambda", lambda),
                           Parameter("omega", omega));
  //####################
  //    variables     //
  //####################
  const auto& center_x = 0.;
  const auto& center_y = 0.;
  const auto& a_x = 1.;
  const auto& a_y = 0.;
  const auto& thickness = 1.e-10;
  const auto& radius = 1.e-3;
  auto initial_condition = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, thickness, radius);
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y, a_x, a_y, epsilon, radius);

  auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, initial_condition));

  //###########################################
  //###########################################
  //     Post-processing                     //
  //###########################################
  //###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, params, vars);
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(main_folder_path, "Problem1", &spatial, frequency, level_of_detail);
  PB problem1("Problem 1", oper, vars, pst, TimeScheme::EulerImplicit, convergence, params);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", std::move(problem1));

  //###########################################
  //###########################################
  //           Time-integration              //
  //###########################################
  //###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.4;
  const auto& dt = 0.1;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, std::move(cc));

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
