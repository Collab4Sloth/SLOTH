/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 1D AllenCahn problem along a radius
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

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Couplings/Coupling.hpp"
#include "Integrators/AllenCahnNLFormIntegrator.hpp"
#include "Operators/PhaseFieldOperator.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/Problem.hpp"
#include "Profiling/Profiling.hpp"
#include "Spatial/Spatial.hpp"
#include "Time/Time.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Utils/UtilsForDebug.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

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
  const auto DIM = 1;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using NLFI2 = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Explicit,
                                          ThermodynamicsPotentials::W, Mobility::Constant>;
  using NLFI3 = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::SemiImplicit,
                                          ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI>;
  using OPE2 = PhaseFieldOperator<FECollection, DIM, NLFI2>;
  using OPE3 = PhaseFieldOperator<FECollection, DIM, NLFI3>;
  using PB = Problem<OPE, VAR, PST>;
  using PB2 = Problem<OPE2, VAR, PST>;
  using PB3 = Problem<OPE3, VAR, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  SpatialDiscretization<FECollection, DIM> spatial("InlineLineWithSegments", 1, refinement_level,
                                                   std::make_tuple(30, 1.e-3));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("epsilon", epsilon),
                           Parameter("mobility", mob), Parameter("sigma", sigma),
                           Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################
  const auto& center_x = 0.;
  const auto& a_x = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = 5.e-4;

  auto user_func = std::function<double(const mfem::Vector&, double)>(
      [center_x, a_x, radius, thickness](const mfem::Vector& x, double time) {
        const auto xx = a_x * (x[0] - center_x);
        const auto r = xx;
        const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);
  auto analytical_solution = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent,
                                                      center_x, a_x, epsilon, radius);
  auto vars = VAR(
      Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, initial_condition, analytical_solution));
  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
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

  // Problem 2:
  const auto crit_cvg_2 = 1.e-12;
  OPE2 oper2(&spatial, params, vars);
  PhysicalConvergence convergence2(ConvergenceType::RELATIVE_MAX, crit_cvg_2);
  auto pst2 = PST(main_folder_path, "Problem2", &spatial, frequency, level_of_detail);
  PB2 problem2("Problem 2 ", oper2, vars, pst2, TimeScheme::EulerExplicit, convergence2, params);

  // Problem 3:
  const auto crit_cvg_3 = 1.e-12;
  OPE3 oper3(&spatial, params, vars);
  PhysicalConvergence convergence3(ConvergenceType::RELATIVE_MAX, crit_cvg_3);
  auto pst3 = PST(main_folder_path, "Problem3", &spatial, frequency, level_of_detail);
  PB3 problem3("Problem 3 ", oper3, vars, pst3, TimeScheme::RungeKutta4, convergence3, params);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", std::move(problem1), std::move(problem2), std::move(problem3));

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 50.0;
  const auto& dt = 0.01;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, std::move(cc));

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
