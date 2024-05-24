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
#include "Coefficients/AnalyticalFunctions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Integrators/AllenCahnNLFormIntegrator.hpp"
#include "Operators/PhaseFieldOperator.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Spatial/Spatial.hpp"
#include "Time/Time.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  const auto DIM = 1;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI>;
  using TIME = TimeDiscretization<PST, OPE, VAR>;
  //###########################################
  //###########################################
  //        Spatial Discretization           //
  //###########################################
  //###########################################
  //##############################
  //          Meshing           //
  //##############################
  auto refinement_level = 0;
  SpatialDiscretization<FECollection, DIM> spatial("InlineLineWithSegments", 1, refinement_level,
                                                   std::make_tuple(30, 1.e-3));
  //##############################
  //    Boundary conditions     //
  //##############################
  auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  //###########################################
  //###########################################
  //           Physical models               //
  //###########################################
  //###########################################
  //####################
  //    parameters    //
  //####################
  // Interface thickness
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
  //####################
  //    variables     //
  //####################
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
  //####################
  //    operators     //
  //####################
  OPE oper(&spatial, params, vars);

  //###########################################
  //###########################################
  //     Post-processing                     //
  //###########################################
  //###########################################
  const std::string& main_folder_path = "Paraview";
  const std::string& calculation_path = "MainPST";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  auto pst = PST("Paraview", "MainPST", &spatial, frequency, level_of_detail);

  //###########################################
  //###########################################
  //           Time-integration              //
  //###########################################
  //###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 100.;
  const auto& dt = 0.25;
  auto time_params =
      Parameters(Parameter("initial_time", t_initial), Parameter("final_time", t_final),
                 Parameter("time_step", dt), Parameter("compute_error", true),
                 Parameter("compute_energies", true));
  auto time = TIME("EulerImplicit", oper, time_params, vars, pst);
  time.execute();
  return 0;
}
