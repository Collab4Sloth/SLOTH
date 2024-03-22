/*
 * Copyright © CEA 2022
 *
 * \brief Main program for the PF-MFEM short application
 * \file main.cpp
 * \author ci230846
 * \date 11/01/2022
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/AnalyticalFunctions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Integrators/AllenCahnMeltingNLFormIntegrator.hpp"
#include "Operators/ConductionOperator.hpp"
#include "Operators/PhaseFieldOperatorMelting.hpp"
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
  const auto DIM = 3;
  using NLFI =
      AllenCahnMeltingNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                       ThermodynamicsPotentials::W, ThermodynamicsPotentials::H,
                                       Mobility::Constant, PhaseChange::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperatorMelting<FECollection, DIM, NLFI>;
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
  SpatialDiscretization<mfem::H1_FECollection, DIM> spatial(
      "GMSH", 1, refinement_level, "../../../../Mesh-examples/camembert3D.msh", false);
  //##############################
  //    Boundary conditions     //
  //##############################
  // 2D y
  //    |_x
  auto boundaries = {
      Boundary("InterPelletPlane", 1, "Neumann", 0.), Boundary("MidPelletPlane", 2, "Neumann", 0.),
      Boundary("FrontSurface", 3, "Neumann", 0.), Boundary("BehindSurface", 4, "Neumann", 0.),
      Boundary("ExternalSurface", 0, "Neumann", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  //###########################################
  //###########################################
  //           Physical models               //
  //###########################################
  //###########################################
  //####################
  //    parameters    //
  //####################
  // Melting factor
  const auto& alpha(7.e3);
  // Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params =
      Parameters(Parameter("epsilon", epsilon), Parameter("epsilon", epsilon),
                 Parameter("mobility", mob), Parameter("sigma", sigma), Parameter("lambda", lambda),
                 Parameter("omega", omega), Parameter("melting_factor", alpha));
  //####################
  //    variables     //
  //####################
  const auto& pellet_radius = 0.00465;
  const auto& pellet_height = 0.01;
  const auto& center_x = 0.;
  const auto& center_y = 0.;
  const auto& center_z = 0.5 * pellet_height;
  const auto& a_x = 1.;
  const auto& a_y = 1.;
  const auto& a_z = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = 1.e-1 * pellet_radius;
  auto vars = VAR(Variable<FECollection, DIM>(
      &spatial, bcs, "phi", "Unconserved", "HyperbolicTangent",
      std::make_tuple(center_x, center_y, center_z, a_x, a_y, a_z, thickness, radius),
      "HyperbolicTangent",
      std::make_tuple(center_x, center_y, center_z, a_x, a_y, a_z, epsilon, radius)));
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
