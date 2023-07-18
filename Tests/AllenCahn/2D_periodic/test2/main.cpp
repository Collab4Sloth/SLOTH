/*
 * Copyright © CEA 2022
 *
 * \brief Main program for the PF-MFEM short application
 * \file main.cpp
 * \author ci230846
 * \date 11/01/2022
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/AnalyticalFunctions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Integrators/AllenCahnNLFormIntegrator.hpp"
#include "Operators/ConductionOperator.hpp"
#include "Operators/PhaseFieldOperator.hpp"
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
  // int DIM1 = 1;
  // mfem::OptionsParser args(argc, argv);
  // args.AddOption(&DIM1, "-d", "--dimension", "dimension");
  // args.ParseCheck();
  // const int DIM = DIM1;
  const auto DIM = 2;
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
  // SpatialDiscretization<FECollection, DIM> spatial(
  //     "GMSH", 1, "../../../../Mesh-examples/periodic_002.msh", true);
  //##############################
  //          Meshing           //
  //##############################

  auto NN = 50;
  auto L = 2.e-3;
  // Create translation vectors defining the periodicity
  mfem::Vector x_translation({L, 0.0});
  // mfem::Vector y_translation({0.0, L});
  std::vector<mfem::Vector> translations = {x_translation};
  SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithQuadrangles", 1,
                                                   std::make_tuple(NN, NN, L, L), translations);

  // // spatial.make_periodic_mesh(translations);
  //##############################
  //    Boundary conditions     //
  //##############################
  // 2D y
  //    |_x
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Periodic", 0.),
                     Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Periodic", 0.)};
  // auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("upper", 2, "Neumann", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  //###########################################
  //###########################################
  //           Physical models               //
  //###########################################
  //###########################################
  //####################
  //    parameters    //
  //####################
  // Cahn number
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
  auto vars = VAR(Variable<FECollection, DIM>(
      &spatial, bcs, "phi", "Unconserved", "HyperbolicTangent",
      std::make_tuple(center_x, center_y, a_x, a_y, thickness, radius), "HyperbolicTangent",
      std::make_tuple(center_x, center_y, a_x, a_y, epsilon, radius)));
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
  const std::string& calculation_path = "EulerImplicit";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  auto pst = PST(main_folder_path, calculation_path, &spatial, frequency, level_of_detail);

  //###########################################
  //###########################################
  //           Time-integration              //
  //###########################################
  //###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 50.;
  const auto& dt = 0.1;
  auto time_params =
      Parameters(Parameter("initial_time", t_initial), Parameter("final_time", t_final),
                 Parameter("time_step", dt), Parameter("compute_error", true),
                 Parameter("compute_energies", true));
  auto time = TIME("EulerImplicit", oper, time_params, vars, pst);
  time.execute();
  return 0;
}
