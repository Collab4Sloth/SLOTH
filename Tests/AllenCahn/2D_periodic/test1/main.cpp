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

/*!
 * \mainpage
 *
 * \tableofcontents
 * \section __quick Quick started
 * \subsection _main_sec0 Installation
 * \subsubsection __main_sub0 MFEM
 * \subsubsection __main_sub1 Post-processing tools
 *
 * \subsection _main_sec1 Development of a phase-field application?
 * A phase-field application is roughly a C++ main program composed of four parts:
 * \arg \ref __spatial "Spatial discretization"
 * \arg \ref __physical "Physical modeling"
 * \arg \ref __postprocessing "Post-processing directives"
 * \arg \ref __time "Time integration"
 *
 *
 */
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
                                         ThermodynamicsPotentials::F, Mobility::Constant>;
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
  // SpatialDiscretization<FECollection, DIM> spatial("GMSH", 1, "Mesh-examples/periodic.msh");

  auto NN = 32;
  auto L = 2. * M_PI;
  SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithQuadrangles", 1,
                                                   std::make_tuple(NN, NN, L, L));
  // Create translation vectors defining the periodicity
  mfem::Vector x_translation({L, 0.0});
  mfem::Vector y_translation({0.0, L});
  std::vector<mfem::Vector> translations = {x_translation, y_translation};

  spatial.make_periodic_mesh(translations);
  //##############################
  //    Boundary conditions     //
  //##############################
  // 2D y
  //    |_x
  auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                     Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
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
  const auto& epsilon(0.3);
  // Two-phase mobility
  const auto& mob(1.);
  const auto& lambda = 1.;
  const auto& omega = 1. / (epsilon * epsilon);
  auto params = Parameters(Parameter("mobility", mob), Parameter("lambda", lambda),
                           Parameter("omega", omega));
  //####################
  //    variables     //
  //####################
  auto vars =
      VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", "Unconserved", "Sinusoide",
                                      std::make_tuple(1.), "Sinusoide", std::make_tuple(1.)));
  //####################
  //    operators     //
  //####################
  // OPE oper(&spatial, params, vars);
  OPE oper(&spatial, params, vars, "Sinusoide2", omega);

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
  const auto& t_final = 1.;
  const auto& dt = 0.005;
  auto time_params =
      Parameters(Parameter("initial_time", t_initial), Parameter("final_time", t_final),
                 Parameter("time_step", dt), Parameter("compute_error", true));
  auto time = TIME("EulerImplicit", oper, time_params, vars, pst);
  time.execute();
  return 0;
}
