/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in 3D piece of pellet fragment
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

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------

  mfem::Mpi::Init(argc, argv);
  int size = mfem::Mpi::WorldSize();
  int rank = mfem::Mpi::WorldRank();
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  const auto DIM = 3;
  using NLFI = AllenCahnMeltingNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                                ThermodynamicsPotentials::W, Mobility::Constant,
                                                ThermodynamicsPotentials::H, PhaseChange::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, PhaseFieldOperatorBase>;
  using PB = Problem<OPE, VAR, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  SpatialDiscretization<mfem::H1_FECollection, DIM> spatial("GMSH", 1, refinement_level,
                                                            "camembert3D.msh", false);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {
      Boundary("InterPelletPlane", 1, "Neumann", 0.), Boundary("MidPelletPlane", 2, "Neumann", 0.),
      Boundary("FrontSurface", 3, "Neumann", 0.), Boundary("BehindSurface", 4, "Neumann", 0.),
      Boundary("ExternalSurface", 0, "Neumann", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Melting factor
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
  // ####################
  //     variables     //
  // ####################
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

  auto initial_condition =
      AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y,
                               center_z, a_x, a_y, a_z, thickness, radius);

  auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, initial_condition));

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
  OPE oper(&spatial, params, TimeScheme::EulerImplicit);

  auto nl_params = Parameters(Parameter("description", "Newton Algorithm"),
                              Parameter("abs_tol", 1.e-20), Parameter("rel_tol", 1.e-20));

  oper.overload_nl_solver(NLSolverType::NEWTON, nl_params);

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(main_folder_path, "Problem1", &spatial, frequency, level_of_detail);
  PB problem1("AllenCahn", oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("Default Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.25;
  const auto& dt = 0.25;
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
