/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a 2D periodic domain
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
  const auto DIM = 2;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::F, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, PhaseFieldOperatorBase>;
  using PB = Problem<OPE, VAR, PST>;
  using PB1 = MPI_Problem<VAR, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  // SpatialDiscretization<FECollection, DIM> spatial("GMSH", 1, 1, "periodic.msh", true);

  std::vector<int> vect_NN{16};  // 16, 32, 64};
  std::vector<std::string> vect_TimeScheme{"EulerImplicit", "EulerExplicit"};

  auto refinement_level = 0;
  for (const auto& time_scheme : vect_TimeScheme) {
    for (const auto& NN : vect_NN) {
      auto L = 2. * M_PI;
      // Create translation vectors defining the periodicity
      mfem::Vector x_translation({L, 0.0});
      mfem::Vector y_translation({0.0, L});
      std::vector<mfem::Vector> translations = {x_translation, y_translation};
      SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithQuadrangles", 1,
                                                       refinement_level,
                                                       std::make_tuple(NN, NN, L, L), translations);

      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                         Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
      auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

      // ###########################################
      // ###########################################
      //            Physical models               //
      // ###########################################
      // ###########################################
      // ####################
      //     parameters    //
      // ####################
      const auto& epsilon(0.3);
      const auto& mob(1.);
      const auto& lambda = 1.;
      const auto& omega = 1. / (epsilon * epsilon);
      auto params = Parameters(Parameter("mobility", mob), Parameter("lambda", lambda),
                               Parameter("omega", omega));
      // ####################
      //     variables     //
      // ####################
      auto analytical_solution = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Sinusoide, 1.);

      auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, analytical_solution,
                                                  analytical_solution));

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
      auto source_terme = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Sinusoide2, omega);
      OPE oper(&spatial, params, TimeScheme::from(time_scheme), source_terme);

      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      auto pst =
          PST(main_folder_path, "Problem1_" + time_scheme, &spatial, frequency, level_of_detail);
      PB problem1("AllenCahn", oper, vars, pst, convergence);

      auto user_func = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, double time) { return 0.; });

      auto initial_rank = AnalyticalFunctions<DIM>(user_func);
      auto vars1 = VAR(Variable<FECollection, DIM>(&spatial, bcs, "MPI rank", 2, initial_rank));
      auto pst1 =
          PST(main_folder_path, "ProblemMPI_" + time_scheme, &spatial, frequency, level_of_detail);
      PB1 problem2("MPI", vars1, pst1, convergence);
      // Coupling 1
      auto cc = Coupling("AllenCahn-MPI Coupling", problem2, problem1);

      // ###########################################
      // ###########################################
      //            Time-integration              //
      // ###########################################
      // ###########################################
      const auto& t_initial = 0.0;
      const auto& t_final = 0.5;
      const auto& dt = 1. / std::pow(static_cast<double>(NN), 2.);
      auto time_params = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);

      // time.get_tree();
      time.solve();
      //---------------------------------------
      // Profiling stop
      //---------------------------------------
      Profiling::getInstance().print();
    }
  }
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
