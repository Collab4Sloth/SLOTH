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

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Couplings/Coupling.hpp"
#include "Integrators/DiffusionNLFormIntegrator.hpp"
#include "Operators/DiffusionOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Spatial/Spatial.hpp"
#include "Time/Time.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
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
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  const auto DIM = 2;
  using NLFI = DiffusionNLFormIntegrator<DiffusionCoefficientDiscretization::Implicit,
                                         DiffusionCoefficients::Linear>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = DiffusionOperator<FECollection, DIM, NLFI, SteadyPhaseFieldOperatorBase>;
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
  SpatialDiscretization<mfem::H1_FECollection, DIM> spatial(
      "InlineSquareWithQuadrangles", 1, refinement_level, std::make_tuple(200, 200, 1., 1.));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Dirichlet", 0.), Boundary("right", 1, "Dirichlet", 0.),
                     Boundary("upper", 2, "Dirichlet", 0.), Boundary("left", 3, "Dirichlet", 0.)};

  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  // kappa + alpha u
  const auto& alpha(0.);
  const auto& kappa(-1.);
  auto params = Parameters(Parameter("kappa", kappa), Parameter("alpha", alpha));
  // ####################
  //     variables     //
  // ####################

  auto user_func =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& v, double time) {
        const auto func = 0.;
        return func;
      });
  auto initial_condition = AnalyticalFunctions<DIM>(user_func);

  auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "c", 2, initial_condition));

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

  auto user_func_source_term =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& v, double time) {
        const auto func = 1.;
        return func;
      });

  auto src_term = AnalyticalFunctions<DIM>(user_func_source_term);
  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, params, src_term);
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(main_folder_path, "Problem1", &spatial, frequency, level_of_detail);
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
  const auto& dt = 0.1;
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
