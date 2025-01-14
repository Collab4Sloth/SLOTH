/**
 * @file main.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Comparaison analytical and numerical solution
 * @version 0.1
 * @date 2024-11-28
 *
 * Copyright CEA (c) 2024
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

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
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  const auto DIM = 1;
  using NLFI =
      // ThermoDiffusionNLFormIntegrator<CoefficientDiscretization::Implicit, Diffusion::Constant>;
      ThermoDiffusionNLFormIntegrator<CoefficientDiscretization::Explicit, Diffusion::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = DiffusionOperator<FECollection, DIM, NLFI, Density::Constant>;
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
  double L = 1e-3;
  int NN = 50;
  SpatialDiscretization<FECollection, DIM> spatial("InlineLineWithSegments", 1, refinement_level,
                                                   std::make_tuple(NN, L));
  // ##############################
  //     Boundary conditions     //
  // // ##############################
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
  const auto& diffusionCoeff(1e-8);
  // const auto& diffusionCoeff(0.);
  //  ####################
  //      variables     //
  //  ####################

  auto user_func =
      std::function<double(const mfem::Vector&, double)>([L](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        const auto epsilon = 1e-4;
        auto func = (0.5 + 0.3 * std::tanh((xx - L / 2) / epsilon));

        const auto noiseLevel = 0.;
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::normal_distribution<> d(0, noiseLevel);

        return func + d(gen);
      });

  auto user_func_analytical = std::function<double(const mfem::Vector&, double)>(
      [L, diffusionCoeff](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        const auto epsilon = 1e-3;
        // auto func = 0.5 * (1 + std::tanh((xx - L / 2) / epsilon));
        auto func = 0.5 * (1 + std::erf((xx - L / 2) / std::sqrt(4 * diffusionCoeff * time)));

        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);
  auto analytical_solution = AnalyticalFunctions<DIM>(user_func_analytical);

  auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "c", 2, initial_condition));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves_Nx_" + std::to_string(NN);
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, TimeScheme::EulerImplicit);
  oper.overload_diffusion(Parameters(Parameter("D", diffusionCoeff)));

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst);

  PB problem1("Problem 1", oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.2;
  const auto& dt = 1e-4;
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
