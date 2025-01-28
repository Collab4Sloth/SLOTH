/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief 1D AllenCahn problem along a radius
 * @version 0.1
 * @date 2025-01-28
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  constexpr int DIM = Test<1>::dim;
  using FECollection = Test<1>::FECollection;
  using VARS = Test<1>::VARS;
  using VAR = Test<1>::VAR;
  using PSTCollection = Test<1>::PSTCollection;
  using PST = Test<1>::PST;
  using SPA = Test<1>::SPA;
  using BCS = Test<1>::BCS;
  /////////////////////////
  using NLFI = AllenCahnConstantMeltingNLFormIntegrator<
      VARS, ThermodynamicsPotentialDiscretization::Implicit, ThermodynamicsPotentials::W,
      Mobility::Constant, ThermodynamicsPotentials::H>;

  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;

  using PB = Problem<OPE, VARS, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  mfem::real_t L = 1e-3;
  size_t NN = 75;
  SPA spatial("InlineLineWithSegments", 1, refinement_level, std::make_tuple(NN, L));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

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
  //  Interface thickness
  mfem::real_t epsilon = 4 * L / NN;
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                           Parameter("lambda", lambda), Parameter("omega", omega),
                           Parameter("melting_factor", alpha));
  // ####################
  //     variables     //
  // ####################
  const auto& center_x = 0.;
  const auto& a_x = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = L / 4;

  auto user_func = std::function<double(const mfem::Vector&, double)>(
      [center_x, a_x, radius, epsilon](const mfem::Vector& x, double time) {
        const auto xx = a_x * (x[0] - center_x);
        const auto r = xx;
        const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / epsilon);
        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);
  auto vars = VARS(VAR(&spatial, bcs, "phi1", 2, initial_condition));
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  mfem::real_t iso_val = 0.5;
  // ####################
  //     operators     //
  // ####################
  std::string calculation_path = "Problem1";
  mfem::real_t iso = 0.5;
  auto p_pst1 = Parameters(
      Parameter("main_folder_path", main_folder_path),
      Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
      Parameter("level_of_detail", level_of_detail), Parameter("iso_val_to_compute", iso));

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst1);
  PB problem1(oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 50.0;
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
