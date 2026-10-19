/**
 * @file main.cpp
 * @brief 2D advection-diffusion in a closed domain with a non-solenoidal velocity:
 *        conservation of the integral of the advected field
 * @version 0.1
 * @date 2026-10-17
 *
 * @copyright Copyright (c) 2026
 *
 * dT/dt + div(v T) - div(k grad T) = 0   on [0,1]^2,  k = 1e-2
 *   v = (sin(pi x), sin(pi y))   (v.n = 0 on all the boundaries, div v != 0)
 *   homogeneous Neumann on all the boundaries, no Inflow coefficient
 * The total flux through the boundary is zero: int(T) must stay constant (machine precision
 * and solver tolerances), although the compressive velocity accumulates T in the upper-right
 * corner.
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "./Coefficient.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Verbose);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  using namespace Sloth2D;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  const std::string& elem_type = "InlineSquareWithQuadrangles";
  const int order_fe = 1;
  const int refinement_level = 0;
  const int NN = 32;
  const std::tuple<int, int, double, double>& tuple_of_dimensions =
      std::make_tuple(NN, NN, 1., 1.);
  SPA spatial(elem_type, order_fe, refinement_level, tuple_of_dimensions);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                     Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  Coefficient density(Glossary::Concentration, 1.);
  Coefficient heat_capacity(Glossary::Cp, 1.);
  Coefficient conductivity(Glossary::Conductivity, 1.e-2);
  Coefficient velocity_x(Glossary::Velocity_X, Scheme::Explicit, VelocityX());
  Coefficient velocity_y(Glossary::Velocity_Y, Scheme::Explicit, VelocityY());

  // ####################
  //     variables     //
  // ####################
  // Initial condition: Gaussian blob centered in the domain
  auto blob = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& x, [[maybe_unused]] double time) {
        const double r2 = (x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5);
        return std::exp(-r2 / 0.02);
      });
  auto vars =
      VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, AnalyticalFunctions<DIM>(blob)));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
      {"T", {-1.e30, 1.e30}}};
  auto p_pst = Parameters(Parameter("main_folder_path", std::string("Saves")),
                          Parameter("calculation_path", std::string("Problem1")), Parameter("frequency", 10),
                          Parameter("level_of_detail", 1),
                          Parameter("integral_to_compute", map_threshold_integral),
                          Parameter("enable_compute_energies", false));
  auto pst = PST(&spatial, p_pst);

  // ####################
  //     Problems      //
  // ####################
  Coefficients coef_pb(density, heat_capacity, conductivity, velocity_x, velocity_y);
  std::vector<SPA*> spatials{&spatial};
  TransientOPE oper(spatials, {"Fourier", "Advection"}, TimeScheme::EulerImplicit,
                    "HeatTimeDerivative");
  oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-14)));
  TransientPB problem("Advection", oper, vars, {coef_pb}, pst);

  auto cc = Coupling("Advection in a closed domain", problem);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  auto time_params = Parameters(Parameter("initial_time", 0.), Parameter("final_time", 0.5),
                                Parameter("time_step", 0.01));
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
