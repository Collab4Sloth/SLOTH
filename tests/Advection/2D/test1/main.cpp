/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Transport of a tanh interface without diffusion (inflow and outflow boundaries):
 * implicit Euler vs Crank-Nicolson, comparison with the exact solution
 *
 *   dphi/dt + div(v phi) = 0,   v = (1, 0),
 *   phi(x, y, t) = (1 - tanh((x - x0 - t) / w)) / 2,   x0 = 0.25,  w = 0.1,
 * on [0, 1]^2, t in [0, 0.5]: the interface goes from x = 0.25 to x = 0.75.
 * Boundaries (Inflow coefficient attached to all the boundaries, inflow/outflow switch from the
 * sign of v.n): left (x = 0) inflow with the exact solution, right (x = 1) outflow, lower and upper
 * walls (v.n = 0).
 * Implicit Euler adds a numerical diffusion of order dt (the interface is smeared).
 * Crank-Nicolson (theta-scheme with theta = 1/2, parameter "Advection.theta", applied also to the
 * inflow/outflow terms) is of order 2 and non-dissipative.
 *
 * Fixed mesh: Q3, 128 x 2 elements (solution independent of y), spatial error below the smallest
 * time error. Time step divided by 2 from 0.01 to 0.00125.
 * Expected slopes of the L2 error on phi at t = 0.5: EulerImplicit: 1, CrankNicolson: 2
 * (see tests/tools/time_convergence.py).
 * @version 0.1
 * @date 2026-11-28
 *
 * Copyright CEA (C) 2026
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "./Coefficient.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Verbose);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  Profiling::getInstance().enable();
  using namespace Sloth2D;

  // ##############################
  //  Fixed spatial discretization
  // ##############################
  const int refinement_level = 0;
  const std::string elem_type = "InlineSquareWithQuadrangles";
  const int order_fe = 3;
  const int NN = 128;      // elements in x
  const int NY = 2;        // elements in y (solution independent of y)
  const double x0 = 0.25;  // initial position of the interface
  const double w = 0.1;    // width parameter of the tanh profile

  // ##############################
  //  Time schemes and time steps
  // ##############################
  // Name (folder Saves_<name>_dt_<dt>), theta of the theta-scheme (implicit Euler time scheme)
  const std::vector<std::pair<std::string, double>> vect_scheme{{"EulerImplicit", 1.},
                                                                {"CrankNicolson", 0.5}};
  const std::vector<double> vect_dt{0.01, 0.005, 0.0025, 0.00125};

  for (const auto& [scheme_name, theta] : vect_scheme) {
    for (const auto dt : vect_dt) {
      SPA spatial(elem_type, order_fe, refinement_level, std::make_tuple(NN, NY, 1., 1.));

      // ##############################
      //     Boundary conditions     //
      // ##############################
      // No condition for the variable: inflow/outflow through the Inflow coefficient
      auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                         Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
      auto bcs = BCS(&spatial, boundaries);

      // ####################
      //     variables     //
      // ####################
      // Exact solution: interface transported at the velocity 1 in x
      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [x0, w](const mfem::Vector& x, double time) {
            return 0.5 * (1. - std::tanh((x[0] - x0 - time) / w));
          });
      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);
      auto vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition,
                           analytical_solution));

      // ####################
      //     coefficients  //
      // ####################
      // Imposed velocity and inflow value (exact solution, evaluated at t^{n+1/2} for
      // Crank-Nicolson)
      Coefficient velocity_x(Glossary::Velocity_X, Scheme::Explicit, VelocityX());
      Coefficient velocity_y(Glossary::Velocity_Y, Scheme::Explicit, VelocityY());
      Coefficient inflow(Glossary::Inflow, Scheme::Explicit, InflowCoefficient());
      inflow.set_bdr_index_coef(std::vector<int>{0, 1, 2, 3});
      Coefficients coefs(velocity_x, velocity_y, inflow);

      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // Folder parsed by time_convergence.py: Saves_<scheme>_dt_<dt>
      std::ostringstream dt_str;
      dt_str << std::setprecision(10) << dt;
      const std::string& main_folder_path = "Saves_" + scheme_name + "_dt_" + dt_str.str();

      const int level_of_detail = 1;
      // Fields saved every 0.05 (10 snapshots of the transported interface)
      const int frequency = static_cast<int>(std::lround(0.05 / dt));
      const std::string calculation_path = "Problem1";
      std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
          {"phi", {-10.1, 10.1}}};
      auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                              Parameter("calculation_path", calculation_path),
                              Parameter("frequency", frequency),
                              Parameter("level_of_detail", level_of_detail),
                              Parameter("integral_to_compute", map_threshold_integral),
                              Parameter("enable_save_specialized_at_iter", true),
                              Parameter("enable_compute_energies", false));
      auto pst = PST(&spatial, p_pst);

      // ####################
      //     operators     //
      // ####################
      std::vector<SPA*> spatials{&spatial};
      // Crank-Nicolson from the first step (no implicit Euler first step): same results as the
      // references
      const auto theta_params = Parameters(Parameter("Advection.theta", theta),
                                           Parameter("theta_implicit_first_step", false));
      TransientOPE oper(spatials, {"Advection"}, theta_params, TimeScheme::EulerImplicit,
                        "TimeDerivative");
      oper.overload_nl_solver(
          NLSolverType::NEWTON,
          Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                     Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
      oper.overload_solver(HypreSolverType::HYPRE_GMRES,
                           Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
      oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU,
                                   Parameters(Parameter("type", 1)));
      TransientPB problem1(oper, vars, {coefs}, pst);
      auto cc = Coupling("Advection coupling", problem1);

      // ###########################################
      //            Time-integration              //
      // ###########################################
      auto time_params = Parameters(Parameter("initial_time", 0.), Parameter("final_time", 0.5),
                                    Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);
      time.solve();
      Profiling::getInstance().print();
    }
  }
  MPI_Finalize();
  return 0;
}
