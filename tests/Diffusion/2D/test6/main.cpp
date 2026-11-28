/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Time convergence test (MMS) of Fick: implicit Euler vs BDF2 vs Crank-Nicolson
 *
 * Time convergence test (MMS) of the diffusion integrator Fick:
 *   dc/dt = div(D grad c) + S,  D = 1,   c(x, y, t) = exp(t) sin(pi x)
 * on [0, 1]^2 (Dirichlet c = 0 on x = 0 and x = 1, Neumann on y = 0 and y = 1),
 * S = (1 + pi^2) exp(t) sin(pi x) (see coefficient.json). The non-polynomial time dependence
 * gives a time discretization error which does not vanish.
 * Same problem as HeatTransfer/2D/test10 (Fourier).
 *
 * Fixed mesh: Q2, 160 x 2 elements (solution independent of y), spatial error below the
 * smallest time error. Time step divided by 2 from 0.1 to 0.0125.
 * Schemes: implicit Euler, BDF2, Crank-Nicolson (theta-scheme with theta = 1/2: parameters
 * "Fick.theta" and "SourceTerm.theta", the source term being evaluated at t^{n+1/2}).
 * Crank-Nicolson from the first step (CrankNicolson) or with a first step in implicit Euler
 * (CrankNicolsonStart, default of the theta-scheme: parameter "theta_implicit_first_step").
 * Expected slopes of the L2 error on c at t = 1: EulerImplicit: 1, BDF2: 2, CrankNicolson: 2,
 * CrankNicolsonStart: 2
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
#include <iomanip>
#include <map>
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
  const int order_fe = 2;
  const int NN = 160;  // elements in x
  const int NY = 2;    // elements in y (solution independent of y)

  // ##############################
  //  Time schemes and time steps
  // ##############################
  // Name (folder Saves_<name>_dt_<dt>), time scheme, theta of the theta-scheme, first time step
  // with implicit Euler (default of the theta-scheme)
  struct TimeSchemeCase {
    std::string name;
    TimeScheme::value scheme;
    double theta;
    bool implicit_first_step;
  };
  const std::vector<TimeSchemeCase> vect_scheme{
      {"EulerImplicit", TimeScheme::EulerImplicit, 1., false},
      {"BDF2", TimeScheme::BDF2, 1., false},
      {"CrankNicolson", TimeScheme::EulerImplicit, 0.5, false},
      {"CrankNicolsonStart", TimeScheme::EulerImplicit, 0.5, true}};
  const std::vector<double> vect_dt{0.1, 0.05, 0.025, 0.0125};

  for (const auto& [scheme_name, scheme, theta, implicit_first_step] : vect_scheme) {
    for (const auto dt : vect_dt) {
      SPA spatial(elem_type, order_fe, refinement_level, std::make_tuple(NN, NY, 1., 1.));

      // ##############################
      //     Boundary conditions     //
      // ##############################
      // u = 0 on x = 0 and x = 1 for all t (sin(pi x) = 0): time-independent conditions
      auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 0.),
                         Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
      auto bcs = BCS(&spatial, boundaries);

      // ####################
      //     variables     //
      // ####################
      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, double time) {
            return std::exp(time) * std::sin(M_PI * x[0]);
          });
      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);
      auto vars = VARS(VAR(&spatial, bcs, "c", Glossary::MoleFraction, 2, initial_condition,
                           analytical_solution));

      // ####################
      //     coefficients  //
      // ####################
      Coefficient diffusivity(Glossary::Diffusivity, 1.);
      // Source term at the time of the step (t^{n+1}), at t^{n+1/2} for Crank-Nicolson
      Coefficient source(Glossary::SourceTerm, Scheme::Explicit, SourceCoefficient());
      Coefficients coefs(diffusivity, source);

      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // Folder parsed by time_convergence.py: Saves_<scheme>_dt_<dt>
      std::ostringstream dt_str;
      dt_str << std::setprecision(10) << dt;
      const std::string& main_folder_path = "Saves_" + scheme_name + "_dt_" + dt_str.str();

      const int level_of_detail = 1;
      const int frequency = 1000000;  // no field output, only the specialized values
      const std::string calculation_path = "Problem1";
      std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
          {"c", {-10.1, 10.1}}};
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
      // CrankNicolson: Crank-Nicolson from the first step; CrankNicolsonStart: first step with
      // implicit Euler (default)
      const auto theta_params =
          Parameters(Parameter("Fick.theta", theta), Parameter("SourceTerm.theta", theta),
                     Parameter("theta_implicit_first_step", implicit_first_step));
      TransientOPE oper(spatials, {"Fick", "SourceTerm"}, theta_params, scheme, "TimeDerivative");
      oper.overload_nl_solver(
          NLSolverType::NEWTON,
          Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                     Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
      oper.overload_solver(HypreSolverType::HYPRE_GMRES,
                           Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
      oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU,
                                   Parameters(Parameter("type", 1)));
      TransientPB problem1(oper, vars, {coefs}, pst);
      auto cc = Coupling("Fick coupling", problem1);

      // ###########################################
      //            Time-integration              //
      // ###########################################
      auto time_params = Parameters(Parameter("initial_time", 0.), Parameter("final_time", 1.),
                                    Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);
      time.solve();
      Profiling::getInstance().print();
    }
  }
  MPI_Finalize();
  return 0;
}
