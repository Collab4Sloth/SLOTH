/**
 * @file test_SteadyAndExplicit.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Complete problems: steady operator (SteadyOPE / SteadyPB) and explicit time scheme
 * @version 0.1
 * @date 2026-10-22
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
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <functional>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;

namespace {
using namespace Sloth2D;

/// @brief Integral of a variable over the domain (all ranks)
double integral(VAR& v) {
  mfem::ParLinearForm b(v.get_fespace());
  mfem::ConstantCoefficient one(1.);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  b.Assemble();
  return b(v.get_ref_gf());
}

/// @brief Largest value of a variable (all ranks)
double max_value(VAR& v) {
  const mfem::Vector u = v.get_unknown();
  double local = u.Size() > 0 ? u.Max() : -1.e300, global = 0.;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

Parameters time_parameters(double t_final, double dt) {
  return Parameters(Parameter("initial_time", 0.), Parameter("final_time", t_final),
                    Parameter("time_step", dt));
}
}  // namespace

// -----------------------------------------------------------------------------
// Steady heat conduction: - div(k grad T) = 0, T = 1 on x = 0, T = 2 on x = 1,
// insulated on y = 0 and y = 1. Exact solution T = 1 + x, exact in Q1.
// -----------------------------------------------------------------------------
TEST_CASE("Problems: steady heat conduction (SteadyOPE, SteadyPB)", "[Problems][Operators]") {
  SPA spatial("InlineSquareWithQuadrangles", 1, 0, std::make_tuple(8, 8, 1., 1.));
  auto bcs = BCS(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 2.),
                 Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.));
  auto vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 1.5));  // initial guess

  Coefficients coefs(Coefficient(Glossary::Conductivity, 3.));
  std::vector<SPA*> spatials{&spatial};
  SteadyOPE oper(spatials, {"Fourier"});
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", std::string("Newton")),
                                     Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12)));
  SteadyPB problem("Heat", oper, vars, {coefs});
  auto coupling = Coupling("Heat", problem);
  auto time = TimeDiscretization(time_parameters(0.1, 0.1), coupling);
  time.solve();

  mfem::FunctionCoefficient exact([](const mfem::Vector& p) { return 1. + p[0]; });
  CHECK(vars[0].get_ref_gf().ComputeL2Error(exact) < 1.e-10);
}

// -----------------------------------------------------------------------------
// Diffusion with the explicit Euler scheme, insulated square:
//   dc/dt = div(D grad c), c0 = 0.5 + 0.2 cos(pi x)
// The integral of c is conserved exactly (lumped mass, rows of the stiffness
// matrix summing to zero) and the perturbation decays.
// -----------------------------------------------------------------------------
TEST_CASE("Problems: explicit Euler diffusion conserves the integral", "[Problems][Operators]") {
  SPA spatial("InlineSquareWithQuadrangles", 1, 0, std::make_tuple(8, 8, 1., 1.));
  auto bcs = BCS(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                 Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann"));
  auto c0 = AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& p, double) { return 0.5 + 0.2 * std::cos(M_PI * p[0]); }));
  auto vars = VARS(VAR(&spatial, bcs, "c", Glossary::Concentration, 2, c0));

  const double initial_integral = integral(vars[0]);
  const double initial_max = max_value(vars[0]);

  Coefficients coefs(Coefficient(Glossary::Diffusivity, 1.));
  std::vector<SPA*> spatials{&spatial};
  TransientOPE oper(spatials, {"Fick"}, TimeScheme::EulerExplicit, "TimeDerivative");
  TransientPB problem("Diffusion", oper, vars, {coefs});
  auto coupling = Coupling("Diffusion", problem);
  // h = 1/8, D = 1: explicit stability requires dt of order h^2 / 4 = 4e-3
  auto time = TimeDiscretization(time_parameters(5.e-3, 1.e-3), coupling);
  time.solve();

  CHECK_THAT(integral(vars[0]), WithinRel(initial_integral, 1.e-10));
  CHECK(max_value(vars[0]) < initial_max);
}

// -----------------------------------------------------------------------------
// SteadyOperator::Mult is not used: the steady solve goes through
// SteadyPhaseFieldReducedOperator (Mult is required by mfem::Operator).
// Calling it is an error.
// -----------------------------------------------------------------------------
TEST_CASE("SteadyOperator: Mult is not used", "[Problems][Operators][errors]") {
  SPA spatial("InlineSquareWithQuadrangles", 1, 0, std::make_tuple(2, 2, 1., 1.));
  std::vector<SPA*> spatials{&spatial};
  SteadyOPE oper(spatials, {"Fourier"});
  const mfem::Vector k(4);
  mfem::Vector y(4);
  SLOTH_CHECK_MFEM_ERROR_WITH(oper.Mult(k, y), "SteadyOperator::Mult is not used");
}

// -----------------------------------------------------------------------------
// Time of the coefficients at each stage of the time scheme.
// du/dt = f(t) = t, u(0) = 0, uniform field (no spatial error), t in [0, 1], dt = 0.25.
// The MFEM ODE solvers set the time of each stage (TransientOperator::SetTime):
//   implicit Euler:   u_{n+1} = u_n + dt f(t_{n+1})  ->  u(1) = dt^2 N (N + 1) / 2 = 0.625
//   explicit Euler:   u_{n+1} = u_n + dt f(t_n)      ->  u(1) = dt^2 N (N - 1) / 2 = 0.375
//   SDIRK23, SDIRK33, RK4 (order >= 2): exact for f linear in t -> u(1) = 1/2
// -----------------------------------------------------------------------------
namespace {
/// @brief Coefficient equal to the time: f(t) = t
class TimeValue final : public FunctionCoefficient {
 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [this](const std::span<const double>&, const std::span<const double>&,
                  const std::span<const double>&, const unsigned int) { return this->get_time(); };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  GradientF() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&,
              const unsigned int) { return std::vector<double>(u.size(), 0.); };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  HessianF() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&,
              const unsigned int) { return std::vector<double>(u.size() * u.size(), 0.); };
  }
};
}  // namespace

TEST_CASE("Problems: coefficients evaluated at the time of each stage", "[Problems][Operators]") {
  const auto [time_scheme, expected] =
      GENERATE(std::pair{TimeScheme::EulerImplicit, 0.625},
               std::pair{TimeScheme::EulerExplicit, 0.375}, std::pair{TimeScheme::SDIRK23, 0.5},
               std::pair{TimeScheme::SDIRK33, 0.5}, std::pair{TimeScheme::RungeKutta4, 0.5});
  CAPTURE(static_cast<int>(time_scheme), expected);

  SPA spatial("InlineSquareWithQuadrangles", 1, 0, std::make_tuple(2, 2, 1., 1.));
  BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann"));
  auto vars = VARS(VAR(&spatial, bcs, "u", Glossary::Concentration, 2, 0.));
  Coefficients coefs(Coefficient(Glossary::SourceTerm, Scheme::Explicit, TimeValue()));
  std::vector<SPA*> spatials{&spatial};
  TransientOPE oper(spatials, {"SourceTerm"}, time_scheme, "TimeDerivative");
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", std::string("Newton")),
                                     Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12)));
  TransientPB problem("u", oper, vars, {coefs});
  auto coupling = Coupling("u", problem);
  auto time =
      TimeDiscretization(Parameters(Parameter("initial_time", 0.), Parameter("final_time", 1.),
                                    Parameter("time_step", 0.25)),
                         coupling);
  time.solve();

  const mfem::Vector u = vars[0].get_unknown();
  for (int i = 0; i < u.Size(); ++i) {
    CAPTURE(i);
    CHECK_THAT(u(i), Catch::Matchers::WithinAbs(expected, 1.e-10));
  }
}
