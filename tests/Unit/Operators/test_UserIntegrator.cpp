/**
 * @file test_UserIntegrator.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the integrators defined outside the kernel (Extensions/UserIntegrator)
 * @version 0.1
 * @date 2026-11-09
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
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

namespace {
using namespace Sloth2D;

/**
 * @brief Heat conduction on the unit square (4 x 4 quadrangles, Q1), conductivity 1.
 * Boundaries: 0 lower (y = 0), 1 right (x = 1), 2 upper (y = 1), 3 left (x = 0).
 */
struct HeatSquare {
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(4, 4, 1., 1.)};
  std::vector<SPA*> spatials{&spatial};
  Coefficients coefs{Coefficient(Glossary::Conductivity, 1.)};
  Parameters newton{Parameter("description", std::string("Newton")), Parameter("rel_tol", 1.e-12),
                    Parameter("abs_tol", 1.e-12)};

  /// @brief Largest absolute difference between two vectors of true dofs (all ranks)
  static double max_difference(const mfem::Vector& a, const mfem::Vector& b) {
    REQUIRE(a.Size() == b.Size());
    double local = 0., global = 0.;
    for (int i = 0; i < a.Size(); ++i) local = std::max(local, std::abs(a(i) - b(i)));
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
  }
};
}  // namespace

TEST_CASE("UserIntegrator: name of an integrator entry", "[Operators][Extensions]") {
  using namespace Sloth2D;
  const IntegratorEntry<VARS> built_in = "Fourier";
  const IntegratorEntry<VARS> user = user_integrator<FourierNLFormIntegrator, VARS>("UserFourier");
  CHECK(integrator_name<VARS>(built_in) == "Fourier");
  CHECK(integrator_name<VARS>(user) == "UserFourier");

  // Aliases of the Sloth2D namespace
  const IntegratorList list{"AllenCahn", CustomIntegrator<FourierNLFormIntegrator>("Custom")};
  REQUIRE(list.size() == 2);
  CHECK(integrator_name<VARS>(list[0]) == "AllenCahn");
  CHECK(integrator_name<VARS>(list[1]) == "Custom");
}

// Steady heat equation, T = 1 on x = 0 and T = 2 on x = 1: exact solution T = 1 + x (exact in Q1)
// with the built-in Fourier integrator class declared as a custom integrator (Sloth2D alias)
TEST_CASE_METHOD(HeatSquare, "UserIntegrator: right-hand side integrator of a steady operator",
                 "[Operators][Extensions]") {
  BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 2.),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.));
  auto vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 1.5));

  SteadyOPE oper(spatials, {CustomIntegrator<FourierNLFormIntegrator>("UserFourier")});
  oper.overload_nl_solver(NLSolverType::NEWTON, newton);
  SteadyPB problem("Heat", oper, vars, {coefs});
  auto coupling = Coupling("Heat", problem);
  auto time =
      TimeDiscretization(Parameters(Parameter("initial_time", 0.), Parameter("final_time", 1.),
                                    Parameter("time_step", 1.)),
                         coupling);
  time.solve();

  mfem::FunctionCoefficient exact([](const mfem::Vector& p) { return 1. + p[0]; });
  CHECK(vars[0].get_ref_gf().ComputeL2Error(exact) < 1.e-10);
}

// Transient heat equation, T = sin(pi x) at t = 0, T = 0 on x = 0 and x = 1: two implicit Euler
// steps with the built-in time derivative and with the same class used as a user integrator
TEST_CASE_METHOD(HeatSquare, "UserIntegrator: time derivative integrator of a transient operator",
                 "[Operators][Extensions]") {
  BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 0.),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.));
  const auto initial_condition =
      AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& p, double) { return std::sin(M_PI * p[0]); }));

  auto run = [&](const IntegratorEntry<VARS>& time_derivative) {
    auto vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, initial_condition));
    TransientOPE oper(spatials, {"Fourier"}, TimeScheme::EulerImplicit, time_derivative);
    oper.overload_nl_solver(NLSolverType::NEWTON, newton);
    TransientPB problem("Heat", oper, vars, {coefs});
    auto coupling = Coupling("Heat", problem);
    auto time =
        TimeDiscretization(Parameters(Parameter("initial_time", 0.), Parameter("final_time", 0.02),
                                      Parameter("time_step", 0.01)),
                           coupling);
    time.solve();
    return vars[0].get_unknown();
  };

  const mfem::Vector built_in = run("TimeDerivative");
  const mfem::Vector user = run(user_integrator<TimeNLFormIntegrator, VARS>("UserTimeDerivative"));
  CHECK(max_difference(built_in, user) < 1.e-12);
}

// A builder returning no integrator is reported with the name of the user integrator
TEST_CASE_METHOD(HeatSquare, "UserIntegrator: builder returning no integrator",
                 "[Operators][Extensions]") {
  BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann"));
  auto vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 1.));

  const IntegratorEntry<VARS> broken = UserIntegrator<VARS>(
      "Broken",
      [](const UserIntegratorArgs<VARS>&) -> SlothNLFormIntegrator<VARS>* { return nullptr; });
  SteadyOPE oper(spatials, {broken});
  const std::vector<mfem::Vector> u{vars[0].get_unknown()};
  auto build = [&]() { return oper.set_nlfi_ptr(1., broken, u); };
  SLOTH_CHECK_MFEM_ERROR_WITH(build(), "Broken");
}
