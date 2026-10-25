/**
 * @file test_BoundaryConditions.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the boundary conditions of a problem: Dirichlet values given by
 * coefficients, Neumann and Robin coefficients
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
#include <functional>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

namespace {
using namespace Sloth2D;

/// @brief Coefficient of the coordinates only: 1 + x + y
class OnePlusXPlusY final : public FunctionCoefficient {
 public:
  bool uses_coordinates() const override { return true; }

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [this](const std::span<const double>&, const std::span<const double>&,
                  const std::span<const double>&, const unsigned int) {
      return 1. + this->coordinates_[0] + this->coordinates_[1];
    };
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

/// @brief Coefficient attached to some boundaries (indices = attributes - 1)
Coefficient on_boundaries(Coefficient c, const std::vector<int>& bdr) {
  c.set_bdr_index_coef(bdr);
  return c;
}

/**
 * @brief Steady heat conduction - div(grad T) = 0 on the unit square
 * (4 x 4 quadrangles, Q1), solved with the boundary conditions and the
 * coefficients given (conductivity 1 included in them). Every expected solution below is affine,
 * hence exact in Q1. Boundaries: 0 lower (y = 0), 1 right (x = 1), 2 upper (y = 1), 3 left (x = 0).
 */
struct HeatSquare {
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(4, 4, 1., 1.)};
  std::vector<SPA*> spatials{&spatial};

  /// @brief L2 error of the solution against exact (all ranks)
  double solve(BCS& bcs, const Coefficients& coefs, double (*exact)(const mfem::Vector&)) {
    auto vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 1.5));
    SteadyOPE oper(spatials, {"Fourier"});
    oper.overload_nl_solver(NLSolverType::NEWTON,
                            Parameters(Parameter("description", std::string("Newton")),
                                       Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12)));
    SteadyPB problem("Heat", oper, vars, {coefs});
    auto coupling = Coupling("Heat", problem);
    auto time =
        TimeDiscretization(Parameters(Parameter("initial_time", 0.), Parameter("final_time", 1.),
                                      Parameter("time_step", 1.)),
                           coupling);
    time.solve();
    mfem::FunctionCoefficient ex(exact);
    return vars[0].get_ref_gf().ComputeL2Error(ex);
  }
};

const Coefficient conductivity(Glossary::Conductivity, 1.);

double one_plus_x(const mfem::Vector& p) { return 1. + p[0]; }
double one_plus_2x(const mfem::Vector& p) { return 1. + 2. * p[0]; }
double one_plus_x_plus_y(const mfem::Vector& p) { return 1. + p[0] + p[1]; }
}  // namespace

TEST_CASE_METHOD(HeatSquare, "BoundaryConditions: Dirichlet values given by coefficients",
                 "[BCs][Problems]") {
  SECTION("coefficient of the coordinates on the whole boundary: T = 1 + x + y") {
    BCS bcs(&spatial, Boundary("lower", 0, "Dirichlet", 0.), Boundary("right", 1, "Dirichlet", 0.),
            Boundary("upper", 2, "Dirichlet", 0.), Boundary("left", 3, "Dirichlet", 0.));
    Coefficients coefs(conductivity, on_boundaries(Coefficient(Glossary::Dirichlet,
                                                               Scheme::Implicit, OnePlusXPlusY()),
                                                   {0, 1, 2, 3}));
    CHECK(solve(bcs, coefs, one_plus_x_plus_y) < 1.e-10);
  }
  SECTION("constant coefficient: replaces the value given with the boundary") {
    BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 2.),
            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 5.));
    Coefficients coefs(conductivity, on_boundaries(Coefficient(Glossary::Dirichlet, 1.), {3}));
    CHECK(solve(bcs, coefs, one_plus_x) < 1.e-10);
  }
}

// Weak form: (grad T, grad psi) - (g, psi)_Neumann + (a T - b, psi)_Robin = 0,
// i.e. dT/dn = g on a Neumann boundary, dT/dn = b - a T on a Robin boundary.
TEST_CASE_METHOD(HeatSquare, "BoundaryConditions: Neumann and Robin coefficients",
                 "[BCs][Problems][Operators]") {
  SECTION("Neumann g = 2 on x = 1, T = 1 on x = 0: T = 1 + 2x") {
    BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.));
    Coefficients coefs(conductivity, on_boundaries(Coefficient(Glossary::Neumann, 2.), {1}));
    CHECK(solve(bcs, coefs, one_plus_2x) < 1.e-10);
  }
  SECTION("Robin a = 1, b = 3 on x = 1, T = 1 on x = 0: T = 1 + x") {
    BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Robin"),
            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.));
    Coefficients coefs(conductivity, on_boundaries(Coefficient(Glossary::Robin_a, 1.), {1}),
                       on_boundaries(Coefficient(Glossary::Robin_b, 3.), {1}));
    CHECK(solve(bcs, coefs, one_plus_x) < 1.e-10);
  }
}
