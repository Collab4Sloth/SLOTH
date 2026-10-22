/**
 * @file test_SlothBaseCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for SlothBaseCoefficient and Coefficient
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
#include <memory>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/CommonCoefficients.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double tol = 1.e-12;
}

TEST_CASE("Coefficient: constant value", "[Coefficients]") {
  Coefficient c(Glossary::Temperature, 300.);

  CHECK(c.is_scalar());  // constant coefficients have the Scheme::Constant
  CHECK_FALSE(c.is_implicit());
  CHECK(c.get_type() == GlossaryType::Temperature);
  CHECK(c.get_quantity().unit == GlossaryUnit::Kelvin);
  CHECK(c.get_id() == 0);

  CHECK_THAT(c.compute(), WithinRel(300., tol));
  CHECK_THAT(c.compute(std::vector{0.5}, 1), WithinRel(300., tol));
  CHECK_THAT(c.compute_gradient(0, std::vector{0.5}, 1), WithinAbs(0., tol));
}

TEST_CASE("Coefficient: time scheme", "[Coefficients]") {
  // One coefficient per scheme: exactly one predicate is true
  SECTION("Implicit") {
    const Coefficient c(Glossary::Mobility, Scheme::Implicit, W());
    CHECK(c.is_implicit());
    CHECK_FALSE(c.is_explicit());
    CHECK_FALSE(c.is_semi_implicit());
    CHECK_FALSE(c.is_scalar());
    CHECK_FALSE(c.is_file());
  }
  SECTION("Explicit") {
    const Coefficient c(Glossary::Mobility, Scheme::Explicit, W());
    CHECK(c.is_explicit());
    CHECK_FALSE(c.is_implicit());
  }
  SECTION("SemiImplicit") {
    const Coefficient c(Glossary::Mobility, Scheme::SemiImplicit, W());
    CHECK(c.is_semi_implicit());
    CHECK_FALSE(c.is_implicit());
  }
}

TEST_CASE("Coefficient: construction from a FunctionCoefficient", "[Coefficients]") {
  const std::vector<double> x{0.3};
  const double expected = W(2.).eval_f(x, 1);

  SECTION("from a temporary (moved)") {
    Coefficient c(Glossary::Mobility, Scheme::Implicit, W(2.));
    CHECK_THAT(c.compute(x, 1), WithinRel(expected, tol));
  }
  SECTION("from an lvalue (copied)") {
    const W w(2.);
    Coefficient c(Glossary::Mobility, Scheme::Implicit, w);
    CHECK_THAT(c.compute(x, 1), WithinRel(expected, tol));
  }
  SECTION("from a shared_ptr (shared)") {
    auto w = std::make_shared<W>(2.);
    Coefficient c(Glossary::Mobility, Scheme::Implicit, w);
    CHECK_THAT(c.compute(x, 1), WithinRel(expected, tol));
    // The coefficient is shared, not copied: set_time reaches the original
    c.set_time(4.);
    CHECK(w->get_time() == 4.);
  }
  SECTION("a null pointer is an error") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        Coefficient(Glossary::Mobility, Scheme::Implicit, std::shared_ptr<FunctionCoefficient>{}),
        "null coefficient pointer");
  }
}

TEST_CASE("SlothBaseCoefficient: compute overloads, with and without dimension",
          "[Coefficients]") {
  // W ignores the auxiliary and explicit values: all overloads give the same
  // result. Each call goes through a different branch of SlothBaseCoefficient.
  W w(2.);
  Coefficient c(Glossary::Mobility, Scheme::Implicit, W(2.));
  const std::vector<double> x{0.3};
  const std::vector<double> exp_x{0.9};
  const std::vector<double> aux{1000.};

  const double f = w.eval_f(x, 1);
  CHECK_THAT(c.compute(x), WithinRel(f, tol));
  CHECK_THAT(c.compute(x, 1), WithinRel(f, tol));
  CHECK_THAT(c.compute(x, aux), WithinRel(f, tol));
  CHECK_THAT(c.compute(x, aux, 1), WithinRel(f, tol));
  CHECK_THAT(c.compute(x, exp_x, aux), WithinRel(f, tol));
  CHECK_THAT(c.compute(x, exp_x, aux, 1), WithinRel(f, tol));

  const double g = w.eval_gradient(0, x, 1);
  CHECK_THAT(c.compute_gradient(0, x), WithinRel(g, tol));
  CHECK_THAT(c.compute_gradient(0, x, 1), WithinRel(g, tol));
  CHECK_THAT(c.compute_gradient(0, x, aux), WithinRel(g, tol));
  CHECK_THAT(c.compute_gradient(0, x, aux, 1), WithinRel(g, tol));
  CHECK_THAT(c.compute_gradient(0, x, exp_x, aux), WithinRel(g, tol));
  CHECK_THAT(c.compute_gradient(0, x, exp_x, aux, 1), WithinRel(g, tol));

  const double h = w.eval_hessian(0, 0, x, 1);
  CHECK_THAT(c.compute_hessian(0, 0, x), WithinRel(h, tol));
  CHECK_THAT(c.compute_hessian(0, 0, x, 1), WithinRel(h, tol));
  CHECK_THAT(c.compute_hessian(0, 0, x, aux), WithinRel(h, tol));
  CHECK_THAT(c.compute_hessian(0, 0, x, aux, 1), WithinRel(h, tol));
  CHECK_THAT(c.compute_hessian(0, 0, x, exp_x, aux), WithinRel(h, tol));
  CHECK_THAT(c.compute_hessian(0, 0, x, exp_x, aux, 1), WithinRel(h, tol));
}

TEST_CASE("SlothBaseCoefficient: name, boundaries and coordinates", "[Coefficients]") {
  Coefficient c(Glossary::Mobility, Scheme::Implicit, W());

  CHECK(c.get_name().empty());
  c.set_name("M_phi");
  CHECK(c.get_name() == "M_phi");

  CHECK(c.get_bdr_index_coef().empty());
  c.set_bdr_index_coef({1, 3});
  CHECK(c.get_bdr_index_coef() == std::vector<int>{1, 3});

  CHECK_FALSE(c.uses_coordinates());
  CHECK_NOTHROW(c.set_coordinates(std::vector{0.1, 0.2}));
}
