/**
 * @file test_Parameters.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for Parameters
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
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Parameters/Parameters.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE("Parameters: get_param_value", "[Parameters]") {
  const Parameters params(Parameter("nb_iter", 10), Parameter("epsilon", 1.e-3));

  CHECK(params.get_size() == 2);
  CHECK(params.get_param_value<int>("nb_iter") == 10);
  CHECK_THAT(params.get_param_value<double>("epsilon"), WithinRel(1.e-3, 1.e-12));
}

TEST_CASE("Parameters: get_param_value errors", "[Parameters][errors]") {
  const Parameters params(Parameter("nb_iter", 10));

  SECTION("missing parameter") {
    SLOTH_CHECK_MFEM_ERROR_WITH(params.get_param_value<int>("missing"), "Not Found");
  }
  SECTION("wrong type: no implicit int -> double conversion") {
    SLOTH_CHECK_MFEM_ERROR_WITH(params.get_param_value<double>("nb_iter"), "Invalid conversion");
  }
}

TEST_CASE("Parameters: construction from a std::vector", "[Parameters]") {
  SECTION("non-empty") {
    const std::vector<Parameter> v{Parameter("a", 1), Parameter("b", 2.)};
    const Parameters params(v);
    CHECK(params.get_size() == 2);
    REQUIRE(params.get_vector().size() == 2);
    CHECK(params.get_vector().front().get_name() == "a");  // order is kept
  }
  SECTION("empty") {
    const Parameters params(std::vector<Parameter>{});
    CHECK(params.get_size() == 0);
    CHECK(params.get_vector().empty());
  }
}

TEST_CASE("Parameters: has_parameter and get_parameter", "[Parameters]") {
  const Parameters params(Parameter("nb_iter", 10), Parameter("epsilon", 1.e-3));

  CHECK(params.has_parameter("nb_iter"));
  CHECK_FALSE(params.has_parameter("NB_ITER"));  // parameter names are case-sensitive
  CHECK_FALSE(params.has_parameter("missing"));

  const auto found = params.get_parameter("epsilon");
  REQUIRE(found.has_value());
  CHECK(found->get_name() == "epsilon");
  CHECK_FALSE(params.get_parameter("missing").has_value());
}

TEST_CASE("Parameters: get_param_value_or_default", "[Parameters]") {
  const Parameters params(Parameter("nb_iter", 10), Parameter("verbose", true));

  CHECK(params.get_param_value_or_default<int>("nb_iter", 3) == 10);  // present
  CHECK(params.get_param_value_or_default<int>("missing", 3) == 3);   // default
  CHECK(params.get_param_value_or_default<bool>("verbose", false));
  CHECK(params.get_param_value_or_default<std::string>("missing", "none") == "none");

  SECTION("wrong type is an error, even with a default value") {
    SLOTH_CHECK_MFEM_ERROR_WITH(params.get_param_value_or_default<double>("nb_iter", 1.),
                                "Invalid conversion");
  }
}

TEST_CASE("Parameters: operator+ with a Parameter", "[Parameters]") {
  const Parameters params(Parameter("nb_iter", 10), Parameter("epsilon", 1.e-3));

  SECTION("a new name is appended") {
    const auto sum = params + Parameter("theta", 0.5);
    CHECK(sum.get_size() == 3);
    CHECK_THAT(sum.get_param_value<double>("theta"), WithinRel(0.5, 1.e-12));
  }
  SECTION("an existing name is replaced, not duplicated") {
    const auto sum = params + Parameter("nb_iter", 20);
    CHECK(sum.get_size() == 2);
    CHECK(sum.get_param_value<int>("nb_iter") == 20);
  }
  SECTION("the replacement may change the type") {
    const auto sum = params + Parameter("nb_iter", 20.);
    CHECK_THAT(sum.get_param_value<double>("nb_iter"), WithinRel(20., 1.e-12));
  }
  SECTION("the left operand is not modified") {
    [[maybe_unused]] const auto sum = params + Parameter("nb_iter", 20);
    CHECK(params.get_param_value<int>("nb_iter") == 10);
  }
}

TEST_CASE("Parameters: operator+ with Parameters (right operand wins)", "[Parameters]") {
  const Parameters left(Parameter("nb_iter", 10), Parameter("epsilon", 1.e-3));
  const Parameters right(Parameter("nb_iter", 50), Parameter("theta", 0.5));

  const auto sum = left + right;
  CHECK(sum.get_size() == 3);
  CHECK(sum.get_param_value<int>("nb_iter") == 50);  // from the right operand
  CHECK_THAT(sum.get_param_value<double>("epsilon"), WithinRel(1.e-3, 1.e-12));
  CHECK_THAT(sum.get_param_value<double>("theta"), WithinRel(0.5, 1.e-12));

  SECTION("empty operands") {
    const Parameters empty(std::vector<Parameter>{});
    CHECK((left + empty).get_size() == left.get_size());
    CHECK((empty + right).get_size() == right.get_size());
  }
}

TEST_CASE("Parameters: operator-", "[Parameters]") {
  const Parameters params(Parameter("nb_iter", 10), Parameter("epsilon", 1.e-3));

  SECTION("removes by name, whatever the value") {
    const auto diff = params - Parameter("nb_iter", 0);
    CHECK(diff.get_size() == 1);
    CHECK_FALSE(diff.has_parameter("nb_iter"));
  }
  SECTION("removing a missing name changes nothing") {
    const auto diff = params - Parameter("missing", 0);
    CHECK(diff.get_size() == 2);
  }
}
