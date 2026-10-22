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
