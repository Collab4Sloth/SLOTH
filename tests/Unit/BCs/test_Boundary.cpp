/**
 * @file test_Boundary.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for Boundary
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

#include "BCs/Boundary.hpp"
#include "sloth_unit.hpp"

TEST_CASE("Boundary: type without value", "[BCs]") {
  SECTION("Dirichlet") {
    const Boundary b("left", 1, "Dirichlet");
    CHECK(b.get_boundary_index() == 1);
    CHECK(b.get_boundary_type() == BoundaryConditionType::Dirichlet);
    CHECK_FALSE(b.is_periodic_boundary());
    CHECK(b.get_boundary_value() == 0.);  // default value
  }
  SECTION("Neumann and Robin are not periodic") {
    CHECK_FALSE(Boundary("top", 2, "Neumann").is_periodic_boundary());
    CHECK_FALSE(Boundary("top", 2, "Robin").is_periodic_boundary());
  }
  SECTION("Periodic (keyword case-insensitive)") {
    const Boundary b("right", 3, "periodic");
    CHECK(b.get_boundary_type() == BoundaryConditionType::Periodic);
    CHECK(b.is_periodic_boundary());
  }
  SECTION("unknown type") {
    SLOTH_CHECK_MFEM_ERROR_WITH(Boundary("left", 1, "Cauchy"), "not a valid value");
  }
}

TEST_CASE("Boundary: Dirichlet with a constant value", "[BCs]") {
  const Boundary b("left", 1, "Dirichlet", 0.25);
  CHECK(b.get_boundary_value() == 0.25);
  CHECK_FALSE(b.is_periodic_boundary());

  SECTION("only Dirichlet boundaries take a value") {
    SLOTH_CHECK_MFEM_ERROR_WITH(Boundary("left", 1, "Neumann", 1.), "only Dirichlet");
    SLOTH_CHECK_MFEM_ERROR_WITH(Boundary("left", 1, "Periodic", 1.), "only Dirichlet");
  }
}
