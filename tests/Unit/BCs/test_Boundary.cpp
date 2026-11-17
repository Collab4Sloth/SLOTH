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
#include <vector>

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

TEST_CASE("Boundary: identified by a boundary attribute set", "[BCs]") {
  const mfem::Array<int> attributes({1, 2, 3, 4});
  mfem::AttributeSets sets(attributes);
  sets.SetAttributeSet("left", mfem::Array<int>({4}));
  sets.SetAttributeSet("horizontal", mfem::Array<int>({3, 1}));

  SECTION("name only: one attribute") {
    Boundary b("left", "Dirichlet", 1.);
    CHECK_FALSE(b.has_user_index());
    CHECK(b.get_boundary_indices().empty());
    b.resolve_indices(sets);
    CHECK(b.get_boundary_indices() == std::vector<int>{3});
    CHECK(b.get_boundary_value() == 1.);
  }
  SECTION("name only: several attributes, sorted indices") {
    Boundary b("horizontal", "Neumann");
    b.resolve_indices(sets);
    CHECK(b.get_boundary_indices() == std::vector<int>{0, 2});
  }
  SECTION("name and index consistent with the set") {
    Boundary b("left", 3, "Neumann");
    b.resolve_indices(sets);
    CHECK(b.get_boundary_indices() == std::vector<int>{3});
  }
  SECTION("name and index inconsistent with the set") {
    Boundary b("left", 0, "Neumann");
    SLOTH_CHECK_MFEM_ERROR_WITH(b.resolve_indices(sets), "make it consistent");
    Boundary h("horizontal", 0, "Neumann");
    SLOTH_CHECK_MFEM_ERROR_WITH(h.resolve_indices(sets), "make it consistent");
  }
  SECTION("name without set: the index given is used") {
    Boundary b("upper", 2, "Neumann");
    b.resolve_indices(sets);
    CHECK(b.get_boundary_indices() == std::vector<int>{2});
  }
  SECTION("name without set and without index") {
    Boundary b("upper", "Neumann");
    SLOTH_CHECK_MFEM_ERROR_WITH(b.resolve_indices(sets), "'horizontal' 'left'");
  }
  SECTION("only Dirichlet boundaries take a value") {
    SLOTH_CHECK_MFEM_ERROR_WITH(Boundary("left", "Robin", 1.), "only Dirichlet");
  }
}
