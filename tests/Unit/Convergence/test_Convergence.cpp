/**
 * @file test_Convergence.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for PhysicalConvergence and Convergence
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
#include <tuple>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;
using sloth_unit::point;

// The criterion is local to the MPI rank: ProblemBase::check_convergence
// reduces it over the ranks (MPI_Allreduce).

TEST_CASE("PhysicalConvergence: ABSOLUTE_MAX = max |x - xn|", "[Convergence]") {
  PhysicalConvergence cvg(ConvergenceType::ABSOLUTE_MAX, 1.e-3);

  SECTION("not converged") {
    const auto [converged, criterion] =
        cvg.getPhysicalConvergence(point({1., 2., 3.}), point({1., 2.1, 3.}));
    CHECK_FALSE(converged);
    CHECK_THAT(criterion, WithinRel(0.1, 1.e-12));
  }
  SECTION("converged") {
    const auto [converged, criterion] =
        cvg.getPhysicalConvergence(point({1., 2., 3.}), point({1., 2.0005, 3.}));
    CHECK(converged);
    CHECK_THAT(criterion, WithinRel(5.e-4, 1.e-9));
  }
  SECTION("the threshold itself is not converged (strict inequality)") {
    PhysicalConvergence half(ConvergenceType::ABSOLUTE_MAX, 0.5);
    CHECK_FALSE(std::get<0>(half.getPhysicalConvergence(point({1.}), point({0.5}))));
  }
}

TEST_CASE("PhysicalConvergence: RELATIVE_MAX = max |x - xn| / max(1e-12, |xn|)",
          "[Convergence]") {
  PhysicalConvergence cvg(ConvergenceType::RELATIVE_MAX, 1.e-2);

  SECTION("relative to the previous value") {
    // |1010 - 1000| / 1000 = 1e-2 (not converged: strict), |2 - 1| / 1 = 1
    const auto [converged, criterion] =
        cvg.getPhysicalConvergence(point({1010., 2.}), point({1000., 1.}));
    CHECK_FALSE(converged);
    CHECK_THAT(criterion, WithinRel(1., 1.e-12));
  }
  SECTION("converged") {
    const auto [converged, criterion] =
        cvg.getPhysicalConvergence(point({1001., -2.002}), point({1000., -2.}));
    CHECK(converged);
    CHECK_THAT(criterion, WithinRel(1.e-3, 1.e-9));
  }
  SECTION("previous value 0: divided by 1e-12, no division by zero") {
    const auto [converged, criterion] =
        cvg.getPhysicalConvergence(point({5.e-15}), point({0.}));
    CHECK(converged);  // 5e-15 / 1e-12 = 5e-3 < 1e-2
    CHECK_THAT(criterion, WithinRel(5.e-3, 1.e-9));
  }
}

TEST_CASE("PhysicalConvergence: no local value (rank without dof)", "[Convergence]") {
  const mfem::Vector empty;
  for (const auto type : {ConvergenceType::ABSOLUTE_MAX, ConvergenceType::RELATIVE_MAX}) {
    CAPTURE(static_cast<int>(type));
    PhysicalConvergence cvg(type, 1.e-3);
    const auto [converged, criterion] = cvg.getPhysicalConvergence(empty, empty);
    CHECK(converged);
    CHECK(criterion == 0.);
  }
}

TEST_CASE("PhysicalConvergence: unknown criterion", "[Convergence][errors]") {
  SLOTH_CHECK_MFEM_ERROR_WITH(
      PhysicalConvergence(static_cast<ConvergenceType::value>(99), 1.e-3), "available");
}

TEST_CASE("Convergence: one criterion per variable, in order", "[Convergence]") {
  SECTION("default: no criterion") { CHECK(Convergence().getPhysicalConvergence().empty()); }

  SECTION("several criteria") {
    const Convergence cvg(PhysicalConvergence(ConvergenceType::ABSOLUTE_MAX, 1.),
                          PhysicalConvergence(ConvergenceType::RELATIVE_MAX, 1.e-6));
    auto criteria = cvg.getPhysicalConvergence();
    REQUIRE(criteria.size() == 2);
    // Same increment 0.5 on a value 10: absolute 0.5 < 1, relative 0.05 > 1e-6
    CHECK(std::get<0>(criteria[0].getPhysicalConvergence(point({10.5}), point({10.}))));
    CHECK_FALSE(std::get<0>(criteria[1].getPhysicalConvergence(point({10.5}), point({10.}))));
  }
}
