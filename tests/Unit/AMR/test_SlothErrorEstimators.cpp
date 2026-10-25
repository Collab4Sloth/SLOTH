/**
 * @file test_SlothErrorEstimators.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for SlothErrorEstimators (AMR)
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
#include <memory>

#include "AMR/SlothErrorEstimators.hpp"
#include "mfem.hpp"
#include "sloth_unit.hpp"

namespace {
double affine(const mfem::Vector& p) { return 1. + 2. * p[0] - 3. * p[1]; }
double quadratic(const mfem::Vector& p) { return p[0] * p[0] + p[1] * p[1]; }

/// @brief Largest local error over all the ranks (collective)
double max_error(mfem::ErrorEstimator& estimator) {
  const mfem::Vector& errors = estimator.GetLocalErrors();
  double local = errors.Size() > 0 ? errors.Max() : 0.;
  double global = 0.;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}
}  // namespace

// Both estimators measure the discontinuity of the flux grad(u) between
// elements: zero for an affine field (constant gradient), not for x^2 + y^2.
TEST_CASE_METHOD(sloth_unit::UnitSquare, "SlothErrorEstimators: Kelly and Zienkiewicz-Zhu",
                 "[AMR]") {
  const auto type = GENERATE(ErrorEstimatorType::KELLY, ErrorEstimatorType::L2_ZIENKIEWICZ_ZHU);
  CAPTURE(static_cast<int>(type));
  mfem::DiffusionIntegrator integ;  // flux: grad(u)
  SlothErrorEstimators estimators(type, &integ);

  SECTION("affine field: zero error") {
    auto u = project(affine);
    CHECK(max_error(*estimators.get_value(u, *fes, *pmesh)) < 1.e-10);
  }
  SECTION("quadratic field: non-zero error") {
    auto u = project(quadratic);
    CHECK(max_error(*estimators.get_value(u, *fes, *pmesh)) > 1.e-4);
  }
  SECTION("flux spaces shared by the fields of the same mesh and order") {
    // Kelly: L2 flux; Zienkiewicz-Zhu: L2 flux and RT smoothed flux
    const std::size_t spaces_per_order = type == ErrorEstimatorType::KELLY ? 1 : 2;
    // Fields created and destroyed one after the other: nothing is kept per field
    for (int k = 0; k < 10; ++k) {
      auto u = project(k % 2 == 0 ? affine : quadratic);
      const double error = max_error(*estimators.get_value(u, *fes, *pmesh));
      CHECK((k % 2 == 0 ? error < 1.e-10 : error > 1.e-4));
    }
    CHECK(estimators.number_of_cached_spaces() == spaces_per_order);

    // Field of order 2 on the same mesh: other flux spaces
    mfem::H1_FECollection fec2(2, 2);
    mfem::ParFiniteElementSpace fes2(pmesh.get(), &fec2);
    mfem::ParGridFunction u2(&fes2);
    mfem::FunctionCoefficient c(affine);
    u2.ProjectCoefficient(c);
    CHECK(max_error(*estimators.get_value(u2, fes2, *pmesh)) < 1.e-10);
    CHECK(estimators.number_of_cached_spaces() == 2 * spaces_per_order);
  }
}

TEST_CASE_METHOD(sloth_unit::UnitSquare, "SlothErrorEstimators: errors", "[AMR][errors]") {
  auto u = project(affine);
  SECTION("Kelly needs an integrator") {
    SlothErrorEstimators estimators(ErrorEstimatorType::KELLY);
    SLOTH_CHECK_MFEM_ERROR_WITH(estimators.get_value(u, *fes, *pmesh), "requires an integrator");
  }
  SECTION("Zienkiewicz-Zhu needs an integrator") {
    SlothErrorEstimators estimators(ErrorEstimatorType::L2_ZIENKIEWICZ_ZHU);
    SLOTH_CHECK_MFEM_ERROR_WITH(estimators.get_value(u, *fes, *pmesh), "requires an integrator");
  }
  SECTION("unknown estimator") {
    SlothErrorEstimators estimators(static_cast<ErrorEstimatorType>(99));
    SLOTH_CHECK_MFEM_ERROR_WITH(estimators.get_value(u, *fes, *pmesh), "unknown estimator");
  }
}
