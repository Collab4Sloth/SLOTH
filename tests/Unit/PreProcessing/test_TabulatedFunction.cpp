/**
 * @file test_TabulatedFunction.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for TabulatedFunction (structured CSV grids)
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
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "PreProcessing/TabulatedFunction.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE("TabulatedFunction: reading and interpolation", "[PreProcessing]") {
  // F = t + 10 x on {0, 1} x {0, 1, 2}; unsorted lines and blank line on purpose
  const sloth_unit::TemporaryFile csv("tab.csv",
                                      "t,x,F\n"
                                      "1,2,21\n"
                                      "0,0,0\n"
                                      "\n"
                                      "0,1,10\n"
                                      "0,2,20\n"
                                      "1,0,1\n"
                                      "1,1,11\n");
  const TabulatedFunction tf(csv.path());
  CHECK(tf.input_dim() == 1);
  CHECK(tf.has_data());
  CHECK_THAT(tf.interpolate(0.5, std::vector{1.5}), WithinRel(15.5, 1.e-12));

  SECTION("mfem::Vector overload") {
    mfem::Vector x(1);
    x[0] = 0.25;
    CHECK_THAT(tf.interpolate(1., x), WithinRel(3.5, 1.e-12));
  }
  SECTION("wrong number of coordinates") {
    SLOTH_CHECK_MFEM_ERROR(tf.interpolate(0., std::vector{0.5, 0.5}));
  }
}

TEST_CASE("TabulatedFunction: invalid files", "[PreProcessing]") {
  SECTION("missing file") {
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction("/nonexistent/sloth_unit.csv"),
                                "cannot open CSV file");
  }
  SECTION("header only") {
    const sloth_unit::TemporaryFile csv("header.csv", "t,x,F\n");
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction(csv.path()), "no data line");
  }
  SECTION("too many input columns (more than 3 spatial coordinates)") {
    const sloth_unit::TemporaryFile csv("cols.csv", "t,a,b,c,d,F\n0,0,0,0,0,1\n");
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction(csv.path()), "not well read");
  }
  SECTION("inconsistent number of columns") {
    const sloth_unit::TemporaryFile csv("ragged.csv", "t,x,F\n0,0,1\n0,1\n");
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction(csv.path()), "line 3 contains 2 columns");
  }
  SECTION("missing grid node (not a structured grid)") {
    const sloth_unit::TemporaryFile csv("holes.csv", "t,x,F\n0,0,0\n0,1,1\n1,0,1\n");
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction(csv.path()), "not a structured grid");
  }
  SECTION("duplicate grid point (and as many lines as grid nodes)") {
    const sloth_unit::TemporaryFile csv("dup.csv", "t,x,F\n0,0,0\n0,0,1\n0,1,1\n1,0,1\n");
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction(csv.path()), "duplicate grid point");
  }
  SECTION("value that is not a number") {
    const sloth_unit::TemporaryFile csv("nan.csv", "t,x,F\n0,0,abc\n");
    SLOTH_CHECK_MFEM_ERROR_WITH(TabulatedFunction(csv.path()), "'abc' is not a valid number");
  }
}
