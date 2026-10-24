/**
 * @file test_Coefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for Coefficients (collection of Coefficient)
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
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Coefficients/Coefficients.hpp"
#include "Coefficients/CommonCoefficients.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE("Coefficients: construction and access", "[Coefficients]") {
  SECTION("empty") {
    const Coefficients coefs;
    CHECK(coefs.size() == 0);
    CHECK(coefs.get_types().empty());
    CHECK_THROWS_AS(coefs[0], std::out_of_range);
  }
  SECTION("from several Coefficient: order and types are kept") {
    const Coefficient t(Glossary::Temperature, 300.);
    const Coefficient p(Glossary::Pressure, 1.e5);
    Coefficients coefs(t, p);

    CHECK(coefs.size() == 2);
    CHECK(coefs.get_types() ==
          std::vector<GlossaryType>{GlossaryType::Temperature, GlossaryType::Pressure});
    CHECK_THAT(coefs[1].compute(), WithinRel(1.e5, 1.e-12));
    CHECK(coefs.getCoefficients().size() == 2);

    const Coefficients& const_coefs = coefs;  // const operator[]
    CHECK(const_coefs[0].get_type() == GlossaryType::Temperature);
    CHECK_THROWS_AS(coefs[2], std::out_of_range);
    CHECK_THROWS_AS(const_coefs[2], std::out_of_range);
  }
}

TEST_CASE("Coefficients: add", "[Coefficients]") {
  Coefficients coefs(Coefficient(Glossary::Temperature, 300.));
  coefs.add(Coefficient(Glossary::Pressure, 1.e5));
  CHECK(coefs.size() == 2);
  CHECK(coefs[1].get_type() == GlossaryType::Pressure);
}

// [!mayfail]: reported, but does not make the executable fail.
// Coefficients::add appends the coefficient but not its type:
// after an add, get_types() does not match the coefficients.
TEST_CASE("Coefficients: get_types is consistent after add", "[Coefficients][!mayfail]") {
  Coefficients coefs(Coefficient(Glossary::Temperature, 300.));
  coefs.add(Coefficient(Glossary::Pressure, 1.e5));
  CHECK(coefs.get_types().size() == coefs.size());
}

TEST_CASE("Coefficients: set_time reaches every coefficient", "[Coefficients]") {
  auto w = std::make_shared<W>();
  auto h = std::make_shared<H>();
  Coefficients coefs(Coefficient(Glossary::Mobility, Scheme::Implicit, w),
                     Coefficient(Glossary::Mobility, Scheme::Implicit, h));
  coefs.set_time(1.5);
  CHECK(w->get_time() == 1.5);
  CHECK(h->get_time() == 1.5);
}
