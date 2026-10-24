/**
 * @file test_Parameter.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for Parameter (one named value)
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
#include <variant>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Parameters/Parameter.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE("Parameter: name and description", "[Parameters]") {
  SECTION("without description: empty description") {
    const Parameter p("alpha", 1.5);
    CHECK(p.get_name() == "alpha");
    CHECK(p.get_description().empty());
  }
  SECTION("with description") {
    const Parameter p("alpha", 1.5, "interface mobility");
    CHECK(p.get_name() == "alpha");
    CHECK(p.get_description() == "interface mobility");
  }
}

TEST_CASE("Parameter: the stored type is the type of the value", "[Parameters]") {
  SECTION("int") {
    const Parameter p("n", 3);
    REQUIRE(std::holds_alternative<int>(p.get_value()));
    CHECK(std::get<int>(p.get_value()) == 3);
  }
  SECTION("double") {
    const Parameter p("x", 2.5);
    REQUIRE(std::holds_alternative<double>(p.get_value()));
    CHECK_THAT(std::get<double>(p.get_value()), WithinRel(2.5, 1.e-12));
  }
  SECTION("bool") {
    const Parameter p("flag", true);
    REQUIRE(std::holds_alternative<bool>(p.get_value()));
    CHECK(std::get<bool>(p.get_value()));
  }
  SECTION("std::string") {
    const Parameter p("name", std::string("phi"));
    REQUIRE(std::holds_alternative<std::string>(p.get_value()));
    CHECK(std::get<std::string>(p.get_value()) == "phi");
  }
  SECTION("a string literal is stored as std::string, not as bool") {
    // C++20 (P0608): std::variant does not convert "phi" (const char*) to bool.
    const Parameter p("name", "phi");
    CHECK(std::holds_alternative<std::string>(p.get_value()));
  }
  SECTION("containers") {
    const Parameter pv("v", vDouble{1., 2., 3.});
    REQUIRE(std::holds_alternative<vDouble>(pv.get_value()));
    CHECK(std::get<vDouble>(pv.get_value()).size() == 3);

    const Parameter pm("m", MapStringDouble{{"U", 0.1}, {"O", 0.9}});
    REQUIRE(std::holds_alternative<MapStringDouble>(pm.get_value()));
    CHECK_THAT(std::get<MapStringDouble>(pm.get_value()).at("O"), WithinRel(0.9, 1.e-12));
  }
}

TEST_CASE("Parameter: print", "[Parameters]") {
  // SLOTH verbosity is Quiet by default: nothing is printed. These checks
  // run every branch of print() (one per type) and check that none fails.
  CHECK_NOTHROW(Parameter("i", 1).print());
  CHECK_NOTHROW(Parameter("d", 1.).print());
  CHECK_NOTHROW(Parameter("s", std::string("a")).print());
  CHECK_NOTHROW(Parameter("b", false).print());
  CHECK_NOTHROW(Parameter("vs", vString{"a", "b"}).print());
  CHECK_NOTHROW(Parameter("vi", vInt{1, 2}).print());
  CHECK_NOTHROW(Parameter("vd", vDouble{1., 2.}).print());
  CHECK_NOTHROW(Parameter("vtss", vTupleStringString{{"a", "b"}}).print());
  CHECK_NOTHROW(Parameter("vtsi", vTupleStringInt{{"a", 1}}).print());
  CHECK_NOTHROW(Parameter("msd", MapStringDouble{{"a", 1.}}).print());
}
