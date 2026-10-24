/**
 * @file test_Utils.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for Utils (strings, verbosity, debug helpers)
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
#include <catch2/matchers/catch_matchers_string.hpp>

#include "Utils/UtilsForData.hpp"
#include "Utils/UtilsForDebug.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::ContainsSubstring;

namespace {
/**
 * @brief Sets the SLOTH verbosity for the duration of a test, then restores it.
 *
 * The verbosity is a global variable: a test that changes it without restoring
 * it would change the behaviour of the tests run after it.
 */
class VerbosityGuard {
 public:
  explicit VerbosityGuard(Verbosity v) : old_(verbosityLevel) { setVerbosity(v); }
  ~VerbosityGuard() { setVerbosity(old_); }
  VerbosityGuard(const VerbosityGuard&) = delete;
  VerbosityGuard& operator=(const VerbosityGuard&) = delete;

 private:
  Verbosity old_;
};
}  // namespace

// -----------------------------------------------------------------------------
// Strings
// -----------------------------------------------------------------------------
TEST_CASE("Utils: case conversion", "[Utils]") {
  CHECK(toUpperCase("AllenCahn_2D") == "ALLENCAHN_2D");
  CHECK(toLowerCase("AllenCahn_2D") == "allencahn_2d");
  CHECK(toUpperCase("").empty());
}

TEST_CASE("Utils: trim", "[Utils]") {
  CHECK(trim("  phi \t\n") == "phi");
  CHECK(trim("phi") == "phi");
  CHECK(trim("a b") == "a b");  // inner spaces are kept
  CHECK(trim("   ").empty());
  CHECK(trim("").empty());
}

TEST_CASE("Utils: replaceSpaces and replacePunctuation", "[Utils]") {
  CHECK(replaceSpaces("molar fraction\tU") == "molar_fraction_U");
  CHECK(replacePunctuation("x(U).O,2") == "x_U__O_2");
}

TEST_CASE("Utils: stringfindInVectorOfString", "[Utils]") {
  const std::vector<std::string> v{"phi", "mu", "T"};
  CHECK(stringfindInVectorOfString(v, "mu"));
  CHECK_FALSE(stringfindInVectorOfString(v, "MU"));  // case-sensitive
  CHECK_FALSE(stringfindInVectorOfString({}, "mu"));
}

// -----------------------------------------------------------------------------
// Verbosity and SlothInfo
// -----------------------------------------------------------------------------
TEST_CASE("Utils: verbosity levels", "[Utils]") {
  const VerbosityGuard guard(Verbosity::Verbose);
  CHECK(verbose_is(Verbosity::Verbose));
  CHECK(verbose_at_least(Verbosity::Quiet));
  CHECK(verbose_at_least(Verbosity::Normal));
  CHECK(verbose_at_least(Verbosity::Verbose));
  CHECK_FALSE(verbose_at_least(Verbosity::Debug));
}

TEST_CASE("Utils: SlothInfo prints according to the verbosity", "[Utils]") {
  SECTION("Quiet: nothing is printed") {
    const VerbosityGuard guard(Verbosity::Quiet);
    sloth_unit::CoutCapture out;
    SlothInfo::print("print");
    SlothInfo::verbose("verbose");
    SlothInfo::debug("debug");
    CHECK(out.str().empty());
  }
  SECTION("Normal: print only") {
    const VerbosityGuard guard(Verbosity::Normal);
    sloth_unit::CoutCapture out;
    SlothInfo::print("x = ", 2);
    SlothInfo::verbose("verbose");
    CHECK(out.str() == "x = 2\n");  // arguments are concatenated
  }
  SECTION("Debug: print, verbose and debug") {
    const VerbosityGuard guard(Verbosity::Debug);
    sloth_unit::CoutCapture out;
    SlothInfo::print("a");
    SlothInfo::verbose("b");
    SlothInfo::debug("c");
    CHECK(out.str() == "a\nb\nc\n");
  }
  SECTION("error messages are printed with Verbosity::Error only") {
    // Verbosity::Error is the LAST level of the enum: SlothInfo::error only
    // prints with Verbosity::Error, not with Normal or Debug.
    {
      const VerbosityGuard guard(Verbosity::Debug);
      sloth_unit::CoutCapture out;
      SlothInfo::error("err");
      CHECK(out.str().empty());
    }
    {
      const VerbosityGuard guard(Verbosity::Error);
      sloth_unit::CoutCapture out;
      SlothInfo::error("err");
      CHECK(out.str() == "err\n");
    }
  }
}

// -----------------------------------------------------------------------------
// Debug helpers
// -----------------------------------------------------------------------------
TEST_CASE("Utils: memory checkpoints", "[Utils]") {
  SECTION("make_memory_checkpoint: the memory footprint is positive") {
    const auto usage = UtilsForDebug::make_memory_checkpoint();
    CHECK(usage.ru_maxrss > 0);
  }
  SECTION("memory_checkpoint prints the message and a size in GB") {
    sloth_unit::CoutCapture out;
    UtilsForDebug::memory_checkpoint("after assembly");
    CHECK_THAT(out.str(), ContainsSubstring("after assembly"));
    CHECK_THAT(out.str(), ContainsSubstring("GB"));
  }
}

TEST_CASE("Utils: throw_if", "[Utils]") {
  CHECK_NOTHROW(throw_if(false, "method", "never raised"));
  SLOTH_CHECK_MFEM_ERROR(throw_if(true, "method", "raised"));
}
