/**
 * @file sloth_unit.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Helpers shared by SLOTH unit tests
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
#pragma once

#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mfem.hpp"

/**
 * @brief Checks that an expression ends in mfem::mfem_error (SLOTH errors, MFEM_ABORT...).
 *
 * SLOTH reports errors with mfem::mfem_error. By default, mfem_error stops the
 * whole program (MPI_Abort). When MFEM is built with MFEM_USE_EXCEPTIONS, the
 * common main makes it throw mfem::ErrorException instead, which Catch2 can
 * catch and check without stopping the other tests.
 *
 * Without MFEM_USE_EXCEPTIONS, these checks are not run: a warning is printed
 * and the rest of the test case goes on.
 *
 *   SLOTH_CHECK_MFEM_ERROR(expr)              expr ends in mfem_error
 *   SLOTH_CHECK_MFEM_ERROR_WITH(expr, "txt")  same, and the message contains "txt"
 */
#ifdef MFEM_USE_EXCEPTIONS

#define SLOTH_CHECK_MFEM_ERROR(...) CHECK_THROWS_AS(__VA_ARGS__, mfem::ErrorException)

// expr is evaluated once: exception type and message are checked together
#define SLOTH_CHECK_MFEM_ERROR_WITH(expr, text)   \
  CHECK_THROWS_MATCHES(expr, mfem::ErrorException, \
                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring(text)))

#else

#define SLOTH_CHECK_MFEM_ERROR(...) \
  WARN("not checked (MFEM built without MFEM_USE_EXCEPTIONS): " #__VA_ARGS__)

#define SLOTH_CHECK_MFEM_ERROR_WITH(expr, text) \
  WARN("not checked (MFEM built without MFEM_USE_EXCEPTIONS): " #expr)

#endif

namespace sloth_unit {

/// @brief Builds an mfem::Vector from a list of values: point({0.1, 0.2})
inline mfem::Vector point(std::initializer_list<double> values) {
  mfem::Vector v(static_cast<int>(values.size()));
  int i = 0;
  for (const double x : values) {
    v[i++] = x;
  }
  return v;
}

/**
 * @brief Captures what is written on std::cout while the object is alive.
 *
 *   sloth_unit::CoutCapture out;
 *   SlothInfo::print("hello");
 *   CHECK(out.str() == "hello\n");
 *
 * std::cout is restored when the object is destroyed, even if a check fails.
 */
class CoutCapture {
 public:
  CoutCapture() : old_(std::cout.rdbuf(buffer_.rdbuf())) {}
  ~CoutCapture() { std::cout.rdbuf(old_); }
  CoutCapture(const CoutCapture&) = delete;
  CoutCapture& operator=(const CoutCapture&) = delete;

  std::string str() const { return buffer_.str(); }

 private:
  std::ostringstream buffer_;
  std::streambuf* old_;
};

/**
 * @brief File written in the system temporary directory for one test, and
 * removed when the object is destroyed (even if a check fails).
 *
 *   sloth_unit::TemporaryFile csv("field.csv", "t,x,F\n0,0,1\n...");
 *   FileCoefficient c(csv.path());
 *
 * The name is made unique by the process id and the MPI rank, so that tests
 * run in parallel (ctest -j, several ranks) do not write the same file.
 */
class TemporaryFile {
 public:
  TemporaryFile(const std::string& name, const std::string& content)
      : path_(std::filesystem::temp_directory_path() /
              ("sloth_unit_" + std::to_string(::getpid()) + "_" +
               std::to_string(mfem::Mpi::WorldRank()) + "_" + name)) {
    std::ofstream out(path_);
    out << content;
  }
  ~TemporaryFile() {
    std::error_code ec;
    std::filesystem::remove(path_, ec);  // never throws in a destructor
  }
  TemporaryFile(const TemporaryFile&) = delete;
  TemporaryFile& operator=(const TemporaryFile&) = delete;

  std::string path() const { return path_.string(); }

 private:
  std::filesystem::path path_;
};

}  // namespace sloth_unit
