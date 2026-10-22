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
