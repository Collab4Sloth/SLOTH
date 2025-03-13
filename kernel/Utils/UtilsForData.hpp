/**
 * @file UtilsForData.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Usefull methods for data
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */
#include <unistd.h>

#include <algorithm>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

using FuncType = std::function<double(const double&, const double&)>;
using FType = std::function<double(const double&)>;
using triplet = std::tuple<std::string, double, std::string>;
using vtriplet = std::vector<triplet>;
using SpecializedValue = std::pair<std::string, double>;

static auto stringfindInVectorOfString = [](const std::vector<std::string> v, const std::string w) {
  return find(v.begin(), v.end(), w) != v.end();
};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Levels of verbosity
 *
 */
enum class Verbosity { Quiet, Normal, Verbose, Debug, Error };

/**
 * @brief The current level of verbosity
 *
 */
static Verbosity verbosityLevel = Verbosity::Quiet;

/**
 * @brief Set the level of verbosity
 *
 * @param verbosity
 */
static void setVerbosity(Verbosity verbosity) { verbosityLevel = verbosity; }

/**
 * @brief Method used to capture the flux coming from an external method and print it depending on
 * the verbosity level
 *
 * @tparam Func
 * @tparam Args
 * @param verbosity
 * @param func
 * @param args
 */
template <typename Func, typename... Args>
static void external_call(Verbosity verbosity, Func func, Args&&... args) {
  if (verbosity <= Verbosity::Verbose) {
    // Save fluxes
    int stdout_backup = dup(fileno(stdout));
    int stderr_backup = dup(fileno(stderr));

    // Redirection towards /dev/null to suppress output
    FILE* stdout_null = freopen("/dev/null", "w", stdout);
    FILE* stderr_null = freopen("/dev/null", "w", stderr);

    if (!stdout_null || !stderr_null) {
      // Handle error: restore original stdout and stderr
      if (stdout_null) {
        fclose(stdout_null);
      }
      if (stderr_null) {
        fclose(stderr_null);
      }
      // Optionally, you can log the error or throw an exception here
      return;
    }
    // Call of the function
    func(std::forward<Args>(args)...);

    // Restore fluxes
    fflush(stdout);
    fflush(stderr);
    dup2(stdout_backup, fileno(stdout));
    dup2(stderr_backup, fileno(stderr));

    // Close temporary fluxes
    close(stdout_backup);
    close(stderr_backup);
  } else {
    // Standard call of the function
    func(std::forward<Args>(args)...);
  }
}

/**
 * @brief Class used to print information depending on the level of verbosity
 *
 */
class SlothInfo {
 public:
  /**
   * @brief Print information on Debug mode
   *
   * @tparam Args
   * @param args
   */
  template <typename... Args>
  static void debug(Args... args) {
    int rank = mfem::Mpi::WorldRank();

    if (Verbosity::Debug <= verbosityLevel && rank == 0) {
      (std::cout << ... << args) << "\n";
    }
  }
  /**
   * @brief Print information on Error mode
   *
   * @tparam Args
   * @param args
   */
  template <typename... Args>
  static void error(Args... args) {
    int rank = mfem::Mpi::WorldRank();
    if (Verbosity::Error <= verbosityLevel && rank == 0) {
      (std::cout << ... << args) << "\n";
    }
  }

  /**
   * @brief Print information on Normal mode
   *
   * @tparam Args
   * @param args
   */
  template <typename... Args>
  static void verbose(Args... args) {
    int rank = mfem::Mpi::WorldRank();
    if (Verbosity::Verbose <= verbosityLevel && rank == 0) {
      (std::cout << ... << args) << "\n";
    }
  }

  /**
   * @brief Print information on Normal mode
   *
   * @tparam Args
   * @param args
   */
  template <typename... Args>
  static void print(Args... args) {
    int rank = mfem::Mpi::WorldRank();
    if (Verbosity::Normal <= verbosityLevel && rank == 0) {
      (std::cout << ... << args) << "\n";
    }
  }
};
