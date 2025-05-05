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
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

using FuncType = std::function<double(const double&, const double&)>;
using FType = std::function<double(const double&)>;
using vTupleStringInt = std::vector<std::tuple<std::string, int>>;
using vTupleStringString = std::vector<std::tuple<std::string, std::string>>;
using vString = std::vector<std::string>;

using MapStringDouble = std::map<std::string, double>;
using vTuple2StringDouble = std::vector<std::tuple<std::string, std::string, double>>;
using Map2String2Double =
    std::map<std::tuple<std::string, std::string>, std::tuple<double, double>>;
using SpecializedValue = std::pair<std::string, double>;

using triplet = std::tuple<std::string, double, std::string>;
using vtriplet = std::vector<triplet>;

/**
 * @brief Search string in a vector of string
 *
 */
static auto stringfindInVectorOfString = [](const std::vector<std::string> v, const std::string w) {
  return find(v.begin(), v.end(), w) != v.end();
};

/**
 * @brief Put string in uppercase
 *
 * @param str
 * @return std::string
 */
static std::string toUpperCase(const std::string& str) {
  std::string result = str;
  std::transform(result.begin(), result.end(), result.begin(), ::toupper);
  return result;
}

/**
 * @brief Put string in lowercase
 *
 * @param str
 * @return std::string
 */
static std::string toLowerCase(const std::string& str) {
  std::string result = str;
  std::transform(result.begin(), result.end(), result.begin(), ::tolower);
  return result;
}

/**
 * @brief Trim leading and trailing whitespace
 *
 * @param str
 * @return std::string
 */
static std::string trim(const std::string& str) {
  auto start = std::ranges::find_if(
      str, [](char c) { return !std::isspace(static_cast<unsigned char>(c)); });
  auto end = std::ranges::find_if(str.rbegin(), str.rend(), [](char c) {
               return !std::isspace(static_cast<unsigned char>(c));
             }).base();

  // Check if the range is empty (if start is greater than or equal to end)
  return (start >= end) ? "" : std::string(start, end);
}

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

// /**
//  * @brief PRint Verbosity level
//  *
//  */
// static void printVerbosityLevel() {
//   switch (verbosityLevel) {
//     case Verbosity::Quiet:
//       std::cout << "Verbosity Level: Quiet" << std::endl;
//       break;
//     case Verbosity::Normal:
//       std::cout << "Verbosity Level: Normal" << std::endl;
//       break;
//     case Verbosity::Verbose:
//       std::cout << "Verbosity Level: Verbose" << std::endl;
//       break;
//     case Verbosity::Debug:
//       std::cout << "Verbosity Level: Debug" << std::endl;
//       break;
//     case Verbosity::Error:
//       std::cout << "Verbosity Level: Error" << std::endl;
//       break;
//   }
// }

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
    if (stdout_backup == -1 || stderr_backup == -1) {
      // Handle error if dup() fails (file descriptor cannot be duplicated)
      return;
    }

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
 * @brief Check existence of a file
 *
 * @param file
 */
static void check_file_existence(const std::string& file) {
  std::filesystem::path filePath = file;
  // Check if the file exists
  if (!std::filesystem::exists(filePath)) {
    const std::string& error_msg = "file " + file + " doesn't exist. Please check your data.";
    mfem::mfem_error(error_msg.c_str());
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
    if (Verbosity::Debug <= verbosityLevel) {
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
    if (Verbosity::Error <= verbosityLevel) {
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
    if (Verbosity::Verbose <= verbosityLevel) {
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
    if (Verbosity::Normal <= verbosityLevel) {
      (std::cout << ... << args) << "\n";
    }
  }
};
