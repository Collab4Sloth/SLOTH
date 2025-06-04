/**
 * @file CalphadOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for Calphad problems
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <limits>
#include <string>
#include <vector>

#include "Utils/Utils.hpp"

#pragma once

namespace CalphadDefaultConstant {

/// @brief Threshold used to identify component to exclude from the initial
/// system
const auto residual_threshold = 1.e-10;
/// @brief Order for sorting nodes as function of the temperature
const std::string temperature_sort_method = "Descending";
/// @brief Order for sorting nodes as function of the pressure
const std::string pressure_sort_method = "No";

/// @brief Element removed when moles fractions initialization is considered
const std::string element_removed_from_ic = "Undefined ";

/// @brief Flag to indicate if the heat capacity must be calculated
const bool compute_heat_capacity = true;

const int error_max = 3;

}  // namespace CalphadDefaultConstant

/**
 * @brief Available outputs for Calphad problems
 */
struct calphad_outputs {
  enum value { mu, dmu, g, gm, h, hm, x, xp, xph, y, dgm, cp, mob, nucleus };
  static value from(const std::string &);
};

calphad_outputs::value calphad_outputs::from(const std::string &v) {
  static PhaseFieldPrivate::mmap<calphad_outputs::value> m{
      {"mu", calphad_outputs::mu},   {"dmu", calphad_outputs::dmu},
      {"x", calphad_outputs::x},     {"xp", calphad_outputs::xp},
      {"xph", calphad_outputs::xph}, {"y", calphad_outputs::y},
      {"g", calphad_outputs::g},     {"gm", calphad_outputs::gm},
      {"h", calphad_outputs::h},     {"hm", calphad_outputs::hm},
      {"dgm", calphad_outputs::dgm}, {"cp", calphad_outputs::cp},
      {"mob", calphad_outputs::mob}, {"nucleus", calphad_outputs::nucleus}};
  return m.find("calphad_outputs", v);
}

struct calphad_nodes_strategies {
  enum value { temperature_sort_method, pressure_sort_method };
  static value from(const std::string &v) {
    static PhaseFieldPrivate::mmap<calphad_nodes_strategies::value> m{
        {"temperature_sort_method", calphad_nodes_strategies::temperature_sort_method},
        {"pressure_sort_method", calphad_nodes_strategies::pressure_sort_method}};
    return m.find("calphad_nodes_strategies", v);
  }
};

struct temperature_sort_method {
  enum value { ascending, descending, no };
  static value from(const std::string &);
};
temperature_sort_method::value temperature_sort_method::from(const std::string &v) {
  static PhaseFieldPrivate::mmap<temperature_sort_method::value> m{
      {"Ascending", temperature_sort_method::ascending},
      {"Descending", temperature_sort_method::descending},
      {"No", temperature_sort_method::no}};
  return m.find("temperature_sort_method", v);
}
struct pressure_sort_method {
  enum value { ascending, descending, no };
  static value from(const std::string &);
};
pressure_sort_method::value pressure_sort_method::from(const std::string &v) {
  static PhaseFieldPrivate::mmap<pressure_sort_method::value> m{
      {"Ascending", pressure_sort_method::ascending},
      {"Descending", pressure_sort_method::descending},
      {"No", pressure_sort_method::no}};
  return m.find("temperature_sort_method", v);
}
