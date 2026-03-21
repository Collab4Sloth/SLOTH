/**
 * @file CalphadOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for Calphad problems
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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
#include <limits>
#include <string>
#include <vector>

#include "Utils/Utils.hpp"

namespace CalphadDefaultConstant {

/// @brief Threshold used to identify component to exclude from the initial
/// system
const auto residual_threshold = 1.e-10;
/// @brief Order for sorting nodes as function of the temperature
const std::string temperature_sort_method = "Descending";
/// @brief Order for sorting nodes as function of the pressure
const std::string pressure_sort_method = "No";

/// @brief Element removed when moles fractions initialization is considered
const std::string element_removed_from_ic = "Undefined";

/// @brief Flag to indicate if the heat capacity must be calculated
const bool compute_heat_capacity = true;

const int error_max = 3;

}  // namespace CalphadDefaultConstant

/**
 * @brief Available outputs for Calphad problems
 */
struct calphad_outputs {
  enum value { mu, dmu, g, gm, h, hm, x, xp, xph, y, dgm, cp, mob, nucleus, error };
  static value from(const std::string &);
};

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

struct pressure_sort_method {
  enum value { ascending, descending, no };
  static value from(const std::string &);
};

struct KKS_nucleation_strategy {
  enum value { liquid_fraction, given_melting_temperature };
  static value from(const std::string &);
};
