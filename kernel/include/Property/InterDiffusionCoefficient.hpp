
/**
 * @file InterDiffusionCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to compute inter-diffusion coefficients
 * @version 0.1
 * @date 2025-09-05
 *
 * Copyright CEA (C) 2025
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
#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Property/PropertyBase.hpp"

/**
 * @brief Compute inter-diffusion coefficients
 *
 * J_i = Sum_j [ M_ij * grad mu_i ] - Sum_j [ M_ij * grad mu_j ]
 *
 * with M_ij = x_i * x_j * (M_i x_j + M_j x_i)
 *
 */
class InterDiffusionCoefficient : public PropertyBase {
 private:
  // List of the chemical elements
  std::set<std::string> list_components_;

  // Flag used to know if the formula must be extended in two-phase
  bool single_phase_ = true;

 protected:
  // Chemical element corresponding to the first element
  std::string first_component_;
  // Chemical element corresponding to the last element (the reference for diffusion chemical
  // potential)
  std::string last_component_;
  // The name of the primary phase
  std::string primary_phase_;
  // The name of the secondary phase (the phase formed during phase change)
  std::string secondary_phase_;
  // The phase-field variable
  std::vector<mfem::Vector> phi_gf_;
  // The molar fraction of the elements
  std::map<std::string, mfem::Vector> x_gf_;
  // The mobilities of the elements for the primary phase
  std::map<std::string, mfem::Vector> primary_mob_gf_;
  // The mobilities of the elements for the secondary phase
  std::map<std::string, mfem::Vector> secondary_mob_gf_;

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;
  void get_property(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;

 public:
  explicit InterDiffusionCoefficient(const Parameters& params);

  ~InterDiffusionCoefficient();
};
