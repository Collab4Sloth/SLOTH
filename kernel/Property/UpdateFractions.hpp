
/**
 * @file UpdateFractions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to update molar fractions
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
#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Property/PropertyBase.hpp"

#pragma once

/**
 * @brief Update molar fractions
 *
 *
 */
class UpdateFractions : public PropertyBase {
 protected:
  double factor_;
  std::vector<mfem::Vector> x_gf_;
  mfem::Vector x2up_gf_;

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;
  void get_property(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;

 public:
  explicit UpdateFractions(const Parameters& params);

  ~UpdateFractions();
};

/**
 * @brief Construct a new UpdateFractions::UpdateFractions object
 *
 */
UpdateFractions::UpdateFractions(const Parameters& params) : PropertyBase(params) {
  this->factor_ = this->params_.template get_param_value<double>("update_factor");
}

/**
 * @brief Check the consistency of the inputs and outputs of the property problem.
 *
 * @param output_system The outputs property problem (primary variables).
 * @param input_system  The inputs of the property problem (auxiliary variables).
 */
void UpdateFractions::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
        output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  // Check variables

  bool has_x = false;

  // Check auxiliary variables
  for (const auto& [aux_infos, aux_gf] : input_system) {
    MFEM_VERIFY(!aux_infos.empty(),
                "UpdateFractions::check_variables_consistency: error while "
                "getting auxiliary variables. Additional informations must be given for each "
                "auxiliary variables.");
    const int vsize = aux_infos.size();
    const std::string lower = toLowerCase(aux_infos.back());
    const std::string_view symbol(lower);

    if (symbol == "x") {
      has_x = true;
      MFEM_VERIFY(vsize == 2,
                  "UpdateFractions::check_variables_consistency: error while "
                  "getting molar fractions.Expected [element, 'x']");
    }
  }

  MFEM_VERIFY(has_x,
              "UpdateFractions::check_variables_consistency: error while "
              "getting x auxiliary variables");
}

/**
 * @brief Get the values of the property
 *
 * The values of the property are stored in the outputs (primary variables of the property
 * problem). They are calculated from the inputs (auxiliary variables of the property problem).
 *
 * @param output_system The outputs property problem (primary variables).
 * @param input_system  The inputs of the property problem (auxiliary variables).
 */
void UpdateFractions::get_property(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
        output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  x_gf_.clear();
  for (const auto& [aux_infos, aux_gf] : input_system) {
    if (!aux_infos.empty()) {
      const int size = aux_infos.size();
      const std::string_view type_info_view(aux_infos.back());
      if (!type_info_view.empty() && type_info_view == "x") {
        const std::string elem_info(aux_infos[size - 2]);
        this->x_gf_.emplace_back(std::move(aux_gf));
      }
    }
  }

  ///
  int j = 0;
  for (auto& [output_infos, output_value] : output_system) {
    mfem::Vector vv(this->x_gf_[0].Size());
    for (unsigned int k = 0; k < this->x_gf_[0].Size(); k++) {
      vv(k) = this->x_gf_[0][k] + this->factor_ * this->x_gf_[1][k];
    }
    output_value.get() = vv;
    j++;
  }
}

/**
 * @brief Destroy the UpdateFractions::UpdateFractions object
 *
 */
UpdateFractions::~UpdateFractions() {}
