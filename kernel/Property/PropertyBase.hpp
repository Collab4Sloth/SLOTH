
/**
 * @file PropertyBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for Property objects
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
#include <string>
#include <tuple>
#include <vector>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"

#pragma once

class PropertyBase {
 protected:
  // Parameters of the property
  const Parameters &params_;

  // Flag used to avoid verification at each time-step
  bool is_checked_{false};

  /**
   * @brief Check the consistency of the inputs and outputs of the property problem.
   *
   * @param output_system The outputs property problem (primary variables).
   * @param input_system  The inputs of the property problem (auxiliary variables).
   */
  virtual void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
          &output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) = 0;

  /**
   * @brief Get the values of the property
   *
   * The values of the property are stored in the outputs (primary variables of the property
   * problem). They are calculated from the inputs (auxiliary variables of the property problem).
   *
   * @param output_system The outputs property problem (primary variables).
   * @param input_system  The inputs of the property problem (auxiliary variables).
   */
  virtual void get_property(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
          &output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) = 0;

 public:
  explicit PropertyBase(const Parameters &params);
  void compute(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
          &output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system);

  virtual ~PropertyBase() = default;
};

/**
 * @brief Construct a new PropertyBase::PropertyBase object
 *
 */
PropertyBase::PropertyBase(const Parameters &params) : params_(params) {}

/**
 * @brief Compute the variables (properties) as functions of auxiliary variables
 *
 * @param output_system The outputs property problem (primary variables).
 * @param input_system  The inputs of the property problem (auxiliary variables).
 */
void PropertyBase::compute(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
        &output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  if (!this->is_checked_) {
    this->check_variables_consistency(output_system, input_system);
    this->is_checked_ = true;
  }

  this->get_property(output_system, input_system);
}
