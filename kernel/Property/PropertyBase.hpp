
/**
 * @file PropertyBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Base class for Property objects
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
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
  const Parameters &params_;
  bool is_checked_{false};
  virtual void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
          &output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) = 0;
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
 * @brief Compute the variables (properties) as functions of auxialiary variables
 *
 * @param vect_unk
 * @param unks_info
 * @param vect_aux_gf
 * @param vect_aux_infos
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
