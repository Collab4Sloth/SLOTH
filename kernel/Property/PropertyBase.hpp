
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

  virtual void check_variables_consistency() = 0;
  virtual void get_property(std::vector<std::unique_ptr<mfem::Vector>> &vect_unk,
                            const std::vector<std::vector<std::string>> &unks_info,
                            const std::vector<mfem::ParGridFunction> &vect_aux_gf,
                            const std::vector<std::vector<std::string>> &vect_aux_infos) = 0;

 public:
  explicit PropertyBase(const Parameters &params);
  void compute(std::vector<std::unique_ptr<mfem::Vector>> &vect_unk,
               const std::vector<std::vector<std::string>> &unks_info,
               const std::vector<mfem::ParGridFunction> &vect_aux_gf,
               const std::vector<std::vector<std::string>> &vect_aux_infos);

  virtual ~PropertyBase();
};

/**
 * @brief Construct a new PropertyBase::PropertyBase object
 *
 */
PropertyBase::PropertyBase(const Parameters &params) : params_(params) {
  this->check_variables_consistency();
}

/**
 * @brief Compute the variables (properties) as functions of auxialiary variables
 *
 * @param vect_unk
 * @param unks_info
 * @param vect_aux_gf
 * @param vect_aux_infos
 */
void PropertyBase::compute(std::vector<std::unique_ptr<mfem::Vector>> &vect_unk,
                           const std::vector<std::vector<std::string>> &unks_info,
                           const std::vector<mfem::ParGridFunction> &vect_aux_gf,
                           const std::vector<std::vector<std::string>> &vect_aux_infos) {
  this->get_property(vect_unk, unks_info, vect_aux_gf, vect_aux_infos);
}

/**
 * @brief Destroy the PropertyBase::PropertyBase object
 *
 */
PropertyBase::~PropertyBase() {}
