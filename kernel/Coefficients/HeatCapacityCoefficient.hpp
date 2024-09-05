/**
 * @file HeatCapacityCoefficient.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class dedicated to HeatCapacity coefficient
 * @version 0.1
 * @date 2024-09-03
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "Coefficients/HeatCapacityFunctions.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

//--------------------------
//--------------------------
template <int ORDER, HeatCapacity NAME>
class HeatCapacityCoefficient : public mfem::Coefficient {
 private:
  HeatCapacityFunctions<ORDER, NAME> property_;
  FType property_function_;
  mfem::ParGridFunction *gfu_;
  double dble_gfu_{std::numeric_limits<double>::max()};

 public:
  HeatCapacityCoefficient(mfem::ParGridFunction *gfu, const Parameters &params);

  HeatCapacityCoefficient(const double gfu, const Parameters &params);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~HeatCapacityCoefficient();
};

/**
 * @brief Construct a new HeatCapacityCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param gfu
 * @param args
 */
template <int ORDER, HeatCapacity NAME>
HeatCapacityCoefficient<ORDER, NAME>::HeatCapacityCoefficient(mfem::ParGridFunction *gfu,
                                                              const Parameters &params)
    : gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

template <int ORDER, HeatCapacity NAME>
HeatCapacityCoefficient<ORDER, NAME>::HeatCapacityCoefficient(const double gfu,
                                                              const Parameters &params)
    : gfu_(nullptr), dble_gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

/**
 * @brief Return the value of the HeatCapacity coefficient at integration point
 *
 * @tparam ORDER
 * @tparam NAME
 * @param T
 * @param ip
 */
template <int ORDER, HeatCapacity NAME>
double HeatCapacityCoefficient<ORDER, NAME>::Eval(mfem::ElementTransformation &T,
                                                  const mfem::IntegrationPoint &ip) {
  auto var_at_ip = this->dble_gfu_;
  if (this->gfu_) {
    var_at_ip = this->gfu_->GetValue(T.ElementNo, ip);
  }

  return this->property_function_(var_at_ip);
}

/**
 * @brief Destroy the HeatCapacityCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, HeatCapacity NAME>
HeatCapacityCoefficient<ORDER, NAME>::~HeatCapacityCoefficient() {}
