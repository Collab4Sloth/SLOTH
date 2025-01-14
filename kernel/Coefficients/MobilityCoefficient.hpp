/**
 * @file MobilityCoefficient.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class dedicated to Mobility coefficient
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

#include "Coefficients/MobilityFunctions.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

//--------------------------
//--------------------------
template <int ORDER, Mobility NAME>
class MobilityCoefficient : public mfem::Coefficient {
 private:
  MobilityFunctions<ORDER, NAME> property_;
  FType property_function_;
  mfem::ParGridFunction *gfu_;
  double dble_gfu_{std::numeric_limits<double>::max()};

 public:
  MobilityCoefficient(mfem::ParGridFunction *gfu, const Parameters &params);

  MobilityCoefficient(const double gfu, const Parameters &params);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~MobilityCoefficient();
};

/**
 * @brief Construct a new MobilityCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param gfu
 * @param args
 */
template <int ORDER, Mobility NAME>
MobilityCoefficient<ORDER, NAME>::MobilityCoefficient(mfem::ParGridFunction *gfu,
                                                      const Parameters &params)
    : gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

template <int ORDER, Mobility NAME>
MobilityCoefficient<ORDER, NAME>::MobilityCoefficient(const double gfu, const Parameters &params)
    : gfu_(nullptr), dble_gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

/**
 * @brief Return the value of the Mobility coefficient at integration point
 *
 * @tparam ORDER
 * @tparam NAME
 * @param T
 * @param ip
 */
template <int ORDER, Mobility NAME>
double MobilityCoefficient<ORDER, NAME>::Eval(mfem::ElementTransformation &T,
                                              const mfem::IntegrationPoint &ip) {
  auto var_at_ip = this->dble_gfu_;
  if (this->gfu_) {
    var_at_ip = this->gfu_->GetValue(T.ElementNo, ip);
  }

  return this->property_function_(var_at_ip);
}

/**
 * @brief Destroy the MobilityCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Mobility NAME>
MobilityCoefficient<ORDER, NAME>::~MobilityCoefficient() {}
