/**
 * @file OmegaCoefficient.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class dedicated to Omega coefficient
 * @version 0.1
 * @date 2025-01-22
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "Coefficients/OmegaFunctions.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

//--------------------------
//--------------------------
template <int ORDER, Omega NAME>
class OmegaCoefficient : public mfem::Coefficient {
 private:
  OmegaFunctions<ORDER, NAME> property_;
  FType property_function_;
  mfem::ParGridFunction *gfu_;
  double dble_gfu_{std::numeric_limits<double>::max()};

 public:
  OmegaCoefficient(mfem::ParGridFunction *gfu, const Parameters &params);

  OmegaCoefficient(const double gfu, const Parameters &params);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~OmegaCoefficient();
};

/**
 * @brief Construct a new OmegaCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param gfu
 * @param args
 */
template <int ORDER, Omega NAME>
OmegaCoefficient<ORDER, NAME>::OmegaCoefficient(mfem::ParGridFunction *gfu,
                                                const Parameters &params)
    : gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

template <int ORDER, Omega NAME>
OmegaCoefficient<ORDER, NAME>::OmegaCoefficient(const double gfu, const Parameters &params)
    : gfu_(nullptr), dble_gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

/**
 * @brief Return the value of the Omega coefficient at integration point
 *
 * @tparam ORDER
 * @tparam NAME
 * @param T
 * @param ip
 */
template <int ORDER, Omega NAME>
double OmegaCoefficient<ORDER, NAME>::Eval(mfem::ElementTransformation &T,
                                           const mfem::IntegrationPoint &ip) {
  auto var_at_ip = this->dble_gfu_;
  if (this->gfu_) {
    var_at_ip = this->gfu_->GetValue(T.ElementNo, ip);
  }

  return this->property_function_(var_at_ip);
}

/**
 * @brief Destroy the OmegaCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Omega NAME>
OmegaCoefficient<ORDER, NAME>::~OmegaCoefficient() {}
