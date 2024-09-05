/**
 * @file DiffusionCoefficient.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class dedicated to Diffusion coefficient
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

#include "Coefficients/DiffusionFunctions.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

//--------------------------
//--------------------------
template <int ORDER, Diffusion NAME>
class DiffusionCoefficient : public mfem::Coefficient {
 private:
  DiffusionFunctions<ORDER, NAME> property_;
  FType property_function_;
  mfem::ParGridFunction *gfu_;
  double dble_gfu_{std::numeric_limits<double>::max()};

 public:
  DiffusionCoefficient(mfem::ParGridFunction *gfu, const Parameters &params);

  DiffusionCoefficient(const double gfu, const Parameters &params);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~DiffusionCoefficient();
};

/**
 * @brief Construct a new DiffusionCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param gfu
 * @param args
 */
template <int ORDER, Diffusion NAME>
DiffusionCoefficient<ORDER, NAME>::DiffusionCoefficient(mfem::ParGridFunction *gfu,
                                                        const Parameters &params)
    : gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

template <int ORDER, Diffusion NAME>
DiffusionCoefficient<ORDER, NAME>::DiffusionCoefficient(const double gfu, const Parameters &params)
    : gfu_(nullptr), dble_gfu_(gfu) {
  this->property_function_ = this->property_.getFunction(params);
}

/**
 * @brief Return the value of the Diffusion coefficient at integration point
 *
 * @tparam ORDER
 * @tparam NAME
 * @param T
 * @param ip
 */
template <int ORDER, Diffusion NAME>
double DiffusionCoefficient<ORDER, NAME>::Eval(mfem::ElementTransformation &T,
                                               const mfem::IntegrationPoint &ip) {
  auto var_at_ip = this->dble_gfu_;
  if (this->gfu_) {
    var_at_ip = this->gfu_->GetValue(T.ElementNo, ip);
  }

  return this->property_function_(var_at_ip);
}

/**
 * @brief Destroy the DiffusionCoefficient object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Diffusion NAME>
DiffusionCoefficient<ORDER, NAME>::~DiffusionCoefficient() {}
