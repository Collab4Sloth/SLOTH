/**
 * @file HomogeneousEnergyCoefficient.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-18
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <numeric>

#include "Coefficients/PhaseFieldPotentials.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

//--------------------------
//--------------------------
//--------------------------
class GradientEnergyCoefficient : public mfem::Coefficient {
 private:
  mfem::ParGridFunction *gfu_;
  double lambda_;

 public:
  GradientEnergyCoefficient(mfem::ParGridFunction *gfu, const double &lambda);
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~GradientEnergyCoefficient() = default;
};

/**
 * @brief Construct a new GradientEnergyCoefficient::GradientEnergyCoefficient
 * object
 *
 * @param gfu
 * @param lambda
 */
DEBILE_INLINE GradientEnergyCoefficient::GradientEnergyCoefficient(mfem::ParGridFunction *gfu,
                                                     const double &lambda)
    : gfu_(gfu), lambda_(lambda) {}

/**
 * @brief Evaluate the GradientEnergyCoefficient at integration point
 *
 * @param T
 * @param ip
 * @return double
 */
DEBILE_INLINE double GradientEnergyCoefficient::Eval(mfem::ElementTransformation &T,
                                       const mfem::IntegrationPoint &ip) {
  mfem::Vector gradu;
  this->gfu_->GetGradient(T, gradu);

  const auto &value = std::inner_product(gradu.begin(), gradu.end(), gradu.begin(), 0.);
  return this->lambda_ * value;
}

/**
 * @brief Destroy the GradientEnergyCoefficient::GradientEnergyCoefficient object
 *
 */
//GradientEnergyCoefficient::~GradientEnergyCoefficient() {}

//--------------------------
//--------------------------
//--------------------------

template <ThermodynamicsPotentials ENERGY>
class HomogeneousEnergyCoefficient : public mfem::Coefficient {
 private:
  PotentialFunctions<0, ThermodynamicsPotentialDiscretization::Implicit, ENERGY> energy_potential_;
  std::function<double(const double &)> energy_function_;
  mfem::ParGridFunction *gfu_;
  double omega_;

 public:
  HomogeneousEnergyCoefficient(mfem::ParGridFunction *gfu, const double &omega);
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~HomogeneousEnergyCoefficient() = default;
};

/**
 * @brief Construct a new HomogeneousEnergyCoefficient< ENERGY>::HomogeneousEnergyCoefficient
 * object
 *
 * @tparam ENERGY
 * @param omega_
 */
template <ThermodynamicsPotentials ENERGY>
HomogeneousEnergyCoefficient<ENERGY>::HomogeneousEnergyCoefficient(mfem::ParGridFunction *gfu,
                                                                   const double &omega)
    : gfu_(gfu), omega_(omega) {
  this->energy_function_ = this->energy_potential_.getPotentialFunction();
}

/**
 * @brief Evaluate the HomogeneousEnergyCoefficient at integration point
 *
 * @tparam ENERGY
 * @param T
 * @param ip
 * @return double
 */
template <ThermodynamicsPotentials ENERGY>
double HomogeneousEnergyCoefficient<ENERGY>::Eval(mfem::ElementTransformation &T,
                                                  const mfem::IntegrationPoint &ip) {
  const auto phi = this->gfu_->GetValue(T.ElementNo, ip);
  return this->omega_ * this->energy_function_(phi);
}
