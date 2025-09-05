/**
 * @file EnergyCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Energy coefficients
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
  ~GradientEnergyCoefficient();
};

/**
 * @brief Construct a new GradientEnergyCoefficient::GradientEnergyCoefficient
 * object
 *
 * @param gfu
 * @param lambda
 */
GradientEnergyCoefficient::GradientEnergyCoefficient(mfem::ParGridFunction *gfu,
                                                     const double &lambda)
    : gfu_(gfu), lambda_(lambda) {}

/**
 * @brief Evaluate the GradientEnergyCoefficient at integration point
 *
 * @param T
 * @param ip
 * @return double
 */
double GradientEnergyCoefficient::Eval(mfem::ElementTransformation &T,
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
GradientEnergyCoefficient::~GradientEnergyCoefficient() {}

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
  ~HomogeneousEnergyCoefficient();
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

/**
 * @brief Destroy the HomogeneousEnergyCoefficient< ENERGY>::HomogeneousEnergyCoefficient
 * object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 */
template <ThermodynamicsPotentials ENERGY>
HomogeneousEnergyCoefficient<ENERGY>::~HomogeneousEnergyCoefficient() {}
