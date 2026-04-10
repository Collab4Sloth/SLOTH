/**
 * @file MfemCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Built a MFEM Coefficient from a SLOTH Coefficient
 * @version 0.1
 * @date 2026-03-30
 *
 * @copyright CEA (C) 2026
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
#pragma once
#include <span>
#include <vector>

#include "Coefficients/Coefficients.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

class MfemCoefficient : public mfem::Coefficient {
 private:
  int id_ = 0;
  Coefficients sloth_coefficient_;
  const mfem::ParGridFunction& gf_;
  const std::vector<mfem::ParGridFunction>& vaux_gf_;

  mutable std::vector<double> vaux_gf_at_ip_;

 public:
  MfemCoefficient(int id, const Coefficients& coef, const mfem::ParGridFunction& vun,
                  const std::vector<mfem::ParGridFunction>& vauxn);
  double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;
  virtual ~MfemCoefficient() = default;
};

/**
 * @brief Construct a new Mfem Coefficient:: Mfem Coefficient object
 *
 * @param coef
 * @param vun
 * @param vauxn
 */
MfemCoefficient::MfemCoefficient(int id, const Coefficients& coef, const mfem::ParGridFunction& vun,
                                 const std::vector<mfem::ParGridFunction>& vauxn)
    : id_(id), sloth_coefficient_(coef), gf_(vun), vaux_gf_(vauxn) {}

/**
 * @brief Eval MfemCoefficient at the integration point
 *
 * @param T
 * @param ip
 * @return double
 */
double MfemCoefficient::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) {
  if (this->sloth_coefficient_[id_].is_scalar()) {
    return this->sloth_coefficient_[id_].compute();
  }

  this->vaux_gf_at_ip_.resize(this->vaux_gf_.size());
  for (size_t k = 0; k < this->vaux_gf_.size(); ++k) {
    this->vaux_gf_at_ip_[k] = this->vaux_gf_[k].GetValue(T, ip);
  }

  double var_at_ip = this->gf_.GetValue(T, ip);
  std::array<double, 1> var = {var_at_ip};
  return this->sloth_coefficient_[id_].compute(std::span<const double>(var),
                                               std::span<const double>(this->vaux_gf_at_ip_));
}
