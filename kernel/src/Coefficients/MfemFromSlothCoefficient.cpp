/**
 * @file MfemFromSlothCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Wraps a SLOTH Coefficient, possibly depending on other variables, as an mfem::Coefficient
 * @version 0.1
 * @date 2026-11-26
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
#include "Coefficients/MfemFromSlothCoefficient.hpp"

#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new MfemFromSlothCoefficient object
 *
 * @param coef SLOTH coefficient to wrap (copied)
 */
MfemFromSlothCoefficient::MfemFromSlothCoefficient(::Coefficient coef)
    : sloth_coefficient_(std::move(coef)) {}

/**
 * @brief True if the wrapped coefficient is a constant
 */
bool MfemFromSlothCoefficient::is_scalar() const { return this->sloth_coefficient_.is_scalar(); }

/**
 * @brief True if the wrapped coefficient is file-based (Scheme::File)
 */
bool MfemFromSlothCoefficient::is_file() const { return this->sloth_coefficient_.is_file(); }

/**
 * @brief True if the wrapped coefficient may depend on other variables, i.e. if it is
 *        neither a constant nor file-based
 */
bool MfemFromSlothCoefficient::requires_dependencies() const {
  return !this->is_scalar() && !this->is_file();
}

/**
 * @brief Set the grid functions whose values are forwarded to the wrapped coefficient
 *
 * @param dependencies Non-owning pointers, in the order expected by the coefficient
 */
void MfemFromSlothCoefficient::bind(std::vector<const mfem::ParGridFunction*> dependencies) {
  for (const auto* gf : dependencies) {
    MFEM_VERIFY(gf != nullptr, "MfemFromSlothCoefficient::bind(): null grid function");
  }
  this->dependencies_ = std::move(dependencies);
}

/**
 * @brief Forget the grid functions set by bind()
 */
void MfemFromSlothCoefficient::unbind() { this->dependencies_.clear(); }

/**
 * @brief Number of grid functions currently bound
 */
std::size_t MfemFromSlothCoefficient::get_dependencies_number() const {
  return this->dependencies_.size();
}

/**
 * @brief Set the time of both the mfem::Coefficient and the wrapped SLOTH coefficient
 *
 * @param t Time
 */
void MfemFromSlothCoefficient::SetTime(mfem::real_t t) {
  mfem::Coefficient::SetTime(t);
  if (!this->is_scalar()) {
    this->sloth_coefficient_.set_time(t);
  }
}

/**
 * @brief Evaluate the wrapped coefficient at an integration point
 *
 * @param T Element transformation
 * @param ip Integration point
 * @return mfem::real_t
 */
mfem::real_t MfemFromSlothCoefficient::Eval(mfem::ElementTransformation& T,
                                            const mfem::IntegrationPoint& ip) {
  if (this->is_scalar()) {
    return this->sloth_coefficient_.compute();
  }

  T.SetIntPoint(&ip);
  T.Transform(ip, this->x_);
  const std::span<const double> x(this->x_.GetData(), this->x_.Size());

  if (this->is_file()) {
    return this->sloth_coefficient_.compute(x, std::span<const double>());
  }

  if (this->sloth_coefficient_.uses_coordinates()) {
    this->sloth_coefficient_.set_coordinates(x);
  }
  this->values_at_ip_.resize(this->dependencies_.size());
  for (std::size_t k = 0; k < this->dependencies_.size(); ++k) {
    this->values_at_ip_[k] = this->dependencies_[k]->GetValue(T, ip);
  }
  return this->sloth_coefficient_.compute(std::span<const double>(this->values_at_ip_));
}
