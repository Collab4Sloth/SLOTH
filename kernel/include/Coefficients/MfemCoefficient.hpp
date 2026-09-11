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

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/Coefficients.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Wraps a SLOTH ::Coefficient as an mfem::Coefficient, so it can be
 *        projected onto a grid function or plugged into an mfem
 *        integrator.
 */
class MfemCoefficient : public mfem::Coefficient {
 private:
  int id_ = 0;
  ::Coefficient sloth_coefficient_;

  const std::vector<mfem::ParGridFunction>& gf_;
  const std::vector<mfem::ParGridFunction>& gfn_;
  const std::vector<mfem::ParGridFunction>& vaux_gf_;

  mutable std::vector<double> vaux_gf_at_ip_;
  mutable std::vector<double> gf_at_ip_;
  mutable std::vector<double> gfn_at_ip_;

 public:
  MfemCoefficient(::Coefficient coef, const std::vector<mfem::ParGridFunction>& vu,
                  const std::vector<mfem::ParGridFunction>& vun,
                  const std::vector<mfem::ParGridFunction>& vauxn);
  double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;

  mfem::ParGridFunction ProjectToGridFunction(mfem::ParFiniteElementSpace& fespace);
  mfem::ParGridFunction ProjectToGridFunction();

  virtual ~MfemCoefficient() = default;
};

/**
 * @brief Construct a new MfemCoefficient object.
 *
 * @param coef SLOTH coefficient to wrap.
 * @param vu Current-time vector of grid function for the coefficient's own variable.
 * @param vun Previous-time vector of grid function for the coefficient's own variable.
 * @param vauxn Grid functions of the auxiliary variables.
 */
MfemCoefficient::MfemCoefficient(::Coefficient coef, const std::vector<mfem::ParGridFunction>& vu,
                                 const std::vector<mfem::ParGridFunction>& vun,
                                 const std::vector<mfem::ParGridFunction>& vauxn)
    : sloth_coefficient_(coef), gf_(vu), gfn_(vun), vaux_gf_(vauxn) {
  MFEM_VERIFY(this->gf_.size() == this->gfn_.size(),
              "Error of consistency between Current-time and Previous-time vectors of grid "
              "function. Please check their sizes.");
}

/**
 * @brief Evaluate the wrapped coefficient at an integration point.
 *
 * Dispatches on the coefficient's scheme (scalar/semi-implicit/
 * explicit/implicit) to call the matching ::Coefficient::compute overload.
 *
 * @param T Element transformation at the integration point.
 * @param ip Integration point.
 * @return Coefficient value.
 */
double MfemCoefficient::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) {
  if (this->sloth_coefficient_.is_scalar()) {
    return this->sloth_coefficient_.compute();
  }

  this->vaux_gf_at_ip_.resize(this->vaux_gf_.size());
  for (size_t k = 0; k < this->vaux_gf_.size(); ++k) {
    this->vaux_gf_at_ip_[k] = this->vaux_gf_[k].GetValue(T, ip);
  }

  this->gf_at_ip_.resize(this->gf_.size());
  this->gfn_at_ip_.resize(this->gfn_.size());
  for (size_t k = 0; k < this->gf_.size(); ++k) {
    this->gf_at_ip_[k] = this->gf_[k].GetValue(T, ip);
    this->gfn_at_ip_[k] = this->gfn_[k].GetValue(T, ip);
  }

  if (this->sloth_coefficient_.is_semi_implicit()) {
    return this->sloth_coefficient_.compute(std::span<const double>(this->gf_at_ip_),
                                            std::span<const double>(this->gfn_at_ip_),
                                            std::span<const double>(this->vaux_gf_at_ip_));
  } else if (this->sloth_coefficient_.is_explicit()) {
    return this->sloth_coefficient_.compute(std::span<const double>(this->gfn_at_ip_),
                                            std::span<const double>(this->vaux_gf_at_ip_));
  } else {
    return this->sloth_coefficient_.compute(std::span<const double>(this->gf_at_ip_),
                                            std::span<const double>(this->vaux_gf_at_ip_));
  }
}
/**
 * @brief Project the coefficient onto a grid function on a given space.
 *
 * @param fespace Finite element space to project onto.
 * @return Grid function holding the projected values.
 */
mfem::ParGridFunction MfemCoefficient::ProjectToGridFunction(mfem::ParFiniteElementSpace& fespace) {
  mfem::ParGridFunction result(&fespace);
  result.ProjectCoefficient(*this);
  return result;
}

/**
 * @brief Project the coefficient onto a grid function, using the finite
 *        element space of the first index in `gf_` (assumption of homogeneous FESpace)
 *
 * @return Grid function holding the projected values.
 */
mfem::ParGridFunction MfemCoefficient::ProjectToGridFunction() {
  return ProjectToGridFunction(*gf_[0].ParFESpace());
}
