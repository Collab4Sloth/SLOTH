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
