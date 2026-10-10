/**
 * @file MfemFileCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Wraps a FileCoefficient as an mfem::Coefficient evaluated at physical points
 * @version 0.1
 * @date 2026-10-01
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
#include <array>
#include <cmath>
#include <span>
#include <utility>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/FileCoefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Wraps a SLOTH FileCoefficient as an mfem::Coefficient.
 *
 * Unlike MfemCoefficient, which forwards the values of the variables at the
 * integration point, this wrapper forwards the physical coordinates of the
 * integration point to the wrapped coefficient.
 *
 * @note The evaluation time is owned by the wrapped FileCoefficient
 *       (set through its set_time()). The mfem::Coefficient time
 *       (SetTime / GetTime) is not used.
 */
class MfemFileCoefficient : public mfem::Coefficient {
 private:
  ::Coefficient file_coefficient_;
  /// Optional grid functions, only used by ProjectToGridFunction() to get a default FE space
  const std::vector<mfem::ParGridFunction>* gf_ = nullptr;

  mfem::Vector x_;  ///< Physical point, reused across Eval calls (no reallocation)

 public:
  explicit MfemFileCoefficient(::Coefficient coef);

  MfemFileCoefficient(::Coefficient coef, const std::vector<mfem::ParGridFunction>& vu);
  /// Forbid temporaries: only a pointer to vu is stored
  MfemFileCoefficient(::Coefficient coef, std::vector<mfem::ParGridFunction>&& vu) = delete;

  double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;

  mfem::ParGridFunction ProjectToGridFunction(mfem::ParFiniteElementSpace& fespace);
  mfem::ParGridFunction ProjectToGridFunction();

  ~MfemFileCoefficient() override = default;
};