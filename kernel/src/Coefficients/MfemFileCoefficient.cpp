/**
 * @file MfemFileCoefficient.cpp
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
#include "Coefficients/MfemFileCoefficient.hpp"

#include <array>
#include <cmath>
#include <span>
#include <utility>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/FileCoefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct from an existing file-based Coefficient, without grid functions.
 *
 * @param coef Coefficient to wrap (copied; the tabulated field is shared). Must be file-based.
 */
MfemFileCoefficient::MfemFileCoefficient(::Coefficient coef) : file_coefficient_(std::move(coef)) {
  MFEM_VERIFY(this->file_coefficient_.is_file(),
              "MfemFileCoefficient: the wrapped coefficient is not file-based (Scheme::File).");
}
/**
 * @brief Construct from an existing file-based Coefficient.
 *
 * @param coef Coefficient to wrap (copied; the tabulated field is shared). Must be file-based.
 * @param vu   Grid functions used only to get a default FE space in ProjectToGridFunction().
 *             Only a pointer is stored: vu must outlive this object.
 */
MfemFileCoefficient::MfemFileCoefficient(::Coefficient coef,
                                         const std::vector<mfem::ParGridFunction>& vu)
    : MfemFileCoefficient(std::move(coef)) {
  this->gf_ = &vu;
}
/**
 * @brief Evaluate the tabulated field at an integration point.
 *
 * Maps the reference integration point to physical coordinates with
 * T.Transform(), then interpolates the field at the time held by the wrapped FileCoefficient.
 *
 * @param T  Element transformation.
 * @param ip Integration point (reference coordinates).
 * @return double Field value at the physical point.
 */
double MfemFileCoefficient::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) {
  T.SetIntPoint(&ip);
  T.Transform(ip, this->x_);

  return this->file_coefficient_.compute(
      std::span<const double>(this->x_.GetData(), this->x_.Size()), std::span<const double>());
}

/**
 * @brief Project the coefficient onto a grid function.
 *
 * @param fespace Finite element space to project onto.
 * @return mfem::ParGridFunction Grid function holding the projected values.
 */
mfem::ParGridFunction MfemFileCoefficient::ProjectToGridFunction(
    mfem::ParFiniteElementSpace& fespace) {
  mfem::ParGridFunction result(&fespace);
  result.ProjectCoefficient(*this);
  return result;
}

/**
 * @brief Project the coefficient onto a grid function, using the finite
 *        element space of the first grid function given at construction
 *        (assumption of homogeneous FESpace).
 *
 * @return Grid function holding the projected values.
 */
mfem::ParGridFunction MfemFileCoefficient::ProjectToGridFunction() {
  MFEM_VERIFY(this->gf_ != nullptr && !this->gf_->empty(),
              "MfemFileCoefficient::ProjectToGridFunction(): no grid function given at "
              "construction. Use ProjectToGridFunction(fespace) instead.");
  return ProjectToGridFunction(*(*this->gf_)[0].ParFESpace());
}