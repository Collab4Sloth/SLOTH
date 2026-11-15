/**
 * @file SlothParBlockNonlinearForm.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Block nonlinear form calling the integrators before each assembly
 * @version 0.1
 * @date 2026-11-14
 *
 * Copyright CEA (C) 2026
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

#include <algorithm>
#include <vector>

#include "Operators/SlothParBlockNonlinearForm.hpp"

template <class VARS>
SlothParBlockNonlinearForm<VARS>::SlothParBlockNonlinearForm(
    mfem::Array<mfem::ParFiniteElementSpace*>& fes)
    : mfem::ParBlockNonlinearForm(fes) {}

/**
 * @brief Adds a domain integrator (owned by the form, as in mfem::ParBlockNonlinearForm)
 */
template <class VARS>
void SlothParBlockNonlinearForm<VARS>::AddDomainIntegrator(SlothNLFormIntegrator<VARS>* nlfi) {
  this->integrators_.push_back({nlfi, nullptr});
  mfem::ParBlockNonlinearForm::AddDomainIntegrator(nlfi);
}

/**
 * @brief Adds a boundary integrator restricted to bdr_marker (owned by the form, as in
 * mfem::ParBlockNonlinearForm)
 */
template <class VARS>
void SlothParBlockNonlinearForm<VARS>::AddBoundaryIntegrator(SlothNLFormIntegrator<VARS>* nlfi,
                                                             mfem::Array<int>& bdr_marker) {
  this->integrators_.push_back({nlfi, &bdr_marker});
  mfem::ParBlockNonlinearForm::AddBoundaryIntegrator(nlfi, bdr_marker);
}

/**
 * @brief Residual: integrators prepared at x, then assembly of mfem::ParBlockNonlinearForm
 */
template <class VARS>
void SlothParBlockNonlinearForm<VARS>::Mult(const mfem::Vector& x, mfem::Vector& y) const {
  Catch_Time_Section("SlothParBlockNonlinearForm::Mult");
  const bool prepared = this->prepare_quadrature_points(x, Assembly::Residual);
  mfem::ParBlockNonlinearForm::Mult(x, y);
  if (prepared) this->clear_quadrature_points();
}

/**
 * @brief Jacobian: integrators prepared at x, then assembly of mfem::ParBlockNonlinearForm
 */
template <class VARS>
mfem::BlockOperator& SlothParBlockNonlinearForm<VARS>::GetGradient(const mfem::Vector& x) const {
  Catch_Time_Section("SlothParBlockNonlinearForm::GetGradient");
  const bool prepared = this->prepare_quadrature_points(x, Assembly::Jacobian);
  mfem::BlockOperator& gradient = mfem::ParBlockNonlinearForm::GetGradient(x);
  if (prepared) this->clear_quadrature_points();
  return gradient;
}

/**
 * @brief Gives the local dofs of each block (same values as the element vectors of the
 * assembly) to the integrators using the quadrature points of the mesh.
 *
 * @param x True dofs of each block
 * @param assembly Residual or Jacobian
 * @return false if no integrator uses the quadrature points (nothing is done)
 */
template <class VARS>
bool SlothParBlockNonlinearForm<VARS>::prepare_quadrature_points(const mfem::Vector& x,
                                                                 const Assembly assembly) const {
  const auto uses_quadrature_points = [](const Entry& entry) {
    return entry.nlfi->uses_quadrature_points();
  };
  if (std::ranges::none_of(this->integrators_, uses_quadrature_points)) {
    return false;
  }

  // As in mfem::ParBlockNonlinearForm::Mult: x is not modified
  this->x_true_.Update(const_cast<mfem::Vector&>(x), this->block_trueOffsets);
  this->x_local_.Update(this->block_offsets);
  for (int s = 0; s < this->fes.Size(); ++s) {
    this->fes[s]->GetProlongationMatrix()->Mult(this->x_true_.GetBlock(s),
                                                this->x_local_.GetBlock(s));
  }
  for (const Entry& entry : this->integrators_) {
    if (uses_quadrature_points(entry)) {
      entry.nlfi->prepare_quadrature_points(this->fes, this->x_local_, assembly,
                                            entry.bdr_marker);
    }
  }
  return true;
}

/**
 * @brief End of the assembly: the prepared values are no longer used by the integrators
 */
template <class VARS>
void SlothParBlockNonlinearForm<VARS>::clear_quadrature_points() const {
  for (const Entry& entry : this->integrators_) {
    if (entry.nlfi->uses_quadrature_points()) {
      entry.nlfi->clear_quadrature_points();
    }
  }
}
