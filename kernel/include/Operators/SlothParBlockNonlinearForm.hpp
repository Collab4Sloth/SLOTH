/**
 * @file SlothParBlockNonlinearForm.hpp
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

#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief mfem::ParBlockNonlinearForm whose integrators can prepare the values at the quadrature
 * points of the mesh before each assembly.
 *
 * Before the assembly of the residual (Mult) and of the Jacobian (GetGradient), the integrators
 * whose uses_quadrature_points() is true receive the current values of the unknowns at the local
 * dofs of each block (SlothNLFormIntegrator::prepare_quadrature_points). The assembly itself is
 * the one of mfem::ParBlockNonlinearForm.
 *
 * GetLocalGradient (not virtual in MFEM, not used in SLOTH) does not prepare the integrators.
 *
 * Without integrator using the quadrature points, nothing is done before the assembly.
 *
 * @tparam VARS Template parameter defining the variables
 */
template <class VARS>
class SlothParBlockNonlinearForm : public mfem::ParBlockNonlinearForm {
 public:
  explicit SlothParBlockNonlinearForm(mfem::Array<mfem::ParFiniteElementSpace*>& fes);

  void AddDomainIntegrator(SlothNLFormIntegrator<VARS>* nlfi);
  void AddBoundaryIntegrator(SlothNLFormIntegrator<VARS>* nlfi, mfem::Array<int>& bdr_marker);

  void Mult(const mfem::Vector& x, mfem::Vector& y) const override;
  mfem::BlockOperator& GetGradient(const mfem::Vector& x) const override;

 private:
  void prepare_quadrature_points(const mfem::Vector& x) const;

  std::vector<SlothNLFormIntegrator<VARS>*> integrators_;  ///< domain and boundary (not owned)
  mutable mfem::BlockVector x_true_;   ///< true dofs of each block (view on x)
  mutable mfem::BlockVector x_local_;  ///< local dofs of each block
};

#include "Operators/SlothParBlockNonlinearForm.tpp"
