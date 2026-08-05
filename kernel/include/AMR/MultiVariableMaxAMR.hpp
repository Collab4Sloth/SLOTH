/**
 * @file MultiVariableMaxAMR.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief AMR strategy combining several variables into a single, shared
 *        mesh mutation: for each element, the local error is taken as the
 *        maximum across all variables in the collection, and refinement/
 *        derefinement decisions are then made once from this combined
 *        criterion (rather than mutating the mesh once per variable).
 *        Suitable for coupled multiphase-field problems where every
 *        variable must stay resolved on the mesh, even if only some of
 *        them are responsible for a given refinement decision.
 * @version 0.1
 * @date 2026-08-04
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
#include "AMR/AMRBase.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class VAR>
class MultiVariableMaxAMR : public AMRBase<VAR> {
 public:
  MultiVariableMaxAMR(mfem::ParMesh& mesh, bool is_nc_simplices);
  virtual ~MultiVariableMaxAMR() = default;

  bool Refine(VAR& vars) override;
  bool Derefine(VAR& vars) override;
};

#include "AMR/MultiVariableMaxAMR.tpp"
