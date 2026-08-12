/**
 * @file SingleVariableAMR.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief AMR strategy driven by a single, designated variable: error
 *        estimation, refinement and derefinement decisions are all based
 *        exclusively on that one variable's field, using
 *        mfem::ThresholdRefiner/mfem::ThresholdDerefiner directly (which
 *        handle error estimation and mesh mutation together in one call).
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
class SingleVariableAMR : public AMRBase<VAR> {
 private:
  unsigned int var_id_;

 public:
  SingleVariableAMR(mfem::ParMesh& mesh, bool is_nc_simplices, unsigned int var_id);
  virtual ~SingleVariableAMR() = default;

  bool Refine(VAR& vars) final;
  bool Derefine(VAR& vars) final;
};

#include "AMR/SingleVariableAMR.tpp"
