/**
 * @file AMRBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class orchestrating adaptive mesh refinement (AMR): defines
 *        the common criteria (error estimator integrator, thresholds,
 *        NC limit) and the shared refine/derefine/update orchestration
 *        (InitialRefine, StepRefine, StepDerefine), delegating the actual
 *        per-strategy error estimation and mesh mutation to derived
 *        classes via the Refine()/Derefine() pure virtual methods.
 * @version 0.1
 * @date 2026-08-04
 *
 * @copyright CEA (C) 2026
 *
 * @anchor AMR
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
#include "AMR/SlothErrorEstimators.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class VAR>
class AMRBase {
 protected:
  mfem::ParMesh& par_mesh_;
  bool is_nc_simplices_;
  double amr_max_elem_error_{0.0};
  int amr_nc_limit_{0};
  int amr_max_preref_cycles_{1};
  int amr_max_level_{0};
  SlothErrorEstimators* amr_estimator_{nullptr};

  virtual bool Refine(VAR& vars) = 0;
  virtual bool Derefine(VAR& vars) = 0;

 public:
  AMRBase(mfem::ParMesh& mesh, bool is_nc_simplices);
  virtual ~AMRBase() = default;

  void SetCriteria(SlothErrorEstimators* estimator, double max_elem_error, int amr_max_level,
                   int nc_limit, int max_preref_cycles);

  void EnsureCriteriaSet() const;

  void InitialRefine(VAR& vars, std::vector<VAR*> auxvars);
  void StepRefine(VAR& vars, std::vector<VAR*> auxvars);
  void StepDerefine(VAR& vars, std::vector<VAR*> auxvars);

  mfem::Array<mfem::Refinement> FilterRefinedCandidates(
      const mfem::Array<mfem::Refinement>& candidates);

 protected:
  void EnsureNCMeshIfNeeded();
};

#include "AMR/AMRBase.tpp"
