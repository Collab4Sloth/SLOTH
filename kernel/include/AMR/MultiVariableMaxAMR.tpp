/**
 * @file MultiVariableMaxAMR.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  AMR strategy combining several variables into a single, shared
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
#include "AMR/MultiVariableMaxAMR.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new MultiVariableMaxAMR<VAR>::MultiVariableMaxAMR
 *        object.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param mesh Shared parallel mesh this AMR object will refine/derefine.
 *            The caller retains ownership; this object only keeps a
 *            reference and must not outlive it.
 * @param is_nc_simplices Whether nonconforming refinement is allowed on
 *                        simplex elements of `mesh`. Has no effect on
 *                        quadrilateral/hexahedral meshes.
 */
template <class VAR>
MultiVariableMaxAMR<VAR>::MultiVariableMaxAMR(mfem::ParMesh& mesh, bool is_nc_simplices)
    : AMRBase<VAR>(mesh, is_nc_simplices) {}

/**
 * @brief Refine the shared mesh once, based on the per-element maximum
 *        error across all variables in `vars`.
 *
 * @details For each variable, builds a fresh estimator via
 *          `amr_estimator_->get_value()` and reads its per-element local
 *          errors (this only estimates the error; it does not mutate the
 *          mesh). These are combined element-wise into `combined_error`
 *          via `max()`, so that an element is a refinement candidate as
 *          soon as ANY variable needs it resolved there. Candidates are
 *          filtered by `FilterRefinedCandidates()` (max refinement
 *          depth), and the decision to mutate is based on a global
 *          reduction (`mesh.ReduceInt()`) rather than the local
 *          candidate count, to stay correct under MPI when some ranks
 *          have no local candidate. The mesh is then mutated exactly
 *          once, via `GeneralRefinement()` — unlike calling
 *          `mfem::ThresholdRefiner::Apply()` separately per variable,
 *          which would mutate the mesh once per variable and require
 *          resynchronizing every other variable's finite element space
 *          in between.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object.
 * @return true  At least one element was marked and the mesh was refined.
 * @return false No element exceeded the error threshold; the mesh is
 *               unchanged.
 */
template <class VAR>
bool MultiVariableMaxAMR<VAR>::Refine(VAR& vars) {
  this->EnsureNCMeshIfNeeded();
  this->EnsureCriteriaSet();
  const int ne = this->par_mesh_.GetNE();
  mfem::Vector combined_error(ne);
  combined_error = 0.0;

  auto accumulate_error = [&](auto& var) {
    auto estimator =
        this->amr_estimator_->get_value(var.get_ref_gf(), *var.get_fespace(), this->par_mesh_);
    const mfem::Vector& local_errors = estimator->GetLocalErrors();
    for (int e = 0; e < ne; e++) {
      combined_error(e) = std::max(combined_error(e), local_errors(e));
    }
  };

  for (unsigned int iv = 0; iv < vars.get_variables_number(); iv++) {
    accumulate_error(vars[iv]);
  }

  mfem::Array<mfem::Refinement> refinements;
  for (int e = 0; e < ne; e++) {
    if (combined_error(e) > this->amr_max_elem_error_) {
      refinements.Append(mfem::Refinement(e));
    }
  }

  refinements = this->FilterRefinedCandidates(refinements);
  int num_marked = this->par_mesh_.ReduceInt(refinements.Size());
  bool did_refine = (num_marked != 0);

  if (did_refine) {
    this->par_mesh_.GeneralRefinement(refinements, -1, this->amr_nc_limit_);
  }

  return did_refine;
}

/**
 * @brief Derefine the shared mesh once, based on the per-element maximum
 *        error across all variables in `vars`.
 *
 * @details Mirrors `Refine()`: for each variable, builds a fresh
 *          estimator via `amr_estimator_->get_value()`, and the
 *          per-element error is combined via `max()` across all
 *          variables, so that an element is a derefinement candidate
 *          only if its error stays low for EVERY variable (a
 *          conservative criterion, since derefining would otherwise risk
 *          losing resolution needed by a variable that did not
 *          contribute to the decision). The combined error vector is
 *          then passed once to `mfem::ParMesh::DerefineByError()`, which
 *          internally aggregates by SUM the errors of sibling elements
 *          sharing the same parent before comparing to
 *          `derefine_threshold` — hence the threshold is scaled down
 *          (`this->amr_scale_down_factor_ *`) relative to `amr_max_elem_error_` used for
 *          refinement, and may need further empirical recalibration
 *          (e.g. multiplied by the number of siblings) depending on the
 *          observed derefinement behavior.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object.
 * @param auxvars Auxiliary variable collections whose error also
 *               contributes to `combined_error`, if any (empty by
 *               default). Only taken into account by the caller when
 *               `AMRBase::SetIncludeAuxVars(true)` was set.
 * @return true  At least one group of elements was derefined.
 * @return false No group of elements met the derefinement criterion; the
 *               mesh is unchanged.
 */
template <class VAR>
bool MultiVariableMaxAMR<VAR>::Derefine(VAR& vars, std::vector<VAR*> auxvars) {
  this->EnsureNCMeshIfNeeded();
  this->EnsureCriteriaSet();

  const int ne = this->par_mesh_.GetNE();
  mfem::Vector combined_error(ne);
  combined_error = 0.0;

  auto accumulate_error = [&](auto& var) {
    auto estimator =
        this->amr_estimator_->get_value(var.get_ref_gf(), *var.get_fespace(), this->par_mesh_);
    const mfem::Vector& local_errors = estimator->GetLocalErrors();
    for (int e = 0; e < ne; e++) {
      combined_error(e) = std::max(combined_error(e), local_errors(e));
    }
  };
  for (unsigned int iv = 0; iv < vars.get_variables_number(); iv++) {
    accumulate_error(vars[iv]);
  }
  for (auto* auxvar : auxvars) {
    for (std::size_t i = 0; i < auxvar->get_variables_number(); i++) {
      accumulate_error((*auxvar)[i]);
    }
  }

  const double derefine_threshold = this->amr_scale_down_factor_ * this->amr_max_elem_error_;

  bool did_derefine =
      this->par_mesh_.DerefineByError(combined_error, derefine_threshold, this->amr_nc_limit_);

  return did_derefine;
}
