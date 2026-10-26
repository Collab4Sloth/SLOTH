/**
 * @file SingleVariableAMR.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief AMR strategy driven by a single, designated variable: error
 *        estimation, refinement and derefinement decisions are all based
 *        exclusively on that one variable's field, using
 *        mfem::ThresholdRefiner/mfem::ThresholdDerefiner directly (which
 *        handle error estimation and mesh mutation together in one call).
 *
 * @warning Choosing this over MultiVariableMaxAMR is only safe if no
 *          other variable sharing the mesh — principal or auxiliary —
 *          needs its own resolution requirements taken into account.
 *          If in doubt, prefer MultiVariableMaxAMR (with a single
 *          variable, it behaves almost identically, minus the preference
 *          for conforming refinement where geometrically possible).
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
#include "AMR/SingleVariableAMR.hpp"
#include "AMR/SlothErrorEstimators.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new SingleVariableAMR<VAR>::SingleVariableAMR object.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param mesh Shared parallel mesh this AMR object will refine/derefine.
 *            The caller retains ownership; this object only keeps a
 *            reference and must not outlive it.
 * @param is_nc_simplices Whether nonconforming refinement is allowed on
 *                        simplex elements of `mesh`. Has no effect on
 *                        quadrilateral/hexahedral meshes.
 * @param var_id Index, within `VAR`, of the single variable whose field
 *              drives all refinement/derefinement decisions.
 */
template <class VAR>
SingleVariableAMR<VAR>::SingleVariableAMR(mfem::ParMesh& mesh, bool is_nc_simplices,
                                          unsigned int var_id)
    : AMRBase<VAR>(mesh, is_nc_simplices), var_id_(var_id) {}

/**
 * @brief Construct a new SingleVariableAMR object whose pilot variable is
 *        designated by its name (the variable associated with the
 *        refinement criterion in the data set).
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param mesh Shared parallel mesh this AMR object will refine/derefine.
 * @param is_nc_simplices Whether nonconforming refinement is allowed on
 *                        simplex elements of `mesh`.
 * @param var_name Name of the pilot variable.
 */
template <class VAR>
SingleVariableAMR<VAR>::SingleVariableAMR(mfem::ParMesh& mesh, bool is_nc_simplices,
                                          const std::string& var_name)
    : AMRBase<VAR>(mesh, is_nc_simplices), var_name_(var_name) {}

/**
 * @brief Refine the mesh once, based solely on the error estimated on the
 *        pilot variable (`var_id_`).
 *
 * @details Builds a fresh estimator via `amr_estimator_->get_value()`,
 *          configures a `mfem::ThresholdRefiner` with a purely local
 *          error goal (`SetTotalErrorFraction(0.0)`,
 *          `SetLocalErrorGoal(amr_max_elem_error_)`) and conforming
 *          refinement preferred where possible
 *          (`PreferConformingRefinement()`), subject to `amr_nc_limit_`.
 *          Candidates are computed via `MarkWithoutRefining()` (without
 *          mutating the mesh), filtered by `FilterRefinedCandidates()`
 *          (max refinement depth), and the decision to mutate is based
 *          on a global reduction (`mesh.ReduceInt()`) rather than the
 *          local candidate count, to stay correct under MPI when some
 *          ranks have no local candidate. The flux spaces of the estimator
 *          are kept by `amr_estimator_` and reused at the next pass
 *          (updated after a mesh change).
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object; only `vars[var_id_]` is consulted.

 * @return true  The mesh was refined.
 * @return false No element exceeded the local error goal; the mesh is
 *               unchanged.
 */
template <class VAR>
bool SingleVariableAMR<VAR>::Refine(VAR& vars) {
  this->EnsureNCMeshIfNeeded();
  this->EnsureCriteriaSet();
  auto& var = *this->PilotVariables(vars).front();
  auto estimator =
      this->amr_estimator_->get_value(var.get_ref_gf(), *var.get_fespace(), this->par_mesh_);

  mfem::ThresholdRefiner refiner(*estimator);
  refiner.SetTotalErrorFraction(0.0);
  refiner.SetLocalErrorGoal(this->amr_max_elem_error_);
  refiner.PreferConformingRefinement();
  refiner.SetNCLimit(this->amr_nc_limit_);

  mfem::Array<mfem::Refinement> refinements;
  refiner.MarkWithoutRefining(this->par_mesh_, refinements);
  refinements = this->FilterRefinedCandidates(refinements);
  int num_marked = this->par_mesh_.ReduceInt(refinements.Size());

  bool did_refine = (num_marked != 0);
  if (did_refine) {
    this->par_mesh_.GeneralRefinement(refinements, -1, this->amr_nc_limit_);
  }

  return did_refine;
}

/**
 * @brief Derefine the mesh once, based on the error estimated on the pilot
 *        variable (`var_id_`) and on the pilot variables of the other AMR
 *        objects coupled to the owning Problem (see
 *        `AMRBase::DerefineByConsensus()`).
 *
 * @details With no coupled pilot variable, this is exactly the
 *          derefinement of `mfem::ThresholdDerefiner` with the threshold
 *          `scale_down_factor * max_elem_error` (sum of the sibling errors).
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object.
 * @param auxvars Auxiliary variables of the owning Problem.
 * @return true  The mesh was derefined.
 * @return false No element met the derefinement threshold; the mesh is
 *               unchanged.
 */
template <class VAR>
bool SingleVariableAMR<VAR>::Derefine(VAR& vars, std::vector<VAR*> auxvars) {
  this->EnsureNCMeshIfNeeded();
  this->EnsureCriteriaSet();
  return this->DerefineByConsensus(vars, auxvars);
}

/**
 * @brief Pilot variable of this AMR object: the variable named `var_name_`
 *        if given, `vars[var_id_]` otherwise.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Variables of the owning Problem.
 */
template <class VAR>
std::vector<typename AMRBase<VAR>::Variable_t*> SingleVariableAMR<VAR>::PilotVariables(
    VAR& vars) {
  if (!this->var_name_.empty()) {
    return {AMRBase<VAR>::FindVariable(vars, this->var_name_)};
  }
  return {&vars[this->var_id_]};
}
