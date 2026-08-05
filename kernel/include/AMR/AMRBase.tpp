/**
 * @anchor AMR
 * @file AMRBase.tpp
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

/**
 * @brief Construct a new AMRBase<VAR>::AMRBase object.
 *
 * @details Does not activate nonconforming mode on `mesh` by itself
 *          (see `EnsureNCMeshIfNeeded()`, called lazily on the first
 *          refine/derefine attempt).
 *
 * @tparam VAR Variables container type this AMR strategy operates on.
 * @param mesh Shared parallel mesh this AMR object will refine/derefine.
 *            The caller retains ownership; this object only keeps a
 *            reference and must not outlive it.
 * @param is_nc_simplices Whether nonconforming refinement is allowed on
 *                        simplex elements of `mesh` (see
 *                        `SpatialDiscretization::is_nc_simplices()`). Has
 *                        no effect on quadrilateral/hexahedral meshes.
 */
template <class VAR>
AMRBase<VAR>::AMRBase(mfem::ParMesh& mesh, bool is_nc_simplices)
    : par_mesh_(mesh), is_nc_simplices_(is_nc_simplices) {}

/**
 * @brief Set the error-estimation criteria and refinement parameters used
 *        by this AMR strategy.
 *
 * @details Must be called before `InitialRefine()`, `StepRefine()`, or
 *          `StepDerefine()` are used; `EnsureCriteriaSet()` verifies this
 *          precondition and aborts with a clear message if it was not
 *          respected. The caller retains ownership of `integ` and must
 *          keep it alive for as long as this AMR object is used.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param integ Integrator used by the underlying error estimator
 *             (`mfem::L2ZienkiewiczZhuEstimator`) to compute per-element
 *             flux/error, e.g. a `mfem::DiffusionIntegrator`.
 * @param max_elem_error Local error threshold above which an element is
 *                       marked for refinement (used directly for
 *                       refinement; scaled down internally for
 *                       derefinement).
 * @param nc_limit Maximum allowed refinement level difference between
 *                neighboring elements (2:1 balance constraint); `0` means
 *                unlimited.
 * @param max_preref_cycles Maximum number of refinement cycles applied to
 *                          the initial condition by `InitialRefine()`,
 *                          before the time-stepping loop starts.
 */
template <class VAR>
void AMRBase<VAR>::SetCriteria(mfem::BilinearFormIntegrator* integ, double max_elem_error,
                               int nc_limit, int max_preref_cycles) {
  this->amr_integ_ = integ;
  this->amr_max_elem_error_ = max_elem_error;
  this->amr_nc_limit_ = nc_limit;
  this->amr_max_preref_cycles_ = max_preref_cycles;
}

/**
 * @brief
 *
 * @tparam VAR
 */
template <class VAR>
void AMRBase<VAR>::EnsureCriteriaSet() const {
  MFEM_VERIFY(this->amr_integ_ != nullptr,
              "AMR criteria not set: call SetCriteria() on this AMR object before "
              "using Refine()/Derefine().");
}

/**
 * @brief Verify that SetCriteria() has already been called on this AMR
 *        object.
 *
 * @details Called internally at the start of every concrete Refine()/
 *          Derefine() implementation, before `amr_integ_` (or any other
 *          criteria member) is dereferenced. Aborts with a clear error
 *          message if the precondition is not met, rather than letting a
 *          null `amr_integ_` be dereferenced silently.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 */
template <class VAR>
void AMRBase<VAR>::EnsureNCMeshIfNeeded() {
  if (!this->par_mesh_.Nonconforming()) {
    this->par_mesh_.EnsureNCMesh(is_nc_simplices_);
  }
}

/**
 * @brief Activate nonconforming (NC) mode on the mesh if it is not already
 *        active.
 *
 * @details Called internally at the start of every concrete Refine()/
 *          Derefine() implementation. In practice this should be a no-op
 *          in most workflows: `EnsureNCMesh()` is expected to have already
 *          been called on the serial mesh, before it was distributed into
 *          a `mfem::ParMesh` (see `SpatialDiscretization`), since a
 *          conforming `ParMesh` cannot be converted to nonconforming
 *          after construction in parallel. This call remains here as a
 *          defensive fallback for the sequential case, where the
 *          conversion is still possible after the fact.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Primary variables sharing the mesh managed by this AMR
 *            object.
 * @param auxvars Auxiliary variable collections sharing the same mesh,
 *               if any (empty by default).
 */
template <class VAR>
void AMRBase<VAR>::InitialRefine(VAR& vars, std::vector<VAR*> auxvars) {
  int cycle = 0;
  bool did_refine = false;
  do {
    did_refine = this->Refine(vars);

    if (did_refine) {
      for (std::size_t i = 0; i < vars.get_variables_number(); i++) {
        auto& var = vars[i];
        var.UpdateAndRebalance();
        var.setInitialCondition();
      }
      for (auto* auxvar_container : auxvars) {
        for (std::size_t i = 0; i < auxvar_container->get_variables_number(); i++) {
          (*auxvar_container)[i].UpdateAndRebalance();
        }
      }
    }
    cycle++;
  } while (did_refine && cycle < this->amr_max_preref_cycles_);
  SlothInfo::verbose("Pre-refinement: ", cycle, " cycles applied on initial condition.");
}

/**
 * @brief Apply one refinement pass during the time-stepping loop.
 *
 * @details Calls the derived class's `Refine()` implementation, and if it
 *          reports that the mesh was actually refined, resynchronizes
 *          every variable in `vars` via `Variable::UpdateAndRebalance()` —
 *          not just the variable(s) that may have driven the refinement
 *          decision, since all variables share the same mesh and must
 *          stay consistent with its current state.
 *
 * @note Intended to be called once per time step, typically paired with a
 *       preceding `StepDerefine()` call (see `Problem::save()`), never
 *       interleaved with another mesh-mutating call without an
 *       `UpdateAndRebalance()` pass in between.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object.
 * @param auxvars Auxiliary variable collections sharing the same mesh,
 *               if any (empty by default).
 */
template <class VAR>
void AMRBase<VAR>::StepRefine(VAR& vars, std::vector<VAR*> auxvars) {
  bool did_refine = this->Refine(vars);

  if (did_refine) {
    for (std::size_t i = 0; i < vars.get_variables_number(); i++) {
      auto& var = vars[i];
      var.UpdateAndRebalance();
    }
    for (auto* auxvar_container : auxvars) {
      for (std::size_t i = 0; i < auxvar_container->get_variables_number(); i++) {
        (*auxvar_container)[i].UpdateAndRebalance();
      }
    }
  }
}

/**
 * @brief Apply one derefinement pass during the time-stepping loop.
 *
 * @details Calls the derived class's `Derefine()` implementation, and if
 *          it reports that the mesh was actually derefined, resynchronizes
 *          every variable in `vars` via `Variable::UpdateAndRebalance()` —
 *          not just the variable(s) that may have driven the derefinement
 *          decision, since all variables share the same mesh and must
 *          stay consistent with its current state.
 *
 * @note Intended to be called once per time step, typically followed by a
 *       `StepRefine()` call (see `Problem::save()`): derefine first, to
 *       free up room for refinement in areas that newly need it.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object.
 * @param auxvars Auxiliary variable collections sharing the same mesh,
 *               if any (empty by default).
 */
template <class VAR>
void AMRBase<VAR>::StepDerefine(VAR& vars, std::vector<VAR*> auxvars) {
  bool did_derefine = this->Derefine(vars);

  if (did_derefine) {
    for (std::size_t i = 0; i < vars.get_variables_number(); i++) {
      auto& var = vars[i];
      var.UpdateAndRebalance();
    }
    for (auto* auxvar_container : auxvars) {
      for (std::size_t i = 0; i < auxvar_container->get_variables_number(); i++) {
        (*auxvar_container)[i].UpdateAndRebalance();
      }
    }
  }
}
