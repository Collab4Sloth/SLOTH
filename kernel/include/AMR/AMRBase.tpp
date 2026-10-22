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
#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "AMR/AMRBase.hpp"
#include "AMR/SlothErrorEstimators.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AMRBase<VAR>::AMRBase object.
 *
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
 * @brief Destroy the AMRBase object, removing its pilot variables from the
 *        registry.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 */
template <class VAR>
AMRBase<VAR>::~AMRBase() {
  auto& registry = PilotRegistry();
  auto it_mesh = registry.find(&this->par_mesh_);
  if (it_mesh == registry.end()) {
    return;
  }
  auto& owners = it_mesh->second;
  for (auto it = owners.begin(); it != owners.end();) {
    it = (it->second == this) ? owners.erase(it) : std::next(it);
  }
  if (owners.empty()) {
    registry.erase(it_mesh);
  }
}

/**
 * @brief Registry of the pilot variables of every AMR object, per mesh.
 *
 * @details A variable takes part in the AMR as soon as its Problem
 *          associates it with a refinement criterion (pilot variable of
 *          the AMR object of this Problem). The registry gives, for each
 *          mesh, the AMR object driven by each pilot variable, so that the
 *          derefinement decision of an AMR object can require the
 *          agreement of the pilot variables of the other Problems.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 */
template <class VAR>
typename AMRBase<VAR>::PilotRegistry_t& AMRBase<VAR>::PilotRegistry() {
  static PilotRegistry_t registry;
  return registry;
}

/**
 * @brief Find a variable of `vars` by its name.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Variables of the owning Problem.
 * @param name Name of the variable (as given in the data set).
 * @return Pointer to the variable (error if not found).
 */
template <class VAR>
typename AMRBase<VAR>::Variable_t* AMRBase<VAR>::FindVariable(VAR& vars, const std::string& name) {
  for (std::size_t i = 0; i < vars.get_variables_number(); i++) {
    if (vars[i].getVariableName() == name) {
      return &vars[i];
    }
  }
  mfem::mfem_error(("AMR: pilot variable '" + name + "' not found in the variables of the problem.")
                       .c_str());
  return nullptr;
}

/**
 * @brief Register the pilot variables of this AMR object.
 *
 * @details Called at each AMR stage (cheap); all the AMR objects of a
 *          coupling are registered by `InitialRefine()`, before any
 *          derefinement.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Variables of the Problem owning this AMR object.
 */
template <class VAR>
void AMRBase<VAR>::RegisterPilots(VAR& vars) {
  auto& owners = PilotRegistry()[&this->par_mesh_];
  for (auto* pilot : this->PilotVariables(vars)) {
    owners[static_cast<const void*>(pilot)] = this;
  }
}

/**
 * @brief Accumulate (max) into `combined_error` the local errors of `var`,
 *        estimated with the estimator of this AMR object and normalized by
 *        its derefinement threshold.
 *
 * @details Each pilot variable is thus judged with the criterion of the AMR
 *          object it drives (estimator, `max_elem_error`,
 *          `scale_down_factor`): a normalized error lower than 1 means
 *          "derefinement allowed" for this variable.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param var Pilot variable of this AMR object.
 * @param combined_error Per-element normalized error, updated in place.
 */
template <class VAR>
void AMRBase<VAR>::AccumulateNormalizedError(Variable_t& var, mfem::Vector& combined_error) {
  this->EnsureCriteriaSet();
  const double threshold = this->amr_scale_down_factor_ * this->amr_max_elem_error_;
  MFEM_VERIFY(threshold > 0., "AMRBase: derefinement threshold must be strictly positive.");
  auto estimator =
      this->amr_estimator_->get_value(var.get_ref_gf(), *var.get_fespace(), this->par_mesh_);
  const mfem::Vector& local_errors = estimator->GetLocalErrors();
  const int ne = combined_error.Size();
  MFEM_VERIFY(local_errors.Size() == ne, "AMRBase: error vector size mismatch.");
  const double inv_threshold = 1. / threshold;
  for (int e = 0; e < ne; e++) {
    combined_error(e) = std::max(combined_error(e), local_errors(e) * inv_threshold);
  }
}

/**
 * @brief Derefine the mesh once, only where every pilot variable concerned
 *        agrees.
 *
 * @details The variables concerned are:
 *          - the pilot variables of this AMR object;
 *          - the auxiliary variables of the owning Problem (`coupled_vars`)
 *            that are pilot variables of an AMR object on the same mesh
 *            (i.e. associated with a refinement criterion in their own
 *            Problem). Each of them is judged with the criterion of that
 *            AMR object.
 *          The other variables (auxiliary variables without refinement
 *          criterion, variables of non coupled Problems, MPI rank field...)
 *          are ignored. The per-element normalized errors are combined via
 *          `max()`, then `mfem::ParMesh::DerefineByError()` is applied once
 *          with a unit threshold (sum of the sibling errors, as
 *          `mfem::ThresholdDerefiner`): for a single variable, this is
 *          exactly the derefinement with the threshold
 *          `scale_down_factor * max_elem_error`.
 *
 *          The variables are visited in the order of `vars` then
 *          `coupled_vars`, identical on every rank (the estimators may
 *          exchange data between ranks).
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Variables of the Problem owning this AMR object.
 * @param coupled_vars Auxiliary variables of the owning Problem.
 * @return true if the mesh was derefined.
 */
template <class VAR>
bool AMRBase<VAR>::DerefineByConsensus(VAR& vars, const std::vector<VAR*>& coupled_vars) {
  this->RegisterPilots(vars);
  const int ne = this->par_mesh_.GetNE();
  mfem::Vector combined_error(ne);
  combined_error = 0.0;

  std::vector<const void*> done;
  auto add = [&done, &combined_error](Variable_t& var, AMRBase<VAR>* owner) {
    const void* key = static_cast<const void*>(&var);
    if (std::find(done.begin(), done.end(), key) != done.end()) {
      return;
    }
    done.push_back(key);
    owner->AccumulateNormalizedError(var, combined_error);
  };

  for (auto* pilot : this->PilotVariables(vars)) {
    add(*pilot, this);
  }
  const auto& owners = PilotRegistry()[&this->par_mesh_];
  for (auto* container : coupled_vars) {
    if (container == nullptr) {
      continue;
    }
    for (std::size_t i = 0; i < container->get_variables_number(); i++) {
      auto& var = (*container)[i];
      if (!this->IsOnMesh(var.get_fespace())) {
        continue;
      }
      auto it = owners.find(static_cast<const void*>(&var));
      if (it != owners.end()) {
        add(var, it->second);
      }
    }
  }

  return this->par_mesh_.DerefineByError(combined_error, 1.0, this->amr_nc_limit_);
}

/**
 * @brief Set the error-estimation criteria and refinement parameters used
 *        by this AMR strategy.
 *
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param amr_params Collection of Parameter objects for AMR :
 *             - max_elem_error (mandatory): Local error threshold above which an element is
 *                       marked for refinement (used directly for
 *                       refinement; scaled down internally for
 *                       derefinement).
 *             - amr_max_level (optional): Maximum refinement depth allowed for any single
 *                      element. `0` (or negative) means no limit — an
 *                      element could then be refined at every time step
 *                      with no upper bound.
 *             - nc_limit (optional): Maximum allowed refinement level difference between
 *                neighboring elements (2:1 balance constraint); `0` means
 *                unlimited.
 *             - max_preref_cycles (optional):  Maximum number of refinement cycles applied to
 *                          the initial condition by `InitialRefine()`,
 *                          before the time-stepping loop starts. Must be
 *                          strictly positive to enable refinement cycles.
 *             - scale_down_factor (optional): threshold scaled down relative to
 *                          `amr_max_elem_error_` used for refinement, to provide
 *                          hysteresis and avoid oscillating refine/derefine cycles
 *                          on the same elements
 *             - rebalance (optional, default false): rebalance the parallel mesh
 *                          (`mfem::ParMesh::Rebalance()`, nonconforming meshes only)
 *                          after each AMR pass that modified the mesh
 */
template <class VAR>
void AMRBase<VAR>::SetCriteria(SlothErrorEstimators* estimator, const Parameters& params) {
  // Mandatory
  this->amr_estimator_ = estimator;

  this->amr_max_elem_error_ = params.template get_param_value<double>("max_elem_error");

  // Optional
  this->amr_nc_limit_ = params.template get_param_value_or_default<int>("nc_limit", 0);
  this->amr_max_level_ = params.template get_param_value_or_default<int>("amr_max_level", 0);
  this->amr_max_preref_cycles_ =
      params.template get_param_value_or_default<int>("max_preref_cycles", 1);
  this->amr_scale_down_factor_ =
      params.template get_param_value_or_default<double>("scale_down_factor", 0.25);
  this->amr_rebalance_ = params.template get_param_value_or_default<bool>("rebalance", false);
}

/**
 * @brief Set the auxiliary variables of the Problem owning this AMR object.
 *
 * @details Called by the Problem before each AMR pass. Among these
 *          variables, those associated with a refinement criterion in their
 *          own Problem (pilot variables of another AMR object on the same
 *          mesh, e.g. the other grains of a partitioned polycrystal) take
 *          part in the derefinement decision (see `DerefineByConsensus()`).
 *          They are not to be confused with the variables of every Problem
 *          of the coupling (`auxvars` of the Step* methods), which are only
 *          resynchronized after a mesh modification.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param coupled_vars Auxiliary variables of the owning Problem.
 */
template <class VAR>
void AMRBase<VAR>::SetCoupledVariables(const std::vector<VAR*>& coupled_vars) {
  this->coupled_vars_ = coupled_vars;
}

/**
 * @brief Check whether a finite element space is built on the mesh managed
 *        by this AMR object.
 *
 * @details The auxiliary variables of a coupling may be defined on other
 *          meshes (e.g. one mesh per grain): they are neither modified by
 *          this object nor relevant for its (de)refinement decisions.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param fes Finite element space to check.
 * @return true if `fes` is built on `par_mesh_`.
 */
template <class VAR>
bool AMRBase<VAR>::IsOnMesh(const mfem::ParFiniteElementSpace* fes) const {
  return fes != nullptr && fes->GetParMesh() == &this->par_mesh_;
}

/**
 * @brief Verify that SetCriteria() has already been called on this AMR
 *        object.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 */
template <class VAR>
void AMRBase<VAR>::EnsureCriteriaSet() const {
  MFEM_VERIFY(this->amr_estimator_ != nullptr,
              "AMR criteria not set: call SetCriteria() on this AMR object before "
              "using Refine()/Derefine().");
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
 */
template <class VAR>
void AMRBase<VAR>::EnsureNCMeshIfNeeded() {
  if (!this->par_mesh_.Nonconforming()) {
    this->par_mesh_.EnsureNCMesh(is_nc_simplices_);
  }
}

/**
 * @brief Resynchronize every variable sharing the mesh with its current
 *        state, after one mesh modification (refinement, derefinement or
 *        rebalancing).
 *
 * @details Variables are gathered from `vars` and `auxvars` without
 *          duplicates (the variables of the Problem owning this AMR object
 *          are usually also part of `auxvars`), keeping only those defined
 *          on the mesh managed by this object, and several variables may
 *          share the same finite element space (same SpatialDiscretization).
 *          The update is thus done in four stages, so that each space is
 *          updated once and its transfer operator stays alive until every
 *          grid function defined on it has been transferred:
 *          1. snapshot of the multi-timestep history of each variable
 *             (`Variable::prepareMeshTransfer()`), before any space update;
 *          2. `Update()` of each distinct finite element space;
 *          3. transfer of each grid function and of its history, with the
 *             same transfer operator (interpolation or rebalancing matrix);
 *          4. `UpdatesFinished()` of each distinct finite element space.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Primary variables sharing the mesh managed by this AMR object.
 * @param auxvars Auxiliary variable collections sharing the same mesh.
 */
template <class VAR>
void AMRBase<VAR>::UpdateVariables(VAR& vars, const std::vector<VAR*>& auxvars) {
  using Var = std::remove_reference_t<decltype(vars[0])>;

  std::vector<Var*> unique_vars;
  auto collect = [this, &unique_vars](VAR& container) {
    for (std::size_t i = 0; i < container.get_variables_number(); i++) {
      Var* var = &container[i];
      // Variables defined on another mesh are not affected by this mesh modification
      if (!this->IsOnMesh(var->get_fespace())) {
        continue;
      }
      if (std::find(unique_vars.begin(), unique_vars.end(), var) == unique_vars.end()) {
        unique_vars.push_back(var);
      }
    }
  };
  collect(vars);
  for (auto* auxvar : auxvars) {
    if (auxvar != nullptr) {
      collect(*auxvar);
    }
  }

  // 1. Snapshots, before any space update
  std::vector<mfem::ParFiniteElementSpace*> fespaces;
  std::vector<std::map<int, mfem::ParGridFunction>> histories(unique_vars.size());
  for (std::size_t k = 0; k < unique_vars.size(); k++) {
    unique_vars[k]->prepareMeshTransfer(histories[k]);
    auto* fes = unique_vars[k]->get_fespace();
    if (std::find(fespaces.begin(), fespaces.end(), fes) == fespaces.end()) {
      fespaces.push_back(fes);
    }
  }
  // 2. One update per space
  for (auto* fes : fespaces) {
    fes->Update();
  }
  // 3. Transfer of the solution and of its history with the same operator
  for (std::size_t k = 0; k < unique_vars.size(); k++) {
    auto* var = unique_vars[k];
    auto& gf = var->get_ref_gf();
    gf.Update();
    var->finalizeMeshTransfer(histories[k]);
    mfem::Vector unk;
    gf.GetTrueDofs(unk);
    var->updateAfterMeshChange(unk);
  }
  // 4. Free the transfer operators
  for (auto* fes : fespaces) {
    fes->UpdatesFinished();
  }
}

/**
 * @brief Rebalance the parallel mesh among the MPI ranks, then
 *        resynchronize every variable sharing it.
 *
 * @details No-op if rebalancing is disabled (`rebalance` parameter), if
 *          the computation runs on a single rank, or if the mesh is
 *          conforming (`mfem::ParMesh::Rebalance()` is only available for
 *          nonconforming meshes). Rebalancing is a mesh modification of
 *          its own: it must be followed by a dedicated update of the
 *          spaces and grid functions (it cannot be merged with the update
 *          following a refinement). The estimator flux spaces do not need
 *          any update here: the estimator updates them itself.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Primary variables sharing the mesh managed by this AMR object.
 * @param auxvars Auxiliary variable collections sharing the same mesh.
 * @return true if the mesh was rebalanced.
 */
template <class VAR>
bool AMRBase<VAR>::Rebalance(VAR& vars, const std::vector<VAR*>& auxvars) {
  if (!this->amr_rebalance_ || this->par_mesh_.GetNRanks() == 1) {
    return false;
  }
  if (!this->par_mesh_.Nonconforming()) {
    SlothInfo::verbose("AMRBase::Rebalance: conforming mesh, rebalancing skipped.");
    return false;
  }
  this->par_mesh_.Rebalance();
  this->UpdateVariables(vars, auxvars);
  // The estimator flux spaces are updated by the estimator itself at the next AMR pass.
  return true;
}

/**
 * @brief Print (on rank 0, Verbose mode) the distribution of the mesh
 *        elements among the MPI ranks: min, max and mean number of
 *        elements per rank, and imbalance factor max/mean.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param stage Label of the AMR stage that has just been applied.
 */
template <class VAR>
void AMRBase<VAR>::PrintLoadBalance(const std::string& stage) const {
  MPI_Comm comm = this->par_mesh_.GetComm();
  const int nranks = this->par_mesh_.GetNRanks();
  const long long local_ne = this->par_mesh_.GetNE();
  long long min_ne = 0;
  long long max_ne = 0;
  long long sum_ne = 0;
  MPI_Allreduce(&local_ne, &min_ne, 1, MPI_LONG_LONG, MPI_MIN, comm);
  MPI_Allreduce(&local_ne, &max_ne, 1, MPI_LONG_LONG, MPI_MAX, comm);
  MPI_Allreduce(&local_ne, &sum_ne, 1, MPI_LONG_LONG, MPI_SUM, comm);
  if (this->par_mesh_.GetMyRank() == 0) {
    const double mean_ne = static_cast<double>(sum_ne) / nranks;
    const double imbalance = mean_ne > 0. ? static_cast<double>(max_ne) / mean_ne : 1.;
    SlothInfo::verbose("AMR load balance (", stage, "): ", sum_ne, " elements on ", nranks,
                       " ranks, min/max/mean = ", min_ne, "/", max_ne, "/", mean_ne,
                       ", imbalance (max/mean) = ", imbalance);
  }
}

/**
 * @brief Apply the initial refinement loop on the initial condition,
 *        before the time-stepping loop starts.
 *
 * @details After each mesh refinement cycle (up to `amr_max_preref_cycles_`):
 *          - `Refine()` decides whether to refine, based on `vars`;
 *          - every variable of `vars` and `auxvars` is resynchronized
 *            (`UpdateVariables()`) — all of them share the same mesh and
 *            must stay consistent with its current state;
 *          - the mesh is rebalanced if requested (`Rebalance()`), so that
 *            the next cycle starts from a balanced partition;
 *          - each variable is then reset to its initial condition via
 *            `Variable::setInitialCondition()`.
 *
 * @tparam VAR Variables container type this AMR strategy operates on.
 * @param vars Primary variables sharing the mesh managed by this AMR
 *            object.
 * @param auxvars Auxiliary variable collections sharing the same mesh,
 *               if any (empty by default).
 */
template <class VAR>
void AMRBase<VAR>::InitialRefine(VAR& vars, std::vector<VAR*> auxvars) {
  this->RegisterPilots(vars);
  int cycle = 0;
  bool did_refine = false;
  while (cycle < this->amr_max_preref_cycles_) {
    did_refine = this->Refine(vars);

    if (!did_refine) {
      break;
    }

    this->UpdateVariables(vars, auxvars);
    const bool did_rebalance = this->Rebalance(vars, auxvars);

    for (std::size_t i = 0; i < vars.get_variables_number(); i++) {
      vars[i].setInitialCondition();
    }
    for (auto* auxvar_container : auxvars) {
      for (std::size_t i = 0; i < auxvar_container->get_variables_number(); i++) {
        (*auxvar_container)[i].setInitialCondition();
      }
    }

    cycle++;
    this->PrintLoadBalance("pre-refinement cycle " + std::to_string(cycle) +
                           (did_rebalance ? ", rebalanced" : ""));
  }
  SlothInfo::verbose("Pre-refinement: ", cycle, " cycles applied on initial condition.");
}

/**
 * @brief Apply one refinement pass during the time-stepping loop.
 *
 * @details - `Refine()` decides based on `vars` only ;
 *          - if the mesh was refined, every variable in `vars` AND
 *            `auxvars` is resynchronized via `UpdateVariables()`, and a
 *            rebalancing is scheduled (see `StepRebalance()`).
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object; drives the refinement decision.
 * @param auxvars Auxiliary variable collections sharing the same mesh;
 *               not consulted for the refinement decision, but always
 *               resynchronized afterwards (empty by default).
 */
template <class VAR>
void AMRBase<VAR>::StepRefine(VAR& vars, std::vector<VAR*> auxvars) {
  bool did_refine = this->Refine(vars);

  if (did_refine) {
    this->UpdateVariables(vars, auxvars);
    this->pending_rebalance_ = true;
  }
}

/**
 * @brief Apply one derefinement pass during the time-stepping loop.
 *
 * @details - `Derefine()` decides based on the pilot variables of this
 *            object and on the auxiliary variables of the owning Problem
 *            that are pilot variables of another AMR object on the same
 *            mesh (see `DerefineByConsensus()`), not on every variable of
 *            `auxvars`;
 *          - if the mesh was derefined, every variable in `vars` AND
 *            `auxvars` is resynchronized via `UpdateVariables()`, and a
 *            rebalancing is scheduled (see `StepRebalance()`).
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Collection of variables sharing the mesh managed by this
 *            AMR object.
 * @param auxvars Variable collections of every Problem of the coupling,
 *               resynchronized after a derefinement (empty by default).
 */
template <class VAR>
void AMRBase<VAR>::StepDerefine(VAR& vars, std::vector<VAR*> auxvars) {
  bool did_derefine = this->Derefine(vars, this->coupled_vars_);

  if (did_derefine) {
    this->UpdateVariables(vars, auxvars);
    this->pending_rebalance_ = true;
  }
}

/**
 * @brief Rebalance the mesh once at the end of an AMR pass, if the mesh
 *        was modified by `StepDerefine()` and/or `StepRefine()` since the
 *        last call.
 *
 * @details A single rebalancing per pass (derefine + refine) avoids two
 *          redistributions of the mesh per time step. The load balance is
 *          printed after each pass that modified the mesh.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param vars Primary variables sharing the mesh managed by this AMR object.
 * @param auxvars Auxiliary variable collections sharing the same mesh.
 */
template <class VAR>
void AMRBase<VAR>::StepRebalance(VAR& vars, std::vector<VAR*> auxvars) {
  if (!this->pending_rebalance_) {
    return;
  }
  this->pending_rebalance_ = false;
  const bool did_rebalance = this->Rebalance(vars, auxvars);
  this->PrintLoadBalance(did_rebalance ? "AMR pass, rebalanced" : "AMR pass");
}

/**
 * @brief Filter refinement candidates that would exceed the configured
 *        maximum refinement depth for their element.
 *
 * @details No-op if amr_max_level_ is 0 (no limit configured, default).
 *          Prevents an element from being refined indefinitely across many
 *          time steps: each StepRefine() call only refines once, so
 *          without this filter an element could accumulate one extra
 *          refinement level per time step with no upper bound.
 *
 * @tparam VAR Variable container type this AMR strategy operates on.
 * @param candidates Refinement candidates before filtering.
 * @return mfem::Array<mfem::Refinement> Filtered candidates, keeping only
 *        elements whose current depth is below amr_max_level_.
 */
template <class VAR>
mfem::Array<mfem::Refinement> AMRBase<VAR>::FilterRefinedCandidates(
    const mfem::Array<mfem::Refinement>& candidates) {
  if (this->amr_max_level_ <= 0) {
    return candidates;
  }

  MFEM_VERIFY(this->par_mesh_.ncmesh != nullptr,
              "FilterRefinedCandidates requires a nonconforming mesh (ncmesh == nullptr).");

  mfem::Array<mfem::Refinement> filtered_candidates;
  for (int i = 0; i < candidates.Size(); i++) {
    int elem = candidates[i].index;
    int depth = this->par_mesh_.ncmesh->GetElementDepth(elem);
    if (depth < this->amr_max_level_) {
      filtered_candidates.Append(candidates[i]);
    }
  }
  return filtered_candidates;
}