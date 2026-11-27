/**
 * @file Variable.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Variable class used to build and manage variables of model
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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
#include <any>
#include <concepts>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/Coefficient.hpp"
#include "Coefficients/FileCoefficient.hpp"
#include "Coefficients/MfemFromSlothCoefficient.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/VariableSolution.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @details The initial condition can be given as:
 *          - a double;
 *          - an AnalyticalFunctions<DIM>;
 *          - an mfem::FunctionCoefficient;
 *          - a SLOTH ::Coefficient (constant, file-based, or depending on other variables:
 *            see set_dependencies());
 *          - a (file name, directory) tuple, for a MFEM grid function file read on each rank.
 *          The analytical solution, if any, can be given as a double, an AnalyticalFunctions<DIM>,
 *          an mfem::FunctionCoefficient or a SLOTH ::Coefficient (never depending on other
 *          variables).
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param spatial Spatial discretization
 * @param bcs Boundary conditions
 * @param variable_name Name of the variable
 * @param type Quantity associated with the variable
 * @param depth Number of stored time levels (at least 2)
 * @param initial_condition Initial condition
 * @param attribute_names Names of the element attribute sets where the initial condition is
 *                        projected (whole domain if empty)
 * @param analytical_solution Analytical solution (none by default)
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const InitialCondition<DIM>& initial_condition,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalSolution<DIM>& analytical_solution)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  if (!attribute_names.empty()) {
    this->set_attributes(spatial, attribute_names);
  }
  this->setInitialCondition(spatial, initial_condition);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object, with an analytical solution but without
 *        attribute names
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param spatial Spatial discretization
 * @param bcs Boundary conditions
 * @param variable_name Name of the variable
 * @param type Quantity associated with the variable
 * @param depth Number of stored time levels (at least 2)
 * @param initial_condition Initial condition
 * @param analytical_solution Analytical solution
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const InitialCondition<DIM>& initial_condition,
                           const AnalyticalSolution<DIM>& analytical_solution)
    : Variable(spatial, bcs, variable_name, type, depth, initial_condition, std::set<std::string>{},
               analytical_solution) {}

/**
 * @brief Associate additional information to the variable
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @tparam Args
 * @param add_var_info
 */
template <class T, int DIM>
template <class... Args>
void Variable<T, DIM>::set_additional_information(Args&&... add_var_info) {
  if constexpr (sizeof...(add_var_info) == 0) {
    this->additional_variable_info_.resize(0);
  } else {
    (add_variable_info(std::forward<Args>(add_var_info)), ...);
  }
}

/**
 * @brief Define the variables the initial condition depends on, and project the initial
 *        condition with their current values.
 *
 * @details Only meaningful when the initial condition is given by a SLOTH ::Coefficient that is
 *          neither a constant nor file-based. The values of the dependencies at each point are
 *          forwarded to the coefficient, in the order of the arguments.
 *
 *          Must be called before any use of the variable, in particular before adding it to a
 *          Variables container (which stores a copy). Without call to set_dependencies(), such
 *          an initial condition is projected without dependency at the first use of the
 *          variable.
 *
 *          Only the names of the dependencies are kept: after a mesh change, the initial
 *          condition is reprojected with setInitialCondition(gridfunctions).
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @tparam Vars Variable<T, DIM>
 * @param dependencies Variables the initial condition depends on
 */
template <class T, int DIM>
template <class... Vars>
  requires(sizeof...(Vars) > 0 && (std::same_as<Vars, Variable<T, DIM>> && ...))
void Variable<T, DIM>::set_dependencies(Vars&... dependencies) {
  MFEM_VERIFY(
      this->ic_dependent_coef_ != nullptr && this->ic_dependent_coef_->requires_dependencies(),
      "set_dependencies(): the initial condition of variable '" + this->variable_name_ +
          "' cannot depend on other variables (only a SLOTH Coefficient that is neither "
          "constant nor file-based can).");
  MFEM_VERIFY(this->ic_pending_,
              "set_dependencies(): the initial condition of variable '" + this->variable_name_ +
                  "' has already been projected. set_dependencies() must be called before any "
                  "use of the variable, in particular before adding it to a Variables container.");

  std::vector<std::string> names{dependencies.getVariableName()...};
  for (const auto& name : names) {
    MFEM_VERIFY(
        name != this->variable_name_,
        "set_dependencies(): variable '" + this->variable_name_ + "' cannot depend on itself.");
  }

  // get_ref_gf() projects the initial condition of the dependencies if not done yet
  std::vector<const mfem::ParGridFunction*> gridfunctions{&dependencies.get_ref_gf()...};

  this->dependencies_ = std::move(names);
  this->projectInitialCondition(gridfunctions);
}

/**
 * @brief Return the names of the variables the initial condition depends on
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::vector<std::string>
 */
template <class T, int DIM>
std::vector<std::string> Variable<T, DIM>::get_dependencies() const {
  return this->dependencies_;
}

//////////////////////////////////
// Set initial conditions
//////////////////////////////////

/**
 * @brief Define an initial condition on the basis of an input file
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param input_file
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(
    mfem::ParMesh* pmesh, const std::tuple<std::string, std::string>& input_file_dir) {
  auto [input_file, input_dir] = input_file_dir;
  int rank = mfem::Mpi::WorldRank();

  std::ostringstream oss;
  oss << std::setfill('0') << std::setw(6) << rank;
  std::filesystem::path inputfile_rank =
      std::filesystem::path(input_dir) / (input_file + "." + oss.str());

  std::ifstream file_stream(inputfile_rank);

  if (!file_stream.is_open()) {
    const std::string& mess = "Cannot open input file: " + inputfile_rank.string();
    mfem::mfem_error(mess.c_str());
  }
  mfem::GridFunction file_gf(pmesh, file_stream);

  this->uh_ = std::move(file_gf);

  this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Define the initial condition
 *
 * @details A SLOTH ::Coefficient that may depend on other variables is not projected here: its
 *          projection is deferred to set_dependencies() or, failing that, to the first use of
 *          the variable.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param spatial Spatial discretization
 * @param initial_condition Initial condition
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(SpatialDiscretization<T, DIM>* spatial,
                                           const InitialCondition<DIM>& initial_condition) {
  if (initial_condition.is_input_file()) {
    MFEM_VERIFY(this->el_attr_.Size() == 0,
                "Variable '" + this->variable_name_ +
                    "': attribute names cannot be used with an initial condition read from a "
                    "file.");
    this->setInitialCondition(spatial->mesh_, initial_condition.get_input_file_dir());
    return;
  }

  this->ic_coef_ = initial_condition.SlothToMfemCoefficient();
  this->ic_dependent_coef_ = std::dynamic_pointer_cast<MfemFromSlothCoefficient>(this->ic_coef_);

  if (initial_condition.requires_dependencies()) {
    this->ic_pending_ = true;
    return;
  }
  this->projectInitialCondition();
}

/**
 * @brief Project the initial condition onto the grid function, without dependency
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::projectInitialCondition() {
  this->projectInitialCondition(std::vector<const mfem::ParGridFunction*>{});
}

/**
 * @brief Project the initial condition onto the grid function, and reset the multi-timestep
 *        history to it
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param dependencies Grid functions of the variables the initial condition depends on, in the
 *                     order of get_dependencies() (empty if none)
 */
template <class T, int DIM>
void Variable<T, DIM>::projectInitialCondition(
    const std::vector<const mfem::ParGridFunction*>& dependencies) {
  MFEM_VERIFY(this->ic_coef_ != nullptr,
              "Variable '" + this->variable_name_ + "': no initial condition to project.");
  MFEM_VERIFY(dependencies.empty() || this->ic_dependent_coef_ != nullptr,
              "Variable '" + this->variable_name_ +
                  "': the initial condition cannot depend on other variables.");
  for (const auto* gf : dependencies) {
    MFEM_VERIFY(gf != nullptr && gf->ParFESpace()->GetParMesh() == this->fespace_->GetParMesh(),
                "Variable '" + this->variable_name_ +
                    "': the variables the initial condition depends on must share its mesh.");
  }

  if (this->ic_dependent_coef_ != nullptr) {
    this->ic_dependent_coef_->bind(dependencies);
  }

  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    mfem::VectorArrayCoefficient vc(1);
    vc.Set(0, this->ic_coef_.get(), false);
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(*this->ic_coef_);
  }
  this->uh_.GetTrueDofs(this->unk_);

  // Only non-owning pointers: never keep them beyond the projection
  if (this->ic_dependent_coef_ != nullptr) {
    this->ic_dependent_coef_->unbind();
  }

  for (auto& [key, vec] : this->map_of_unk_) {
    vec = this->unk_;
  }
  this->ic_pending_ = false;
}

/**
 * @brief Project the initial condition if its projection has been deferred (initial condition
 *        given by a SLOTH ::Coefficient, without call to set_dependencies())
 *
 * @remark const so that it can be called from the const accessors: the projection is part of
 *         the construction of the variable, deferred.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::ensureInitialCondition() const {
  if (this->ic_pending_) {
    const_cast<Variable<T, DIM>*>(this)->projectInitialCondition();
  }
}

/**
 * @brief True if the initial condition can be projected again, e.g. after a mesh change
 *        (false for an initial condition read from a file)
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
bool Variable<T, DIM>::isInitialConditionReprojectable() const {
  return this->ic_coef_ != nullptr;
}

/**
 * @brief Reapply this variable's initial condition to its current grid
 *        function.
 *
 * @details Used to project the initial condition again after a mesh
 *          change (e.g. during `AMRBase::InitialRefine()`), so that the
 *          interpolated solution on the refined mesh is replaced by a
 *          fresh, exact evaluation of the initial condition instead of
 *          the (possibly less accurate) interpolation of the previous,
 *          coarser-mesh values. The multi-timestep history is reset to it,
 *          which is why it can only be called before the first time update.
 *          Aborts with a clear error if the initial condition is not
 *          reprojectable (read from a file) or depends on other variables
 *          (use `setInitialCondition(gridfunctions)` instead).
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition() {
  MFEM_VERIFY(this->nb_updates_ == 0,
              "setInitialCondition(): variable '" + this->variable_name_ +
                  "' has already been advanced in time: its history cannot be reset.");
  MFEM_VERIFY(this->isInitialConditionReprojectable(),
              "setInitialCondition() called without argument, but the initial condition of "
              "variable '" +
                  this->getVariableName() + "' cannot be reprojected (read from a file).");
  MFEM_VERIFY(this->dependencies_.empty(),
              "setInitialCondition() called without argument, but the initial condition of "
              "variable '" +
                  this->getVariableName() +
                  "' depends on other variables: use setInitialCondition(gridfunctions).");
  this->projectInitialCondition();
}

/**
 * @brief Reapply this variable's initial condition to its current grid function, looking up
 *        the variables it depends on (if any) by name.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param gridfunctions Grid functions of the available variables, by name. Must contain all the
 *                      dependencies, already up to date on the current mesh.
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(
    const std::map<std::string, mfem::ParGridFunction*>& gridfunctions) {
  MFEM_VERIFY(this->nb_updates_ == 0,
              "setInitialCondition(): variable '" + this->variable_name_ +
                  "' has already been advanced in time: its history cannot be reset.");
  MFEM_VERIFY(this->isInitialConditionReprojectable(),
              "setInitialCondition(): the initial condition of variable '" +
                  this->getVariableName() + "' cannot be reprojected (read from a file).");

  std::vector<const mfem::ParGridFunction*> dependencies;
  dependencies.reserve(this->dependencies_.size());
  for (const auto& name : this->dependencies_) {
    const auto it = gridfunctions.find(name);
    MFEM_VERIFY(it != gridfunctions.end(), "setInitialCondition(): variable '" + name +
                                               "', on which the initial condition "
                                               "of variable '" +
                                               this->getVariableName() +
                                               "' depends, is not available.");
    dependencies.push_back(it->second);
  }
  this->projectInitialCondition(dependencies);
}

/**
 * @brief Define the analytical solution
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param analytical_solution Analytical solution (possibly undefined)
 */
template <class T, int DIM>
void Variable<T, DIM>::setAnalyticalSolution(const AnalyticalSolution<DIM>& analytical_solution) {
  this->analytical_solution_ = analytical_solution.get_function();
  this->analytical_solution_coef_ = analytical_solution.SlothToMfemCoefficient();
}

/**
 * @brief Return the second-to-last term of the saved variables
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::Vector
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_second_to_last() const {
  this->ensureInitialCondition();
  return std::prev(std::prev(this->map_of_unk_.end()))->second;
}

/**
 * @brief Return the last term of the saved variables
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::Vector
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_last() const {
  this->ensureInitialCondition();
  return std::prev(this->map_of_unk_.end())->second;
}

/**
 * @brief Return a copy of the state of the Variable (current unknown, history, number of updates)
 *
 * @details The copy refers to the current finite element space: it can be restored only if the
 *          mesh has not changed in between.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return State
 */
template <class T, int DIM>
typename Variable<T, DIM>::State Variable<T, DIM>::save_state() const {
  this->ensureInitialCondition();
  return State{this->unk_, this->map_of_unk_, this->nb_updates_};
}

/**
 * @brief Restore a state returned by save_state() (current unknown, grid function, history and
 *        number of updates)
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param state State to restore.
 */
template <class T, int DIM>
void Variable<T, DIM>::restore_state(const State& state) {
  this->ensureInitialCondition();
  MFEM_VERIFY(state.unk.Size() == this->unk_.Size(),
              "Variable::restore_state: the size of the saved state differs from the current one.");
  this->unk_ = state.unk;
  this->uh_.SetFromTrueDofs(this->unk_);
  this->map_of_unk_ = state.history;
  this->nb_updates_ = state.nb_updates;
}

/**
 * @brief update the GridFunction on the basis of its associated unknown vector
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::update(const mfem::Vector& unk) {
  this->ensureInitialCondition();
  MFEM_VERIFY(!this->is_extrapolated_,
              "Variable::update: " << this->variable_name_
                                   << " is in an extrapolated state (restore_state expected).");
  this->saveBeforeUpdate();
  this->unk_ = unk;
  this->uh_.SetFromTrueDofs(this->unk_);

  auto current_solution = std::prev(this->map_of_unk_.end());
  current_solution->second = this->unk_;
}

/**
 * @brief Temporarily set an extrapolated state at the end of the current step.
 *
 * @details The current value becomes u~ and the history is shifted as for a time update
 *          (second-to-last = current value), so that an integrator reading {last,
 *          second-to-last} sees {u~, u^n}. The previous state is saved and must be
 *          restored by restore_state() before any time update, mesh change or
 *          post-processing.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param u_tilde Extrapolated state (true dofs)
 */
template <class T, int DIM>
void Variable<T, DIM>::set_extrapolated_state(const mfem::Vector& u_tilde) {
  this->ensureInitialCondition();
  MFEM_VERIFY(!this->is_extrapolated_,
              "Variable::set_extrapolated_state: " << this->variable_name_
                                                   << " is already in an extrapolated state.");
  MFEM_VERIFY(u_tilde.Size() == this->unk_.Size(),
              "Variable::set_extrapolated_state: wrong size for " << this->variable_name_);
  this->saved_map_of_unk_ = this->map_of_unk_;
  this->saved_unk_ = this->unk_;
  this->is_extrapolated_ = true;

  this->saveBeforeUpdate();
  this->unk_ = u_tilde;
  // Truncation in the physical bounds (projection: never increases the extrapolation error)
  if (this->extrapolation_bounds_.has_value()) {
    const auto [lower, upper] = *this->extrapolation_bounds_;
    mfem::real_t* h_unk = this->unk_.HostReadWrite();
    for (int i = 0; i < this->unk_.Size(); i++) {
      h_unk[i] = std::min(std::max(h_unk[i], lower), upper);
    }
  }
  this->uh_.SetFromTrueDofs(this->unk_);
  std::prev(this->map_of_unk_.end())->second = this->unk_;
}

/**
 * @brief Physical bounds applied only to the extrapolated state
 *        (see set_extrapolated_state)
 * @details   Physical bounds [lower, upper] of the variable, applied ONLY to its extrapolated state
 * (BDF2, lagged auxiliary variable of a partitioned coupling): the extrapolation is truncated  node
 * by node in [lower, upper]. Since the exact solution lies in [lower, upper], the truncation never
 * increases the extrapolation error (order 2 kept). The variable itself (u^{n+1}, outputs) is never
 * truncated.
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param lower lower bound
 * @param upper upper bound
 */
template <class T, int DIM>
void Variable<T, DIM>::set_extrapolation_bounds(const double lower, const double upper) {
  MFEM_VERIFY(lower < upper,
              "Variable::set_extrapolation_bounds: lower bound must be lower than "
              "the upper bound for "
                  << this->variable_name_);
  this->extrapolation_bounds_ = std::make_pair(lower, upper);
}

/**
 * @brief Restore the state saved by set_extrapolated_state (no effect otherwise)
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::restore_state() {
  if (!this->is_extrapolated_) return;
  this->map_of_unk_ = this->saved_map_of_unk_;
  this->unk_ = this->saved_unk_;
  this->uh_.SetFromTrueDofs(this->unk_);
  this->saved_map_of_unk_.clear();
  this->saved_unk_.Destroy();
  this->is_extrapolated_ = false;
}

/**
 * @brief Return the number of time steps the variable has been advanced (counted by
 * ProblemBase::update, not by the initialization nor by a mesh change)
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return int
 */
template <class T, int DIM>
int Variable<T, DIM>::get_number_of_updates() const {
  return this->nb_updates_;
}

/**
 * @brief Count one time update of the variable (called by ProblemBase::update)
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::count_time_update() {
  ++this->nb_updates_;
}

/**
 * @brief Check if variable is extrapolated
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return true
 * @return false
 */
template <class T, int DIM>
bool Variable<T, DIM>::is_extrapolated() const {
  return this->is_extrapolated_;
}

/**
 * @brief Get the extrapolation bounds object
 *
 * @tparam T
 * @tparam DIM
 * @return const Variable<T, DIM>::std::optional<std::pair<double, double>>&
 */
template <class T, int DIM>
const std::optional<std::pair<double, double>>& Variable<T, DIM>::get_extrapolation_bounds() const {
  return this->extrapolation_bounds_;
}

/**
 * @brief Snapshot the multi-timestep history into temporary grid functions,
 *        on the finite element space as it exists BEFORE the mesh changes.
 *
 * @details Must be called before `fespace_->Update()`/`gf.Update()` are
 *          invoked (i.e. before the mesh/fespace transfer operator is
 *          built), and paired with `finalizeMeshTransfer()` afterwards.
 *          Each historical unknown vector in `map_of_unk_` is projected
 *          into a `mfem::ParGridFunction` on the current (pre-refinement)
 *          finite element space, so that it can later be interpolated onto
 *          the new mesh state using the exact same transfer operator as
 *          the current solution.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param tmp_history Output map, keyed the same way as `map_of_unk_`,
 *                    populated with one grid function per historical
 *                    timestep entry.
 */
template <class T, int DIM>
void Variable<T, DIM>::prepareMeshTransfer(std::map<int, mfem::ParGridFunction>& tmp_history) {
  this->ensureInitialCondition();
  for (const auto& [key, vec] : this->map_of_unk_) {
    mfem::ParGridFunction gf_h(this->fespace_);
    gf_h.SetFromTrueDofs(vec);
    tmp_history.emplace(key, std::move(gf_h));
  }
}

/**
 * @brief Interpolate the snapshotted multi-timestep history onto the new
 *        mesh state, and write the result back into `map_of_unk_`.
 *
 * @details Must be called after `fespace_->Update()`/`gf.Update()` (while
 *          the mesh's transfer operator built by that update is still
 *          valid, i.e. before `UpdatesFinished()` is called on the
 *          finite element space), using the grid functions previously
 *          produced by `prepareMeshTransfer()`. Each grid function's
 *          `Update()` reuses the same transfer operator that was just
 *          used to interpolate the current solution, ensuring the history
 *          is remapped consistently with it.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param tmp_history Map of temporary grid functions produced by
 *                    `prepareMeshTransfer()`, one per historical timestep.
 */
template <class T, int DIM>
void Variable<T, DIM>::finalizeMeshTransfer(std::map<int, mfem::ParGridFunction>& tmp_history) {
  for (auto& [key, gf_h] : tmp_history) {
    gf_h.Update();
    mfem::Vector new_vec;
    gf_h.GetTrueDofs(new_vec);
    this->map_of_unk_.at(key) = new_vec;
  }
}

/**
 * @brief Shift the multi-timestep history one slot forward, dropping the
 *        oldest entry and freeing the last slot for the new solution.
 *
 * @details Called at the start of a normal time-stepping update (see
 *          `update()`), NOT after an AMR-triggered mesh change (see
 *          `updateAfterMeshChange()`, which intentionally does not shift
 *          the history). Copies `map_of_unk_[key]` into `map_of_unk_[key-1]`
 *          for every key except the first, so that after this call the
 *          last slot still holds the most recent solution and is ready to
 *          be overwritten by the new one.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::saveBeforeUpdate() {
  auto begin = this->map_of_unk_.begin();
  auto end = this->map_of_unk_.end();
  ++begin;
  for (auto it = begin; it != end; ++it) {
    auto itm = (it->first) - 1;
    this->map_of_unk_.at(itm) = it->second;
  }
}

/**
 * @brief Update the current unknown vector and grid function after a mesh
 *        change (AMR refinement/derefinement), without shifting the
 *        multi-timestep history.
 *
 * @details Distinct from `update()`, which is used for normal time
 *          advancement and calls `saveBeforeUpdate()` to shift the
 *          history. This method must be used instead when `unk` results
 *          from a mesh transfer (interpolation following AMR) rather than
 *          from solving a new time step: the history itself is remapped
 *          separately, via `prepareMeshTransfer()`/`finalizeMeshTransfer()`,
 *          called around the same mesh change. Calling `update()` (with its
 *          history shift) instead of this method after an AMR transfer
 *          would incorrectly treat the mesh change as a new time step.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param unk New unknown vector (true DOFs), already interpolated onto the
 *           current (post-mesh-change) finite element space.
 */
template <class T, int DIM>
void Variable<T, DIM>::updateAfterMeshChange(const mfem::Vector& unk) {
  this->unk_ = unk;
  this->uh_.SetFromTrueDofs(this->unk_);
}

/**
 * @brief return the unkown vector
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::Vector
 *
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_unknown() const {
  this->ensureInitialCondition();
  return this->unk_;
}
template <class T, int DIM>
mfem::Vector& Variable<T, DIM>::get_ref_unknown() {
  this->ensureInitialCondition();
  return this->unk_;
}

/**
 * @brief Return the additionnal information associated to variable
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::vector<std::string>
 */
template <class T, int DIM>
std::vector<std::string> Variable<T, DIM>::get_additional_variable_info() const {
  return this->additional_variable_info_;
}

/**
 * @brief Return a map of unknows
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::map<int, mfem::Vector>
 */
template <class T, int DIM>
std::map<int, mfem::Vector> Variable<T, DIM>::get_map_unknown() {
  this->ensureInitialCondition();
  return this->map_of_unk_;
}

/**
 * @brief return the gridfunction associated to the unknown
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::ParGridFunction
 */
template <class T, int DIM>
mfem::ParGridFunction Variable<T, DIM>::get_gf() const {
  this->ensureInitialCondition();
  return this->uh_;
}
template <class T, int DIM>
mfem::ParGridFunction& Variable<T, DIM>::get_ref_gf() {
  this->ensureInitialCondition();
  return this->uh_;
}

/**
 * @brief return the function associated to the analytical solution
 *
 * @remark Only available when the analytical solution is given by an AnalyticalFunctions:
 *         use get_analytical_solution_coefficient() for the general case.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::function<double(const mfem::Vector&, double)> nullptr otherwise
 */
template <class T, int DIM>
std::shared_ptr<std::function<double(const mfem::Vector&, double)>>
Variable<T, DIM>::get_analytical_solution() {
  return this->analytical_solution_;
}

/**
 * @brief return the analytical solution as an mfem::Coefficient, whatever the way it was given
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::shared_ptr<mfem::Coefficient> nullptr if no analytical solution was given
 */
template <class T, int DIM>
std::shared_ptr<mfem::Coefficient> Variable<T, DIM>::get_analytical_solution_coefficient() {
  return this->analytical_solution_coef_;
}

/**
 * @brief return the boundary condition object associated to the variable
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return BoundaryConditions<T, DIM>
 */
template <class T, int DIM>
BoundaryConditions<T, DIM>* Variable<T, DIM>::get_boundary_conditions() {
  return &this->bcs_;
}

/**
 * @brief Get the name of the Variable
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::string name of the variable
 */
template <class T, int DIM>
std::string Variable<T, DIM>::getVariableName() const {
  return this->variable_name_;
}

/**
 * @brief Set the variable depth and initialize at the given initial condition
 * @remark By default, 2 levels are considered. Not optimal in term of memory?
 * @remark By default, 3 levels in case of adaptive time-step?
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param depth
 */
template <class T, int DIM>
void Variable<T, DIM>::setVariableDepth(const int& depth) {
  this->depth_ = std::max(2, depth);

  for (auto id = 0; id < this->depth_; id++) {
    this->map_of_unk_.insert(std::pair<int, mfem::Vector>(id, this->unk_));
  }
}

/**
 * @brief
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param var
 */
template <class T, int DIM>
void Variable<T, DIM>::add_variable_info(const std::string& var) {
  additional_variable_info_.emplace_back(var);
}

/**
 * @brief Return the pointer towards the FiniteElementSpace
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
mfem::ParFiniteElementSpace* Variable<T, DIM>::get_fespace() {
  return this->fespace_;
}

/**
 * @brief set the attributes for the variable
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param spatial
 */
template <class T, int DIM>
void Variable<T, DIM>::set_attributes(SpatialDiscretization<T, DIM>* spatial,
                                      const std::set<std::string>& attribute_names) {
  const mfem::AttributeSets& elem_attr_sets = spatial->get_elem_attribute_sets();
  std::set<std::string> set_mesh_attr_names = elem_attr_sets.GetAttributeSetNames();
  MFEM_VERIFY(std::includes(set_mesh_attr_names.begin(), set_mesh_attr_names.end(),
                            attribute_names.begin(), attribute_names.end()),
              "Error while setting attributes associated with elements for variable " +
                  this->variable_name_ + ". Please check your data");
  for (const auto& name : attribute_names) {
    this->el_attr_.Append(elem_attr_sets.GetAttributeSet(name));
  }
}

/**
 * @brief Destroy the Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
Variable<T, DIM>::~Variable() {}
