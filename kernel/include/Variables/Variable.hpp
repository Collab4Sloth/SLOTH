/**
 * @file Variable.hpp
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
#include <any>
#include <concepts>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
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

template <class T, int DIM>
class Variable {
 private:
  T* fecollection_;
  BoundaryConditions<T, DIM> bcs_;
  std::string variable_name_;
  GlossaryQuantity variable_type_;
  mfem::ParFiniteElementSpace* fespace_;

  std::map<int, mfem::Vector> map_of_unk_;
  int depth_ = 2;
  mfem::Vector unk_;
  mfem::ParGridFunction uh_;

  // Number of time steps the variable has been advanced (ProblemBase::update), used to know if it
  // has already been advanced in the current step when it is seen as an auxiliary variable by
  // another problem of a partitioned coupling
  int nb_updates_{0};
  // Saved state while an extrapolated state is set (set_extrapolated_state / restore_state)
  bool is_extrapolated_{false};
  std::map<int, mfem::Vector> saved_map_of_unk_;
  mfem::Vector saved_unk_;
  // Physical bounds applied to the extrapolated state only (set_extrapolation_bounds)
  std::optional<std::pair<double, double>> extrapolation_bounds_;

  // Analytical solution: as a std::function (only if given by an AnalyticalFunctions, kept for
  // get_analytical_solution()) and as an mfem::Coefficient (in all cases)
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};
  std::shared_ptr<mfem::Coefficient> analytical_solution_coef_{nullptr};

  // Initial condition as an mfem::Coefficient (nullptr if read from a file (used for save/restart):
  // not reprojectable)
  std::shared_ptr<mfem::Coefficient> ic_coef_{nullptr};
  // Same object as ic_coef_ when the initial condition is given by a SLOTH ::Coefficient
  std::shared_ptr<MfemFromSlothCoefficient> ic_dependent_coef_{nullptr};
  // Names of the variables the initial condition depends on (see set_dependencies())
  std::vector<std::string> dependencies_;
  // True while the initial condition has not been projected yet (deferred projection)
  bool ic_pending_ = false;

  std::vector<std::string> additional_variable_info_;

  mfem::Array<int> el_attr_;

  void setVariableDepth(const int& depth);

  void add_variable_info(const std::string& var);
  void setInitialCondition(SpatialDiscretization<T, DIM>* spatial,
                           const InitialCondition<DIM>& initial_condition);
  void setInitialCondition(mfem::ParMesh* pmesh,
                           const std::tuple<std::string, std::string>& input_file_dir);
  void setAnalyticalSolution(const AnalyticalSolution<DIM>& analytical_solution);
  void projectInitialCondition();
  void projectInitialCondition(const std::vector<const mfem::ParGridFunction*>& dependencies);
  void ensureInitialCondition() const;
  void saveBeforeUpdate();
  void set_attributes(SpatialDiscretization<T, DIM>* spatial,
                      const std::set<std::string>& attribute_names);

 public:
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const InitialCondition<DIM>& initial_condition,
           const std::set<std::string>& attribute_names = {},
           const AnalyticalSolution<DIM>& analytical_solution = {});

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const InitialCondition<DIM>& initial_condition,
           const AnalyticalSolution<DIM>& analytical_solution);

  template <class... Args>
  void set_additional_information(Args&&... add_var_info);

  template <class... Vars>
    requires(sizeof...(Vars) > 0 && (std::same_as<Vars, Variable<T, DIM>> && ...))
  void set_dependencies(Vars&... dependencies);
  std::vector<std::string> get_dependencies() const;

  mfem::Vector get_last() const;
  mfem::Vector get_second_to_last() const;
  std::string getVariableName() const;
  std::vector<std::string> get_additional_variable_info() const;
  void update(const mfem::Vector& unk);

  int get_number_of_updates() const;
  void count_time_update();
  void set_extrapolated_state(const mfem::Vector& u_tilde);
  void restore_state();
  bool is_extrapolated() const;
  void set_extrapolation_bounds(const double lower, const double upper);
  const std::optional<std::pair<double, double>>& get_extrapolation_bounds() const;

  mfem::Vector get_unknown() const;
  mfem::Vector& get_ref_unknown();
  std::map<int, mfem::Vector> get_map_unknown();
  mfem::ParGridFunction get_gf() const;
  mfem::ParGridFunction& get_ref_gf();
  mfem::ParGridFunction get_igf() const;
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> get_analytical_solution();
  std::shared_ptr<mfem::Coefficient> get_analytical_solution_coefficient();
  BoundaryConditions<T, DIM>* get_boundary_conditions();
  mfem::ParFiniteElementSpace* get_fespace();

  bool isInitialConditionReprojectable() const;
  void setInitialCondition();
  void setInitialCondition(const std::map<std::string, mfem::ParGridFunction*>& gridfunctions);
  void updateAfterMeshChange(const mfem::Vector& unk);
  void prepareMeshTransfer(std::map<int, mfem::ParGridFunction>& tmp_history);
  void finalizeMeshTransfer(std::map<int, mfem::ParGridFunction>& tmp_history);

  /**
   * @brief State of a Variable: current unknown, multi-timestep history and number of updates.
   *        Used to restore the Variable when a time-step is rejected.
   */
  struct State {
    mfem::Vector unk;
    std::map<int, mfem::Vector> history;
    int nb_updates{0};
  };
  State save_state() const;
  void restore_state(const State& state);

  ~Variable();
};
#include "Variables/Variable.tpp"
