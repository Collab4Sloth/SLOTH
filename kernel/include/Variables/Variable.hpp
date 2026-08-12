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
#include <algorithm>
#include <any>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class T, int DIM>
class Variable {
 private:
  T* fecollection_;
  BoundaryConditions<T, DIM> bcs_;
  std::string variable_name_;
  GlossaryQuantity variable_type_;
  mfem::ParFiniteElementSpace* fespace_;

  // std::shared_ptr<AnalyticalFunctions<DIM>> ics_;
  std::map<int, mfem::Vector> map_of_unk_;
  int depth_ = 2;
  mfem::Vector unk_;
  mfem::ParGridFunction uh_;

  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};

  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> initial_condition_{nullptr};

  std::vector<std::string> additional_variable_info_;

  mfem::Array<int> el_attr_;

  void setVariableDepth(const int& depth);

  void add_variable_info(const std::string& var);
  std::function<double(const mfem::Vector&, double)> buildAnalyticalFunction(
      const AnalyticalFunctions<DIM>& analytical_function);
  void setInitialCondition(const AnalyticalFunctions<DIM>& initial_condition_name);
  void setInitialCondition(const mfem::FunctionCoefficient& initial_condition_function);
  void setInitialCondition(const double& initial_condition_value);
  void setInitialCondition(mfem::ParMesh* pmesh,
                           const std::tuple<std::string, std::string>& input_file_dir);

  void setAnalyticalSolution(const AnalyticalFunctions<DIM>& analytical_solution_name);
  void setAnalyticalSolution(const mfem::FunctionCoefficient& analytical_solution_function);
  void saveBeforeUpdate();
  void set_attributes(SpatialDiscretization<T, DIM>* spatial,
                      const std::set<std::string>& attribute_names);

 public:
  /////////////////////////////
  // Without attributes names   //
  /////////////////////////////
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const std::tuple<std::string, std::string>& input_file_dir);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const double& initial_condition_value);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const double& initial_condition_value,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const double& initial_condition_value,
           const mfem::FunctionCoefficient& analytical_solution_function);

  /////////////////////////////
  // With attributes names   //
  /////////////////////////////
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const std::set<std::string>& attribute_names);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const std::set<std::string>& attribute_names,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const std::set<std::string>& attribute_names,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::set<std::string>& attribute_names);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::set<std::string>& attribute_names,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::set<std::string>& attribute_names,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const double& initial_condition_value, const std::set<std::string>& attribute_names);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const double& initial_condition_value, const std::set<std::string>& attribute_names,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, GlossaryQuantity type, const int& depth,
           const double& initial_condition_value, const std::set<std::string>& attribute_names,
           const mfem::FunctionCoefficient& analytical_solution_function);

  template <class... Args>
  void set_additional_information(Args&&... add_var_info);

  mfem::Vector get_last() const;
  mfem::Vector get_second_to_last() const;
  std::string getVariableName() const;
  std::vector<std::string> get_additional_variable_info() const;
  // std::shared_ptr<AnalyticalFunctions<DIM>> getInitialCondition();
  void update(const mfem::Vector& unk);
  mfem::Vector get_unknown() const;
  mfem::Vector& get_ref_unknown();
  std::map<int, mfem::Vector> get_map_unknown();
  mfem::ParGridFunction get_gf() const;
  mfem::ParGridFunction& get_ref_gf();
  mfem::ParGridFunction get_igf() const;
  // mfem::ParGridFunction get_analytical_solution();
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> get_analytical_solution();
  BoundaryConditions<T, DIM>* get_boundary_conditions();
  mfem::ParFiniteElementSpace* get_fespace();

  void setInitialCondition();
  void UpdateAndRebalance();
  void updateAfterMeshChange(const mfem::Vector& unk);
  void prepareMeshTransfer(std::map<int, mfem::ParGridFunction>& tmp_history);
  void finalizeMeshTransfer(std::map<int, mfem::ParGridFunction>& tmp_history);

  ~Variable();
};
#include "Variables/Variable.tpp"
