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
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <mfem/fem/gridfunc.hpp>
#include <mfem/general/error.hpp>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

////////////////////////////////
// Without attributes names   //
////////////////////////////////
/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param input_file
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const std::tuple<std::string, std::string>& input_file_dir)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(spatial->mesh_, input_file_dir);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const AnalyticalFunctions<DIM>& initial_condition_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_name);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const AnalyticalFunctions<DIM>& initial_condition_name,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_name);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const AnalyticalFunctions<DIM>& initial_condition_name,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);

  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_name);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const double& initial_condition_value)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}
/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const double& initial_condition_value,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const double& initial_condition_value,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/////////////////////////////
// With attributes names   //
/////////////////////////////
/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const AnalyticalFunctions<DIM>& initial_condition_name,
                           const std::set<std::string>& attribute_names)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_name);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const AnalyticalFunctions<DIM>& initial_condition_name,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_name);

  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const AnalyticalFunctions<DIM>& initial_condition_name,
                           const std::set<std::string>& attribute_names,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);

  const auto dim = spatial->get_dimension();

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_name);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::set<std::string>& attribute_names)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);

  this->setAnalyticalSolution(analytical_solution_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::set<std::string>& attribute_names,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const double& initial_condition_value,
                           const std::set<std::string>& attribute_names)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}
/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const double& initial_condition_value,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           GlossaryQuantity type, const int& depth,
                           const double& initial_condition_value,
                           const std::set<std::string>& attribute_names,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name), variable_type_(type) {
  this->fespace_ = spatial->get_finite_element_space();

  this->uh_.SetSpace(fespace_);

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

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
 * @brief build the function associated to initial_condition_name
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param initial_condition_name
 * @return std::function<double(const mfem::Vector&, double)>
 */
template <class T, int DIM>
std::function<double(const mfem::Vector&, double)> Variable<T, DIM>::buildAnalyticalFunction(
    const AnalyticalFunctions<DIM>& analytical_function_name) {
  // this->ics_ = std::make_shared<AnalyticalFunctions<DIM>>();
  std::shared_ptr<AnalyticalFunctions<DIM>> ics(
      new AnalyticalFunctions<DIM>(analytical_function_name));
  return ics->getFunction();
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
 * @brief Define an initial condition on the basis of an analytical function defined by its name
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param initial_condition_name
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(const AnalyticalFunctions<DIM>& initial_condition_name) {
  this->initial_condition_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      this->buildAnalyticalFunction(initial_condition_name));
  mfem::FunctionCoefficient ic_fc(*this->initial_condition_);
  mfem::VectorArrayCoefficient vc(1);
  vc.Set(0, &ic_fc, false);
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }

  this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Reapply this variable's initial condition to its current grid
 *        function, using the analytical function set at construction.
 *
 * @details Used to project the initial condition again after a mesh
 *          change (e.g. during `AMRBase::InitialRefine()`), so that the
 *          interpolated solution on the refined mesh is replaced by a
 *          fresh, exact evaluation of the analytical function instead of
 *          the (possibly less accurate) interpolation of the previous,
 *          coarser-mesh values. Requires `initial_condition_` to have
 *          been set beforehand via the `AnalyticalFunctions<DIM>`
 *          constructor/overload of `setInitialCondition()`; aborts with a
 *          clear error otherwise.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition() {
  MFEM_VERIFY(this->initial_condition_ != nullptr,
              "setInitialCondition() called without argument, but no analytical initial "
              "condition was previously set for variable '" +
                  this->getVariableName() + "'.");

  mfem::FunctionCoefficient ic_fc(*this->initial_condition_);
  mfem::VectorArrayCoefficient vc(1);
  vc.Set(0, &ic_fc, false);
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }
  this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Define an initial condition on the basis of a FunctionCoefficient
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param initial_condition_function
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(
    const mfem::FunctionCoefficient& initial_condition_function) {
  auto ic_fc = initial_condition_function;
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    mfem::VectorArrayCoefficient vc(1);
    vc.Set(0, &ic_fc, false);
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }
  this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Define an initial condition on the basis of a double value
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param initial_condition_value
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(const double& initial_condition_value) {
  mfem::ConstantCoefficient ic_fc(initial_condition_value);
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    mfem::VectorArrayCoefficient vc(1);
    vc.Set(0, &ic_fc, false);
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }
  this->uh_.GetTrueDofs(this->unk_);
}

//////////////////////////////////
//////////////////////////////////
/**
 * @brief Define an analytical solution
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param analytical_solution_name
 */
template <class T, int DIM>
void Variable<T, DIM>::setAnalyticalSolution(
    const AnalyticalFunctions<DIM>& analytical_solution_name) {
  this->analytical_solution_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      this->buildAnalyticalFunction(analytical_solution_name));
}

/**
 * @brief Define an analytical solution
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param analytical_solution_function
 */
template <class T, int DIM>
void Variable<T, DIM>::setAnalyticalSolution(
    const mfem::FunctionCoefficient& analytical_solution_function) {
  this->analytical_solution_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      analytical_solution_function);
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
  return std::prev(this->map_of_unk_.end())->second;
}

/**
 * @brief update the GridFunction on the basis of its associated unknown vector
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::update(const mfem::Vector& unk) {
  this->saveBeforeUpdate();
  this->unk_ = unk;
  this->uh_.SetFromTrueDofs(this->unk_);

  auto current_solution = std::prev(this->map_of_unk_.end());
  current_solution->second = this->unk_;
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
  for (const auto& [key, vec] : this->map_of_unk_) {
    mfem::ParGridFunction gf_h(this->fespace_);
    gf_h.SetFromTrueDofs(vec);
    tmp_history.emplace(key, std::move(gf_h));
  }
}

/**
 * @brief Resynchronize this variable's finite element space, grid function,
 *        unknown vector and multi-timestep history with the current state
 *        of the (possibly just refined/derefined) shared mesh.
 *
 * @details Must be called after any operation that mutates the underlying
 *          mesh (refinement or derefinement), for every variable attached
 *          to it — regardless of whether that specific variable's own
 *          refinement criterion triggered the change, since all variables
 *          sharing the mesh must stay consistent with its current state.
 *
 *          The multi-timestep history (`map_of_unk_`) is remapped using the
 *          exact same transfer operator as the current solution
 *          (`prepareMeshTransfer`/`finalizeMeshTransfer`, called around the
 *          `fespace_`/`gf` update), so that older timesteps stay usable by
 *          multi-step time schemes after the mesh has changed.
 *
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void Variable<T, DIM>::UpdateAndRebalance() {
  mfem::ParFiniteElementSpace* var_fespace_v = this->get_fespace();
  auto& gf = this->get_ref_gf();

  // Prepare a snapshot before update
  std::map<int, mfem::ParGridFunction> tmp_history;
  this->prepareMeshTransfer(tmp_history);
  // Update
  var_fespace_v->Update();
  gf.Update();
  // Interpolation of the snapshot with the same transfer operator
  this->finalizeMeshTransfer(tmp_history);

  mfem::Vector& unk = this->get_ref_unknown();
  gf.GetTrueDofs(unk);

  this->updateAfterMeshChange(unk);

  var_fespace_v->UpdatesFinished();
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
  return this->unk_;
}
template <class T, int DIM>
mfem::Vector& Variable<T, DIM>::get_ref_unknown() {
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
  return this->uh_;
}
template <class T, int DIM>
mfem::ParGridFunction& Variable<T, DIM>::get_ref_gf() {
  return this->uh_;
}

/**
 * @brief return the function associated to the analytical solution
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::function<double(const mfem::Vector&, double)>
 */
template <class T, int DIM>
std::shared_ptr<std::function<double(const mfem::Vector&, double)>>
Variable<T, DIM>::get_analytical_solution() {
  return this->analytical_solution_;
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
  std::shared_ptr<mfem::AttributeSets> elem_attr_sets = spatial->get_elem_attributes();
  std::set<std::string> set_mesh_attr_names = elem_attr_sets->GetAttributeSetNames();
  MFEM_VERIFY(std::includes(set_mesh_attr_names.begin(), set_mesh_attr_names.end(),
                            attribute_names.begin(), attribute_names.end()),
              "Error while setting attributes associated with elements for variable " +
                  this->variable_name_ + ". Please check your data");
  for (const auto& name : attribute_names) {
    this->el_attr_.Append(elem_attr_sets->GetAttributeSet(name));
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
