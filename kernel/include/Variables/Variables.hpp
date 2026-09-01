/**
 * @file Variables.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Variables class used to build and manage variables
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
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Variables/Variable.hpp"

/**
 * @brief Class used to manage a list of Variable
 *
 */
template <class T, int DIM>
class Variables {
 private:
  T* fecollection_;
  std::vector<Variable<T, DIM>> vect_variables_;

 public:
  template <class... Args>
  explicit Variables(Args... args);

  Variables();

  Variable<T, DIM>& operator[](size_t i);
  void add(Variable<T, DIM> var);
  std::vector<Variable<T, DIM>> getVariables() const;
  size_t get_variables_number() const;
  Variable<T, DIM>& get_variable(const std::string& name);
  std::map<std::string, mfem::ParGridFunction*> get_map_gridfunction();
  std::map<std::string, Variable<T, DIM>> get_map_variable() const;
  ~Variables();
};
#include "Variables/Variables.tpp"
