/**
 * @file Coupling.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class to define a  Coupling object
 * @version 0.1
 * @date 2025-09-05
 *
 * @anchor couplings
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
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Utils/UtilsForDebug.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Trait extracting a single VAR type from a pack of Problem
 *        types, used by Coupling.
 *
 * @tparam Args Problem types to extract the shared VAR type from.
 */
template <class... Args>
struct FirstProblemVarType {
  using type = typename std::tuple_element_t<0, std::tuple<Args...>>::VarType;
};

template <class... Args>
class Coupling {
 private:
  using VAR = typename FirstProblemVarType<Args...>::type;

  std::string name_{"Unnamed Coupling"};
  std::tuple<Args...> problems_;

  std::vector<std::tuple<std::string, std::vector<std::tuple<std::string, bool, double>>>>
      pb_convergence_;

 public:
  std::vector<VAR*> CollectAllVariables();
  using VarType = VAR;
  explicit Coupling(const std::string& name, Args... problems);
  std::string get_name();

  void get_tree();
  void initialize(const int& iter, const double& initial_time, const double time_step,
                  std::vector<VAR*> all_vars);
  void execute(const int& iter, double& next_time, const double& current_time,
               const double& current_time_step);
  void update();
  void post_processing(const int& iter, const double& current_time, std::vector<VAR*> all_vars);
  void finalize();

  std::vector<std::tuple<std::string, std::vector<std::tuple<std::string, bool, double>>>>
  get_convergence();
  ~Coupling();
};
#include "Couplings/Coupling.tpp"
