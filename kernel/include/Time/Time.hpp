/**
 * @file Time.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Main time loop
 * @version 0.1
 * @date 2025-09-05
 *
 * @anchor time
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

#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Problems/Problem.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Trait extracting a single VAR type from a pack of Coupling
 *        types, used by TimeDiscretization.
 *
 * @tparam Args Coupling types to extract the shared VAR type from.
 */
template <class... Args>
struct FirstCouplingVarType {
  using type = typename std::tuple_element_t<0, std::tuple<Args...>>::VarType;
};

template <class... Args>
class TimeDiscretization {
 private:
  using VAR = typename FirstCouplingVarType<Args...>::type;

  const Parameters& params_;
  std::tuple<Args...> couplings_;
  double initial_iter_{0};
  double initial_time_{0.};
  double final_time_;
  double time_step_;
  double current_time_;
  double current_time_step_;
  bool last_step_{false};

  std::vector<std::tuple<
      std::string,
      std::vector<std::tuple<std::string, std::vector<std::tuple<std::string, bool, double>>>>>>
      convergence_;
  bool checkTimeStepConvergence(const auto& convergence);
  void get_parameters();

  void check_data_before_execute();

  void initialize();
  void execute(const int& iter);
  void update();
  void post_processing(const int& iter);
  void finalize();
  void time_management();
  void time_info(const int& iter);
  void check_duplicate_coupling_name();
  void check_duplicate_post_processing_directory();
  std::function<double(double)> time_step_function_;
  std::vector<VAR*> CollectAllVariables();

 public:
  // explicit TimeDiscretization(const Parameters& params, Args&&... couplings);
  explicit TimeDiscretization(const Parameters& params, Args... couplings);
  explicit TimeDiscretization(const std::function<double(double)>& given_time_step,
                              const Parameters& params, Args... couplings);
  void solve();
  void get_tree();
  ~TimeDiscretization();
};

#include "Time/Time.tpp"
