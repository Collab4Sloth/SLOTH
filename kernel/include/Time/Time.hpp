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

#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Problems/Problem.hpp"
#include "Time/TimeStepFailure.hpp"
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
  bool vtk_unified_{false};
  void save_vtk_unified(const int& iter, const double& current_time);

  //////////////////////////////////////////////////////////
  // Adaptive time-stepping (parameter "adaptive_time_step")
  //////////////////////////////////////////////////////////
  // A time-step is rejected if it fails (TimeStepFailure, e.g. Newton not converged): the state
  // of every coupling (variables, internal state of the time schemes) and the time are restored,
  // and the time-step is retried with a smaller time-step. After an accepted time-step, the next
  // time-step depends on the largest number of Newton iterations over all the problems; the
  // ratio between two successive time-steps is limited (time_step_max_ratio, BDF2).
  bool adaptive_{false};
  double time_step_min_{0.};
  double time_step_max_{std::numeric_limits<double>::infinity()};
  // Time-step multiplied by this factor after a rejection
  double cutback_factor_{0.5};
  // Time-step multiplied (divided) by this factor if the number of Newton iterations is lower
  // (greater) than or equal to newton_iterations_low_ (newton_iterations_high_)
  double growth_factor_{1.25};
  int newton_iterations_low_{3};
  int newton_iterations_high_{8};
  // Maximum number of successive rejections of the same time-step
  int max_rejections_{10};
  // Number of successive reductions of the time-step (rejections or Newton criterion), reset
  // when a time-step is accepted without reduction
  int successive_cuts_{0};

  void get_adaptive_parameters();
  void save_state();
  void restore_state();
  std::optional<std::string> try_time_step(const int& iter);
  void reject_time_step(const std::string& reason);
  void cut_time_step(const double& dt_new, const std::string& reason);
  void adapt_time_step(const int& rejections);
  int get_nonlinear_iterations() const;

  // Largest ratio between two successive time-steps
  double max_time_step_ratio_{std::numeric_limits<double>::infinity()};
  // Time at the beginning of the time-step, restored if the time-step is rejected
  double saved_time_{0.};

  double get_max_time_step_ratio() const;

 public:
  explicit TimeDiscretization(const Parameters& params, Args... couplings);
  explicit TimeDiscretization(const std::function<double(double)>& given_time_step,
                              const Parameters& params, Args... couplings);
  void solve();
  void get_tree();
  ~TimeDiscretization();
};

#include "Time/Time.tpp"
