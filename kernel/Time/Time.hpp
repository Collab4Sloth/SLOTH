/**
 * @file Time.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Main time loop 
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
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

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Profiling/Profiling.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <class... Args>
class TimeDiscretization {
 private:
  const Parameters& params_;
  std::tuple<Args...> couplings_;
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
  void post_execute(const int& iter);
  void update();
  void post_processing(const int& iter);
  void finalize();
  void time_management();
  void time_info(const int& iter);

 public:
  // explicit TimeDiscretization(const Parameters& params, Args&&... couplings);
  explicit TimeDiscretization(const Parameters& params, Args... couplings);
  void solve();
  void get_tree();
  ~TimeDiscretization();
};

/**
 * @brief Construct a new Time Discretization< Args...>:: Time Discretization object
 *
 * @tparam Args
 * @param params
 * @param couplings
 */
template <class... Args>
TimeDiscretization<Args...>::TimeDiscretization(const Parameters& params, Args... couplings)
    : params_(params), couplings_(std::make_tuple(std::forward<Args>(couplings)...)) {
  this->get_parameters();
}

/**
 * @brief Set all the global parameters
 *
 * @tparam Args
 * @param params
 */
template <class... Args>
void TimeDiscretization<Args...>::get_parameters() {
  this->initial_time_ =
      this->params_.template get_param_value_or_default<double>("initial_time", 0.);
  this->final_time_ = this->params_.template get_param_value<double>("final_time");
  this->time_step_ = this->params_.template get_param_value<double>("time_step");
  this->current_time_step_ = this->initial_time_;
}

/**
 * @brief Initialization stage
 *
 */
template <class... Args>
void TimeDiscretization<Args...>::initialize() {
  Catch_Time_Section("TimeDiscretization::initialize");

  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose(R"(@@@@@@@@@@@@@@@@%#+--=*%@@@@@@@@@@@@@@@@
@@@@@@@@@@@@%#+-:::--:::-+*%@@@@@@@@@@@@
@@@@@@@@%*+-:::-+#%@@%#*=:::-=*#@@@@@@@@
@@@@@%+-:::-+#%@@@@@@@@@@@#*=::::=@@@@@@
@@@@@=::+#%@@@@@@%%%%@@@@@@@@@%*::+@@@@@
@@@@*::#@@@@@#+--::::-=+*%@@@@@@*::#@@@@
@@@%::+@@@@*-:::::::::::::-+%@@@@=:-%@@@
@@@=:-@@@%=::::::::::::::--::=#@@%-:=@@@
@@*::#@@@-:::::::--:::::+*==:::-*@#::*@@
@#::+@@@#+#=::::+##+::::-%-+:::::-%+::%@
%-:=@@@@#==-%@#::::::::::--*:::::::*=:-@
*::#@@@@@+:-*#*===::::::::#*:::::::-#::*
@-:-@@@@@@#=---:::::::::=#@*:::::::#=:-@
@%::+@@@@@@@*++===+++*#%@@@#::::::**::%@
@@*::#@@@@#-::--:::--=+*%@@@=::::=#::*@@
@@@+::%@@%:::::::::::::::=#@%::::%-:=@@@
@@@@=:=@@@+::::::::::::::::=##::*+:-%@@@
@@@@%::=*%@%+-::::::::--==---##++::#@@@@
@@@@@#=::::=*##*+++*#%@@@@%*=-:::=*@@@@@
@@@@@@@@%*=-:::=+#%@@@#*=:::-=*#@@@@@@@@
@@@@@@@@@@@@%#+-:::-=:::-+*%@@@@@@@@@@@@
@@@@@@@@@@@@@@@@%#+--+#%@@@@@@@@@@@@@@@@
  )");
  }

  const auto& tt = this->initial_time_;
  this->current_time_ = tt;
  const auto& iter = 0;
  std::apply([iter, tt](auto&... coupling) { (coupling.initialize(iter, tt), ...); }, couplings_);
}

/**
 * @brief Call the post_processing method of each coupling
 *
 * @tparam Args
 * @param iter
 */
template <class... Args>
void TimeDiscretization<Args...>::post_processing(const int& iter) {
  Catch_Time_Section("TimeDiscretization::post_processing");

  const auto& current_time = this->current_time_;
  const auto& current_time_step = this->current_time_step_;

  std::apply(
      [iter, current_time, current_time_step](auto&... coupling) {
        (coupling.post_processing(iter, current_time, current_time_step), ...);
      },
      couplings_);
}

/**
 * @brief Call the finalize method of each coupling
 *
 * @tparam Args
 */
template <class... Args>
void TimeDiscretization<Args...>::finalize() {
  Catch_Time_Section("TimeDiscretization::finalize");

  std::apply([](auto&... coupling) { (coupling.finalize(), ...); }, couplings_);
}

/**
 * @brief Executtion of the solving algorithm
 *
 * @tparam Args
 * @param iter
 * @return std::vector<std::vector<std::tuple<bool, double, mfem::Vector>>>
 */
template <class... Args>
void TimeDiscretization<Args...>::execute(const int& iter) {
  Catch_Time_Section("TimeDiscretization::execute");

  auto current_time = this->current_time_;
  auto next_time = this->current_time_;
  const auto& current_time_step = this->current_time_step_;

  this->convergence_.clear();

  std::apply(
      [this, iter, &next_time, current_time, current_time_step](auto&... coupling) {
        double cc_next_time = current_time;
        (
            [&] {
              cc_next_time = current_time;
              coupling.execute(iter, cc_next_time, current_time, current_time_step);
              auto coupling_conv = coupling.get_convergence();
              if (!coupling_conv.empty()) {
                this->convergence_.emplace_back(coupling.get_name(), std::move(coupling_conv));
              }
              next_time = cc_next_time;
            }(),
            ...);
      },
      couplings_);

  // TODO(cci): ici c'est le dernier current_time. A améliorer dans le cadre d'une gestion du pas de
  // temps
  this->current_time_ = next_time;
}

/**
 * @brief Call the post_execute method of each coupling
 *
 * @tparam Args
 * @param iter
 */
template <class... Args>
void TimeDiscretization<Args...>::post_execute(const int& iter) {
  Catch_Time_Section("TimeDiscretization::post_execute");

  const auto& current_time = this->current_time_;
  const auto& current_time_step = this->current_time_step_;
  std::apply(
      [iter, current_time, current_time_step](auto&... coupling) {
        (coupling.post_execute(iter, current_time, current_time_step), ...);
      },
      couplings_);
}

/**
 * @brief Call the update method of each coupling
 *
 * @tparam Args
 */
template <class... Args>
void TimeDiscretization<Args...>::update() {
  Catch_Time_Section("TimeDiscretization::update");

  std::apply([](auto&... coupling) { (coupling.update(), ...); }, couplings_);
}

/**
 * @brief Compute the current value of the time-step
 *
 * @return double
 */
template <class... Args>
void TimeDiscretization<Args...>::time_management() {
  auto dt_real = this->time_step_;
  if ((this->current_time_ + this->time_step_) + 0.5 * this->time_step_ > this->final_time_) {
    this->last_step_ = true;
    dt_real = this->final_time_ - this->current_time_;
  }
  this->current_time_step_ = dt_real;
}

/**
 * @brief Print time-step information
 *
 * @tparam PST
 * @tparam OPE
 * @param iter
 * @param dt_requested
 * @param dt_real
 */
template <class... Args>
void TimeDiscretization<Args...>::time_info(const int& iter) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose("============================== ");
    SlothInfo::verbose("==== Iteration : ", iter);
    SlothInfo::verbose("     - time : ", this->current_time_);
    SlothInfo::verbose("     - time-step requested : ", this->time_step_);
    SlothInfo::verbose("     - time-step accepted  : ", this->current_time_step_);
  }
}

/**
 * @brief Run the calculation
 *
 */
template <class... Args>
void TimeDiscretization<Args...>::solve() {
  Catch_Time_Section("TimeDiscretization::solve");

  //=============================================
  //            INITIALIZATION
  //=============================================
  this->initialize();

  //=============================================
  //           TIME MARCHING
  //=============================================
  for (auto iter = 1; !this->last_step_; iter++) {
    //------------
    // Time-step
    //------------
    this->time_management();
    this->time_info(iter);

    //------------
    // Solve
    //------------
    this->execute(iter);

    //------------
    //  Post Execute
    //------------
    this->post_execute(iter);

    //------------
    //  Update
    //------------
    this->update();

    //-------------------
    // Save current
    //-------------------
    this->post_processing(iter);

    //------------
    // Check convergence
    //------------
    if (checkTimeStepConvergence(this->convergence_)) break;
  }

  //-------------------
  // Save global data
  // Build csv...
  //-------------------
  this->finalize();
  //=============================================
  //          end of loop over time-step
  //=============================================
}

/**
 * @brief Return the numerical coupling tree (recursively over the couplings and the problems)
 *
 * @tparam Args
 */
template <class... Args>
void TimeDiscretization<Args...>::get_tree() {
  // TODO(cci): pas terrible, il faudrait un get_info, get_dod...

  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose(" ============================== ");
    SlothInfo::verbose(" ====  Solving algorithm   ==== ");
    SlothInfo::verbose(" ============================== ");
    SlothInfo::verbose(" >> Time discretization: ");
    SlothInfo::verbose("    - Initial time: ", this->initial_time_);
    SlothInfo::verbose("    - Final time: ", this->final_time_);
    SlothInfo::verbose("    - Maximum time-step: ", this->time_step_);
    std::apply(
        [](auto&... coupling) {
          (SlothInfo::verbose(" >> Coupling: ", coupling.get_name()), ...);
          (coupling.get_tree(), ...);
        },
        couplings_);
  }
}

/**
 * @brief Check if all problems of all couplings have converged for each variable
 *
 * @tparam Args
 * @param convergence
 * @return true
 * @return false
 */
template <class... Args>
bool TimeDiscretization<Args...>::checkTimeStepConvergence(const auto& convergence) {
  if (convergence.empty()) return false;
  for (const auto& [coup_name, coup_vect] : convergence) {
    if (coup_vect.empty()) return false;

    for (const auto& [pb_name, pb_vect] : coup_vect) {
      if (pb_vect.empty()) return false;

      for (const auto& [var_name, has_cvg, criterion] : pb_vect) {
        if (!has_cvg) return false;
      }
    }
  }
  return true;
}

/**
 * @brief Destroy the Time Discretization:: Time Discretization object
 *
 */
template <class... Args>
TimeDiscretization<Args...>::~TimeDiscretization() {}
