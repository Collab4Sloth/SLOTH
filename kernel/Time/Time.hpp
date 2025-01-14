/**
 * @file Time.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Solving algorithm
 * @version 0.1
 * @date 2024-05-31
 *
 * Copyright CEA (c) 2024
 *
 */

#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Problems/Problem.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variable.hpp"
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

  std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>> vect_tup_pb_convergence_;
  void get_parameters();

  void check_data_before_execute();

  void initialize();
  std::vector<std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>>> execute(
      const int& iter);
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
std::vector<std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>>>
TimeDiscretization<Args...>::execute(const int& iter) {
  Catch_Time_Section("TimeDiscretization::execute");

  std::vector<std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>>> results;
  auto current_time = this->current_time_;
  auto next_time = this->current_time_;
  const auto& current_time_step = this->current_time_step_;

  std::apply(
      [iter, &next_time, current_time, current_time_step, &results](auto&... coupling) {
        double cc_next_time = current_time;
        ((cc_next_time = current_time,
          results.emplace_back(
              coupling.execute(iter, cc_next_time, current_time, current_time_step)),
          next_time = cc_next_time),
         ...);
      },
      couplings_);

  // TODO(cci): ici c'est le dernier current_time. A améliorer dans le cadre d'une gestion du pas de
  // temps
  this->current_time_ = next_time;

  return results;
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
    const auto& results = this->execute(iter);
    const auto& tt = this->current_time_;
    std::apply([iter, tt](auto&... coupling) { (coupling.initialize(iter, tt), ...); }, couplings_);
    //------------
    // Check convergence
    //------------
    // TODO(cci): à implémenter pour ne faire l'update qu'à ce moment là

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
 * @brief Destroy the Time Discretization:: Time Discretization object
 *
 */
template <class... Args>
TimeDiscretization<Args...>::~TimeDiscretization() {}
