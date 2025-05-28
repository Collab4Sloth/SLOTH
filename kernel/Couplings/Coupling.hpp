/**
 * @file Coupling.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class to define a  Coupling object
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Utils/UtilsForDebug.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <class... Args>
class Coupling {
 private:
  std::string name_{"Unnamed Coupling"};
  std::tuple<Args...> problems_;

 public:
  explicit Coupling(const std::string& name, Args... problems);
  std::string get_name();

  void get_tree();
  void initialize(const int& iter, const double& initial_time);
  std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>> execute(
      const int& iter, double& next_time, const double& current_time,
      const double& current_time_step);
  void post_execute(const int& iter, const double& current_time, const double& current_time_step);
  void update();
  void post_processing(const int& iter, const double& current_time,
                       const double& current_time_step);
  void finalize();
  ~Coupling();
};

/**
 * @brief Construct a new Coupling< Args...>:: Coupling object
 *
 * @tparam Args
 * @param problems
 */
// Coupling<Args...>::Coupling(const std::string& name, Args&&... problems)
template <class... Args>
Coupling<Args...>::Coupling(const std::string& name, Args... problems)
    : name_(name), problems_(std::make_tuple(std::forward<Args>(problems)...)) {}

/**
 * @brief Return the name of the coupling
 *
 * @tparam Args
 * @return const std::string
 */
template <class... Args>
std::string Coupling<Args...>::get_name() {
  return this->name_;
}

/**
 * @brief List all problems involved in the coupling
 *
 * @tparam Args
 */
template <class... Args>
void Coupling<Args...>::get_tree() {
  std::apply(
      [](auto&... problem) { (SlothInfo::verbose("   - Problem: ", problem.get_name()), ...); },
      this->problems_);
}

/**
 * @brief Initialize all the problems inside the coupling
 *
 * @tparam Args
 * @param initial_time
 */
template <class... Args>
void Coupling<Args...>::initialize(const int& iter, const double& initial_time) {
  std::apply(
      [iter, initial_time](auto&... problem) {
        (problem.initialize(initial_time), ...);
        (void([&problem, iter, initial_time] { problem.save(iter, initial_time); }()), ...);
      },
      problems_);
}

/**
 * @brief Solve all the problems inside the coupling
 *
 * @tparam Args
 * @param iter
 * @param current_time
 * @param current_time_step
 * @return std::vector<std::tuple<bool, double, mfem::Vector>>
 */
template <class... Args>
std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>> Coupling<Args...>::execute(
    const int& iter, double& next_time, const double& current_time,
    const double& current_time_step) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose(" ============================== ");
    SlothInfo::verbose(" ==== Coupling : ", this->name_);
    SlothInfo::verbose(" ============================== ");
  }

  std::vector<std::tuple<bool, double, std::vector<mfem::Vector>>> results;
  std::apply(
      [iter, &next_time, current_time, current_time_step, &results](auto&... problem) {
        double pp_next_time = current_time;
        ((pp_next_time = current_time,
          results.emplace_back(
              problem.execute(iter, pp_next_time, current_time, current_time_step)),
          next_time = pp_next_time, problem.update()),
         ...);
      },
      problems_);

  return results;
}

/**
 * @brief Call recursively the post_execute methods of the problem
 *
 * @tparam Args
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class... Args>
void Coupling<Args...>::post_execute(const int& iter, const double& current_time,
                                     const double& current_time_step) {
  std::apply(
      [iter, current_time, current_time_step](auto&... problem) {
        (problem.post_execute(iter, current_time, current_time_step), ...);
      },
      problems_);
}
/**
 * @brief Update the variables of the problems inside this coupling
 *
 * @tparam Args
 */
template <class... Args>
void Coupling<Args...>::update() {
  // std::apply([](auto&... problem) { (problem.update(), ...); }, problems_);
}

/**
 * @brief Save the variables of the problems inside this coupling
 *
 * @tparam Args
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class... Args>
void Coupling<Args...>::post_processing(const int& iter, const double& current_time,
                                        const double& current_time_step) {
  std::apply(
      [iter, current_time, current_time_step](auto&... problem) {
        (problem.post_processing(iter, current_time, current_time_step), ...);
      },
      problems_);
}

/**
 * @brief Call the finalize methods of each problem
 *
 * @tparam Args
 */
template <class... Args>
void Coupling<Args...>::finalize() {
  std::apply([](auto&... problem) { (problem.finalize(), ...); }, problems_);
}

/**
 * @brief Destroy the Coupling< Args...>:: Coupling object
 *
 * @tparam Args
 */
template <class... Args>
Coupling<Args...>::~Coupling() {}
