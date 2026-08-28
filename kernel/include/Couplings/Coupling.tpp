/**
 * @file Coupling.tpp
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

#include "Utils/UtilsForData.hpp"
#include "Utils/UtilsForDebug.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Coupling< Args...>:: Coupling object
 *
 * @tparam Args Types of the problems stored in the coupling.
 * @param problems
 */
// Coupling<Args...>::Coupling(const std::string& name, Args&&... problems)
template <class... Args>
Coupling<Args...>::Coupling(const std::string& name, Args... problems)
    : name_(name), problems_(std::make_tuple(std::forward<Args>(problems)...)) {
  this->check_duplicate_problem_name();
}

/**
 * @brief Return the name of the coupling
 *
 * @tparam Args Types of the problems stored in the coupling.
 * @return const std::string
 */
template <class... Args>
std::string Coupling<Args...>::get_name() {
  return this->name_;
}

/**
 * @brief Set the name of the coupling
 *
 * @tparam Args Types of the problems stored in the coupling.
 * @param std::string The name of the coupling
 */
template <class... Args>
void Coupling<Args...>::set_name(const std::string& new_name) {
  this->name_ = new_name;
}

/**
 * @brief Checks for duplicate problem names and renames duplicates.
 *
 * @tparam Args Types of the problems stored in the coupling.
 */
template <class... Args>
void Coupling<Args...>::check_duplicate_problem_name() {
  std::set<std::string> problem_names;

  std::apply(
      [&](auto&... problem) {
        (
            [&] {
              const std::string original_name = problem.get_name();
              std::string new_name = original_name;

              int suffix = 0;
              while (problem_names.contains(new_name)) {
                new_name = original_name + std::to_string(++suffix);
              }

              if (new_name != original_name) {
                problem.set_name(new_name);
              }

              problem_names.insert(new_name);
            }(),
            ...);
      },
      problems_);
}

/**
 * @brief Checks for duplicate post-processing directories and add prefix to avoid duplicates in CSV
 * files.
 *
 * @tparam Args Types of the problems stored in the coupling.
 */
template <class... Args>
void Coupling<Args...>::check_duplicate_post_processing_directory(std::size_t nb_couplings) {
  std::map<std::string, int> directory_count;
  std::apply(
      [&](auto&... problem) {
        (
            [&] {
              const std::string pb_dir = problem.get_problem_post_processing_directory();

              if (!pb_dir.empty()) {
                ++directory_count[pb_dir];
              }
            }(),
            ...);
      },
      problems_);

  std::apply(
      [&](auto&... problem) {
        (
            [&] {
              const std::string pb_dir = problem.get_problem_post_processing_directory();
              if (!pb_dir.empty() && directory_count[pb_dir] > 1) {
                std::string name = replaceSpaces(problem.get_name()) + "_time_specialized.csv";
                if (nb_couplings > 1) {
                  name = replaceSpaces(this->get_name()) + "_" + name;
                }
                problem.set_name_save_specialized(name);
              }
            }(),
            ...);
      },
      problems_);
}

/**
 * @brief List all problems involved in the coupling
 *
 * @tparam Args Types of the problems stored in the coupling.
 */
template <class... Args>
void Coupling<Args...>::get_tree() {
  std::apply(
      [](auto&... problem) { (SlothInfo::verbose("   - Problem: ", problem.get_name()), ...); },
      this->problems_);
}

/**
 * @brief Initialize all the problems inside the coupling, then apply AMR
 *        pre-refinement and save the initial state.
 *
 * @details For each Problem, in order: `Problem::initialize()`, then
 *          `Problem::initialize_amr()` (a no-op unless that Problem has
 *          an AMR strategy attached), then `Problem::save()`. `all_vars`
 *          (every Problem's variables of every Coupling).
 *
 * @tparam Args Problem types in this coupling.
 * @param iter Iteration number used when saving the initial state.
 * @param initial_time Initial simulation time.
 * @param time_step Initial time step, forwarded to each Problem's
 *                  `initialize()`.
 */
template <class... Args>
void Coupling<Args...>::initialize(const int& iter, const double& initial_time,
                                   const double time_step, std::vector<VAR*> all_vars) {
  std::apply(
      [iter, initial_time, time_step, all_vars](auto&... problem) {
        (problem.initialize(initial_time, time_step), ...);
        (problem.initialize_amr(all_vars), ...);
        (void([&problem, iter, initial_time] { problem.save(iter, initial_time); }()), ...);
      },
      problems_);
}

/**
 * @brief Collect a reference to the Variables of every Problem in this
 *        coupling (called by TimeDiscretization<Args...>::CollectAllVariables())
 *
 * @tparam Args Problem types in this coupling.
 * @return std::vector<VAR*> One pointer per Problem, referencing its
 *        `variables_` member directly (not a copy) — valid for as long
 *        as this `Coupling` and its Problems are alive.
 */
template <class... Args>
std::vector<typename Coupling<Args...>::VAR*> Coupling<Args...>::CollectAllVariables() {
  std::vector<typename Coupling<Args...>::VAR*> all_vars;
  std::apply([&all_vars](
                 auto&... problem) { (all_vars.push_back(&problem.get_problem_variables()), ...); },
             this->problems_);
  return all_vars;
}

/**
 * @brief Solve all the problems inside the coupling
 *
 * @tparam Args Types of the problems stored in the coupling.
 * @param iter
 * @param next_time
 * @param current_time
 * @param current_time_step
 */
template <class... Args>
void Coupling<Args...>::execute(const int& iter, double& next_time, const double& current_time,
                                const double& current_time_step) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose(" ============================== ");
    SlothInfo::verbose(" ==== Coupling : ", this->name_);
    SlothInfo::verbose(" ============================== ");
  }

  this->pb_convergence_.clear();
  std::apply(
      [this, iter, &next_time, current_time, current_time_step](auto&... problem) {
        double pp_next_time = current_time;
        (
            [&] {
              pp_next_time = current_time;
              problem.execute(iter, pp_next_time, current_time, current_time_step);
              auto pb_conv = problem.get_convergence();
              if (!pb_conv.empty()) {
                this->pb_convergence_.emplace_back(
                    std::make_tuple(problem.get_name(), std::move(pb_conv)));
              }
              next_time = pp_next_time;
              problem.update();
            }(),
            ...);
      },
      problems_);
}

/**
 * @brief Update the variables of the problems inside this coupling
 *
 * @tparam Args Types of the problems stored in the coupling.
 */
template <class... Args>
void Coupling<Args...>::update() {
  // std::apply([](auto&... problem) { (problem.update(), ...); }, problems_);
}

/**
 * @brief Save the variables of the problems inside this coupling, then
 *        apply one AMR derefine/refine pass.
 *
 * @details For each Problem: `Problem::post_processing()` (saves its
 *          variables), then `Problem::save_amr()` (a no-op unless that
 *          Problem has an AMR strategy attached). `all_vars` (every
 *          Problem's variables of every Coupling).
 *
 * @tparam Args Problem types in this coupling.
 * @param iter Current iteration number.
 * @param current_time Current simulation time.
 */
template <class... Args>
void Coupling<Args...>::post_processing(const int& iter, const double& current_time,
                                        std::vector<VAR*> all_vars) {
  std::apply(
      [iter, current_time, all_vars](auto&... problem) {
        (problem.post_processing(iter, current_time), ...);
        (problem.save_amr(all_vars), ...);
      },
      problems_);
}

/**
 * @brief Call the finalize methods of each problem
 *
 * @tparam Args Types of the problems stored in the coupling.
 */
template <class... Args>
void Coupling<Args...>::finalize() {
  std::apply([](auto&... problem) { (problem.finalize(), ...); }, problems_);
}

/**
 * @brief
 *
 * @tparam Args Types of the problems stored in the coupling.
 * @return std::vector<std::tuple<std::string, std::vector<std::tuple<std::string, bool, double>>>>
 */
template <class... Args>
std::vector<std::tuple<std::string, std::vector<std::tuple<std::string, bool, double>>>>
Coupling<Args...>::get_convergence() {
  return this->pb_convergence_;
}

/**
 * @brief Destroy the Coupling< Args...>:: Coupling object
 *
 * @tparam Args Types of the problems stored in the coupling.
 */
template <class... Args>
Coupling<Args...>::~Coupling() {}
