/**
 * @anchor ProblemBase
 * @file ProblemBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for building Problems
 * @version 0.1
 * @date 2025-09-05
 *
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
#include <functional>
#include <list>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/Coefficients.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Concept for VAR types
 *
 * @tparam T
 * @tparam VAR
 */
template <class T, class VAR>
concept PbVar = std::same_as<std::remove_cvref_t<T>, VAR>;

/**
 * @brief Base class for SLOTH problem
 *
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
class ProblemBase {
 private:
  void check_convergence(const std::vector<std::unique_ptr<mfem::Vector>>& unks);

 protected:
  std::string name_{"Unnamed problem"};
  VAR& variables_;
  std::vector<VAR*> auxvariables_;

  std::vector<Coefficients> coefficients_;

  PST& pst_;
  std::vector<mfem::Vector> unknown_;
  std::shared_ptr<Convergence> convergence_;
  std::vector<std::tuple<std::string, bool, double>> var_convergence_;

 public:
  template <PbVar<VAR>... Args>
  ProblemBase(const std::string& name, VAR& variables, const std::vector<Coefficients>& Coeff,
              PST& pst, Args&&... auxvariables);

  template <PbVar<VAR>... Args>
  ProblemBase(const std::string& name, VAR& variables, const std::vector<Coefficients>& Coeff,
              Convergence& convergence, PST& pst, Args&&... auxvariables);

  template <PbVar<VAR>... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst, Args&&... auxvariables);

  template <PbVar<VAR>... Args>
  ProblemBase(const std::string& name, VAR& variables, Convergence& convergence, PST& pst,
              Args&&... auxvariables);

  virtual ~ProblemBase() = default;

  std::string get_name();
  VAR get_problem_variables();

  std::vector<std::tuple<std::string, bool, double>> get_convergence();

  /////////////////////////////////////////////////////

  virtual void initialize([[maybe_unused]] const double& initial_time) {}

  void execute(const int& iter, double& next_time, const double& current_time,
               const double& current_time_step);

  virtual void do_time_step(double& next_time, const double& current_time, double current_time_step,
                            const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                            const std::vector<std::vector<std::string>>& unks_info) = 0;

  void update();

  virtual void post_processing(const int& iter, const double& current_time);

  void save(const int& iter, const double& current_time);

  virtual void finalize();
};

#include "Problems/ProblemBase.tpp"
