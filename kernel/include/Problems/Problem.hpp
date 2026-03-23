/**
 * @file Problem.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to defined a PDE Problem
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
#include <list>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/ProblemBase.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Define a PDE Problem
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 *
 */
template <class OPE, class VAR, class PST>
class Problem : public ProblemBase<VAR, PST> {
 private:
  OPE oper_;
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};
  void save_specialized(bool must_be_saved);

 protected:
  void set_time_coefficients(double time) override;

 public:
  template <PbVar<VAR>... Args>
  Problem(const OPE& oper, VAR& variables, const std::vector<Coefficients>& Coeff,
          Convergence& convergence, PST& pst, Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Problem(const OPE& oper, VAR& variables, const std::vector<Coefficients>& Coeff, PST& pst,
          Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Problem(const std::string& name, const OPE& oper, VAR& variables,
          const std::vector<Coefficients>& Coeff, Convergence& convergence, PST& pst,
          Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Problem(const std::string& name, const OPE& oper, VAR& variables,
          const std::vector<Coefficients>& Coeff, PST& pst, Args&&... auxvariable);

  /////////////////////////////////////////////////////
  void initialize(const double& initial_time) override;

  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  void post_processing(const int& iter, const double& current_time) override;

  void finalize() override;
};

#include "Problems/Problem.tpp"
