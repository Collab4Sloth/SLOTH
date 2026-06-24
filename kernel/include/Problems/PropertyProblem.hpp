/**
 * @file PropertyProblem.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Class used to compute primary variables as a function of auxialiary variables
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
#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "Problems/ProblemBase.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class PROPERTY, class VAR, class PST>
class Property_problem : public ProblemBase<VAR, PST> {
 private:
  // Property object manage by the Property_problem
  PROPERTY* PP_;

  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
  get_output_system(const std::vector<std::vector<std::string>>& unks_info,
                    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk);

  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> get_input_system();

 public:
  template <PbVar<VAR>... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables,
                   const std::vector<Coefficients>& Coeff, Convergence& convergence, PST& pst,
                   Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const Parameters& params, VAR& variables, const std::vector<Coefficients>& Coeff,
                   Convergence& convergence, PST& pst, Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables,
                   const std::vector<Coefficients>& Coeff, PST& pst, Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const Parameters& params, VAR& variables, const std::vector<Coefficients>& Coeff,
                   PST& pst, Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables,
                   Convergence& convergence, PST& pst, Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const Parameters& params, VAR& variables, Convergence& convergence, PST& pst,
                   Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                   Args&&... auxvariable);

  template <PbVar<VAR>... Args>
  Property_problem(const Parameters& params, VAR& variables, PST& pst, Args&&... auxvariable);

  /////////////////////////////////////////////////////
  void initialize(const double& initial_time, const double time_step) override;

  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  ~Property_problem();
};

#include "Problems/PropertyProblem.tpp"
