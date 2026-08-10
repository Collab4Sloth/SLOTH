/**
 * @file Problem.tpp
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

#include "AMR/ListAMR.hpp"
#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/ProblemBase.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a PDE Problem object.
 *
 * Initializes a PDE problem by binding the operator, primary variables,
 * auxiliary variables, coefficients, convergence criteria, and post-processing. The operator is
 * configured with the provided coefficients and energy-computation settings from the problem state.
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 * @tparam Args
 *   Types of auxiliary variables satisfying the `PbVar<VAR>` constraint.
 *
 * @param oper
 *   Nonlinear operator associated with the problem.
 * @param variables
 *   Primary problem variables.
 * @param Coeff
 *   Collection of coefficients used by the operator.
 * @param convergence
 *   Convergence criteria.
 * @param pst
 *   PostProcessing.
 * @param auxvariables
 *   Optional auxiliary variables forwarded to the base problem.
 *
 * @pre
 * - `oper` must be fully constructed.
 * - `Coeff` must be compatible with the operator.
 *
 */
template <class OPE, class VAR, class PST>
template <PbVar<VAR>... Args>
Problem<OPE, VAR, PST>::Problem(const OPE& oper, VAR& variables,
                                const std::vector<Coefficients>& Coeff, Convergence& convergence,
                                PST& pst, Args&&... auxvariables)
    : ProblemBase<VAR, PST>("Default NonLinear problem", variables, Coeff, convergence, pst,
                            auxvariables...),
      oper_(oper) {
  MFEM_VERIFY(variables.get_variables_number() == Coeff.size(),
              "Error: the number of Coefficients objects must be equal to the number of primary "
              "variables. Please check your data.");
  this->oper_.set_coefficients(Coeff, this->pst_.get_enable_compute_energies());
}

/**
 * @brief Construct a PDE Problem object.
 *
 * Initializes a PDE problem by binding the operator, primary variables,
 * auxiliary variables, coefficients,  and post-processing. The operator is
 * configured with the provided coefficients and energy-computation settings from the problem state.
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 * @tparam Args
 *   Types of auxiliary variables satisfying the `PbVar<VAR>` constraint.
 *
 * @param oper
 *   Nonlinear operator associated with the problem.
 * @param variables
 *   Primary problem variables.
 * @param Coeff
 *   Collection of coefficients used by the operator.
 * @param pst
 *   PostProcessing.
 * @param auxvariables
 *   Optional auxiliary variables forwarded to the base problem.
 *
 * @pre
 * - `oper` must be fully constructed.
 * - `Coeff` must be compatible with the operator.
 *
 */
template <class OPE, class VAR, class PST>
template <PbVar<VAR>... Args>
Problem<OPE, VAR, PST>::Problem(const OPE& oper, VAR& variables,
                                const std::vector<Coefficients>& Coeff, PST& pst,
                                Args&&... auxvariables)
    : ProblemBase<VAR, PST>("Default NonLinear problem", variables, Coeff, pst, auxvariables...),
      oper_(oper) {
  MFEM_VERIFY(variables.get_variables_number() == Coeff.size(),
              "Error: the number of Coefficients objects must be equal to the number of primary "
              "variables. Please check your data.");
  this->oper_.set_coefficients(Coeff, this->pst_.get_enable_compute_energies());
}

/**
 * @brief Construct a PDE Problem object.
 *
 * Initializes a PDE problem by binding the name, operator, primary variables,
 * auxiliary variables, coefficients,  convergence criteria, and post-processing. The operator is
 * configured with the provided coefficients and energy-computation settings from the problem state.
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 * @tparam Args
 *   Types of auxiliary variables satisfying the `PbVar<VAR>` constraint.
 *
 * @param name
 *   Name of the problem.
 * @param oper
 *   Nonlinear operator associated with the problem.
 * @param variables
 *   Primary problem variables.
 * @param Coeff
 *   Collection of coefficients used by the operator.
 * @param convergence
 *   Convergence criteria.
 * @param pst
 *   PostProcessing.
 * @param auxvariables
 *   Optional auxiliary variables forwarded to the base problem.
 *
 * @pre
 * - `oper` must be fully constructed.
 * - `Coeff` must be compatible with the operator.
 *
 */
template <class OPE, class VAR, class PST>
template <PbVar<VAR>... Args>
Problem<OPE, VAR, PST>::Problem(const std::string& name, const OPE& oper, VAR& variables,
                                const std::vector<Coefficients>& Coeff, Convergence& convergence,
                                PST& pst, Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, Coeff, convergence, pst, auxvariables...),
      oper_(oper) {
  MFEM_VERIFY(variables.get_variables_number() == Coeff.size(),
              "Error: the number of Coefficients objects must be equal to the number of primary "
              "variables. Please check your data.");
  this->oper_.set_coefficients(Coeff, this->pst_.get_enable_compute_energies());
}

/**
 * @brief Construct a PDE Problem object.
 *
 * Initializes a PDE problem by binding the name, operator, primary variables,
 * auxiliary variables, coefficients,  and post-processing. The operator is
 * configured with the provided coefficients and energy-computation settings from the problem state.
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 * @tparam Args
 *   Types of auxiliary variables satisfying the `PbVar<VAR>` constraint.
 *
 * @param name
 *   Name of the problem.
 * @param oper
 *   Nonlinear operator associated with the problem.
 * @param variables
 *   Primary problem variables.
 * @param Coeff
 *   Collection of coefficients used by the operator.
 * @param pst
 *   PostProcessing.
 * @param auxvariables
 *   Optional auxiliary variables forwarded to the base problem.
 *
 * @pre
 * - `oper` must be fully constructed.
 * - `Coeff` must be compatible with the operator.
 *
 */
template <class OPE, class VAR, class PST>
template <PbVar<VAR>... Args>
Problem<OPE, VAR, PST>::Problem(const std::string& name, const OPE& oper, VAR& variables,
                                const std::vector<Coefficients>& Coeff, PST& pst,
                                Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, Coeff, pst, auxvariables...), oper_(oper) {
  MFEM_VERIFY(variables.get_variables_number() == Coeff.size(),
              "Error: the number of Coefficients objects must be equal to the number of primary "
              "variables. Please check your data.");
  this->oper_.set_coefficients(Coeff, this->pst_.get_enable_compute_energies());
}

/**
 * @brief Initialize the Problem
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 * @param initial_time
 *   Initial time
 * @param time_step
 *   current time_step
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::initialize(const double& initial_time, const double time_step) {
  this->oper_.initialize(initial_time, time_step, this->variables_, this->auxvariables_);
  if (this->amr_ != nullptr) {
    this->amr_->InitialRefine(this->variables_, this->auxvariables_);
  }
}
/**
 * @brief
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::save(const int& iter, const double& current_time) {
  ProblemBase<VAR, PST>::save(iter, current_time);

  if (this->amr_ != nullptr && iter > 0) {
    this->amr_->StepDerefine(this->variables_, this->auxvariables_);
    this->amr_->StepRefine(this->variables_, this->auxvariables_);
  }
}

/**
 * @brief Finalize the Problem
 *
 * Save specialized results in CSV files
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::finalize() {
  // Save specialized values only at the end of calculation if required
  bool must_be_saved = !this->pst_.get_enable_save_specialized_at_iter();
  this->save_specialized(must_be_saved);
}

/**
 * @brief Save specialized results
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 *
 * @param must_be_saved Flag used to indicate if specialized results must be saved
 *
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::save_specialized(bool must_be_saved) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    if (must_be_saved) {
      if (!this->oper_.get_time_specialized().empty()) {
        this->pst_.save_specialized(this->oper_.get_time_specialized());
      }

      if (!this->oper_.get_time_iso_specialized().empty()) {
        for (const auto& [var_name, variable_time_iso_specialized] :
             this->oper_.get_time_iso_specialized()) {
          std::string str = var_name + "_iso_computation.csv";
          this->pst_.save_specialized(variable_time_iso_specialized, str);
        }
      }
    }
  }
  if (must_be_saved) {
    this->oper_.clear_time_specialized();
    this->oper_.clear_iso_time_specialized();
  }
}
/**
 * @brief Saves variables according PostProcessing settings
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 *
 * @param iter Current iteration
 * @param current_time Current time
 *
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::post_processing(const int& iter, const double& current_time) {
  const auto nvars = this->variables_.get_variables_number();
  std::vector<mfem::Vector> u_vect;
  const double current_time_step = this->oper_.get_current_time_step();

  // Isovalue to compute by variable
  const std::map<std::string, double> map_iso_value = this->pst_.get_iso_val_to_compute();

  // Integral value to compute by variable
  const std::map<std::string, std::tuple<double, double>> map_integral =
      this->pst_.get_integral_to_compute();

  for (unsigned int iv = 0; iv < nvars; iv++) {
    auto vv = this->variables_[iv];
    auto solution = vv.get_analytical_solution();
    auto unk = vv.get_unknown();
    auto unk_name = vv.getVariableName();

    // Errors
    if (solution != nullptr) {
      auto solution_func = solution.get();
      this->oper_.ComputeError(iter, current_time, current_time_step, iv, unk_name, unk,
                               *solution_func);
    }

    // Isovalues
    if (!map_iso_value.empty() && map_iso_value.contains(unk_name)) {
      const double iso_value = map_iso_value.at(unk_name);
      this->oper_.ComputeIsoVal(iter, current_time, current_time_step, iv, unk_name, unk,
                                iso_value);
    }

    // Integral
    if (!map_integral.empty() && map_integral.contains(unk_name)) {
      const auto& [lower_bound, upper_bound] = map_integral.at(unk_name);
      this->oper_.ComputeIntegral(iter, current_time, current_time_step, iv, unk_name, unk,
                                  lower_bound, upper_bound);
    }
    u_vect.emplace_back(std::move(unk));
  }
  // Compute Energies is required (True by default)
  if (this->pst_.get_enable_compute_energies()) {
    this->oper_.ComputeEnergies(iter, current_time, current_time_step, u_vect);
  }
  // Save variables for visualization
  ProblemBase<VAR, PST>::post_processing(iter, current_time);

  // Save specialized values at each time-step if required
  bool must_be_saved = this->pst_.get_enable_save_specialized_at_iter();
  this->save_specialized(must_be_saved);
}

/**
 * @brief Do a time-step by calling Step method of the ODE
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *   PostProcessing.
 *
 * @param next_time Next time.
 * @param current_time Current time.
 * @param current_time_step Current time-step.
 * @param iter Current iteration.
 * @param vect_unk Vector of pointers towards the unknowns.
 * @param unks_info Vector containing information associated with the unknowns.
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::do_time_step(
    double& next_time, const double& current_time, double current_time_step, const int iter,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
    [[maybe_unused]] const std::vector<std::vector<std::string>>& unks_info) {
  const size_t unk_size = vect_unk.size();

  // Set time for coefficients, /!\ this is correct only for implicit time scheme
  // TODO: use effective time step according to the time scheme
  this->set_time_coefficients(next_time+current_time_step);

  this->oper_.setGeometry(this->geometry_);
  this->oper_.solve(vect_unk, next_time, current_time, current_time_step, iter);

  // Store the solution into a temporary mfem::Vector that will be used during updating stage, if
  // calculation converges
  for (size_t i = 0; i < unk_size; i++) {
    this->unknown_.emplace_back(*(vect_unk[i]));
  }
}

/**
 * @brief Call OperatorBase method to set the time to coefficients
 *
 * @tparam OPE
 *   Operator type defining the nonlinear PDE.
 * @tparam VAR
 *   Primary variable container type.
 * @tparam PST
 *
 * @param time current time
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::set_time_coefficients(double time) {
  this->oper_.set_time_coefficients(time);
}
