/**
 * @file Problem.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to defined a Problem objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
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

#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/ProblemBase.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class OPE, class VAR, class PST>
class Problem : public ProblemBase<VAR, PST> {
 private:
  OPE oper_;
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};

 public:
  template <class... Args>
  Problem(const OPE& oper, VAR& variables, PST& pst, const PhysicalConvergence& convergence,
          std::list<int> pop_elem, Args&&... auxvariable);

  template <class... Args>
  Problem(const OPE& oper, VAR& variables, PST& pst, const PhysicalConvergence& convergence,
          Args&&... auxvariable);

  template <class... Args>
  Problem(const std::string& name, const OPE& oper, VAR& variables, PST& pst,
          const PhysicalConvergence& convergence, std::list<int> pop_elem, Args&&... auxvariables);

  template <class... Args>
  Problem(const std::string& name, const OPE& oper, VAR& variables, PST& pst,
          const PhysicalConvergence& convergence, Args&&... auxvariable);

  /////////////////////////////////////////////////////
  void initialize(const double& initial_time) override;

  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  void post_execute(const int& iter, const double& current_time,
                    const double& current_time_step) override;
  /////////////////////////////////////////////////////

  void post_processing(const int& iter, const double& current_time,
                       const double& current_time_step) override;

  void finalize() override;

  /////////////////////////////////////////////////////

  ~Problem();
};

/**
 * @brief Construct a new Problem< O P E,  V A R,  P S T>:: Problem object
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param oper
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class OPE, class VAR, class PST>
template <class... Args>
Problem<OPE, VAR, PST>::Problem(const OPE& oper, VAR& variables, PST& pst,
                                const PhysicalConvergence& convergence, std::list<int> pop_elem,
                                Args&&... auxvariables)
    : ProblemBase<VAR, PST>("AllenCahn problem", variables, pst, convergence, pop_elem,
                            auxvariables...),
      oper_(oper) {}

/**
 * @brief Construct a new Problem< O P E,  V A R,  P S T>:: Problem object
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param oper
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class OPE, class VAR, class PST>
template <class... Args>
Problem<OPE, VAR, PST>::Problem(const OPE& oper, VAR& variables, PST& pst,
                                const PhysicalConvergence& convergence, Args&&... auxvariables)
    : ProblemBase<VAR, PST>("AllenCahn problem", variables, pst, convergence, auxvariables...),
      oper_(oper) {}

/**
 * @brief Construct a new Problem< O P E,  V A R,  P S T>:: Problem object
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param oper
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class OPE, class VAR, class PST>
template <class... Args>
Problem<OPE, VAR, PST>::Problem(const std::string& name, const OPE& oper, VAR& variables, PST& pst,
                                const PhysicalConvergence& convergence, std::list<int> pop_elem,
                                Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, pop_elem, auxvariables...),
      oper_(oper) {}

/**
 * @brief Construct a new Problem< O P E,  V A R,  P S T>:: Problem object
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param oper
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class OPE, class VAR, class PST>
template <class... Args>
Problem<OPE, VAR, PST>::Problem(const std::string& name, const OPE& oper, VAR& variables, PST& pst,
                                const PhysicalConvergence& convergence, Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, auxvariables...), oper_(oper) {}

/**
 * @brief Initialize the calculation : operator + ODE
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::initialize(const double& initial_time) {
  this->oper_.initialize(initial_time, this->variables_, this->auxvariables_);
}

/**
 * @brief Compute Error and Energies associated with the problem if required
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::post_execute(const int& iter, const double& current_time,
                                          const double& current_time_step) {}

/**
 * @brief Show calculation information
 *
 * @tparam OPE
 * @tparam VAR
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::finalize() {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    if (!this->pst_.get_enable_save_specialized_at_iter()) {
      this->pst_.save_specialized(this->oper_.get_time_specialized());
    }
    if (!this->pst_.get_iso_val_to_compute().empty()) {
      for (const auto& [var_name, variable_time_iso_specialized] :
           this->oper_.get_time_iso_specialized()) {
        std::string str = var_name + "_iso_computation.csv";
        this->pst_.save_iso_specialized(variable_time_iso_specialized, str);
      }
    }

    SlothInfo::verbose(" ");
    SlothInfo::verbose(" ============================== ");
    SlothInfo::verbose(" Results are saved in the folder : ",
                       this->pst_.get_post_processing_directory());
  }
}

/**
 * @brief  Call the post_execute method of the given problem and saves variables according with
 * PST
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::post_processing(const int& iter, const double& current_time,
                                             const double& current_time_step) {
  ////
  const auto nvars = this->variables_.get_variables_number();
  std::vector<mfem::Vector> u_vect;
  const std::map<std::string, double> map_iso_value = this->pst_.get_iso_val_to_compute();
  const std::map<std::string, double> map_integral = this->pst_.get_integral_to_compute();
  // Errors
  for (auto iv = 0; iv < nvars; iv++) {
    auto vv = this->variables_.getIVariable(iv);
    auto solution = vv.get_analytical_solution();
    auto unk = vv.get_unknown();
    auto unk_name = vv.getVariableName();
    // Isovalues
    if (!map_iso_value.empty() && map_iso_value.contains(unk_name)) {
      const double iso_value = map_iso_value.at(unk_name);
      this->oper_.ComputeIsoVal(iter, current_time, current_time_step, iv, unk_name, unk,
                                iso_value);
    }
    // Errors
    if (solution != nullptr) {
      auto solution_func = solution.get();
      this->oper_.ComputeError(iter, current_time, current_time_step, iv, unk_name, unk,
                               *solution_func);
    }
    // Integral
    if (!map_integral.empty() && map_integral.contains(unk_name)) {
      const double integral_threshold = map_integral.at(unk_name);
      this->oper_.ComputeIntegral(iter, current_time, current_time_step, iv, unk_name, unk,
                                  integral_threshold);
    }
    u_vect.emplace_back(std::move(unk));
  }
  // Energies
  this->oper_.ComputeEnergies(iter, current_time, current_time_step, u_vect);

  // Save for visualization
  ProblemBase<VAR, PST>::post_processing(iter, current_time, current_time_step);
  if (this->pst_.get_enable_save_specialized_at_iter()) {
    this->pst_.save_specialized(this->oper_.get_time_specialized());
  }
}

/**
 * @brief Do a time-step by calling Step method of the ODE
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param next_time
 * @param current_time
 * @param current_time_step
 * @param iter
 * @param vect_unk
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::do_time_step(double& next_time, const double& current_time,
                                          double current_time_step, const int iter,
                                          std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
                                          const std::vector<std::vector<std::string>>& unks_info) {
  // auto& unk = *(vect_unk[0]);

  const size_t unk_size = vect_unk.size();

  // this->oper_.solve(unk, next_time, current_time, current_time_step, iter);
  this->oper_.solve(vect_unk, next_time, current_time, current_time_step, iter);

  // Store the solution into a temporary mfem::Vector that will be used during updating stage, if
  // calculation converges
  // this->unknown_.emplace_back(unk);

  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    this->unknown_.emplace_back(unk_i);
  }
}

/**
 * @brief Destroy the Problem object
 *
 * @tparam OPE
 * @tparam SOLVER
 */
template <class OPE, class VAR, class PST>
Problem<OPE, VAR, PST>::~Problem() {}
