/**
 * @file MPI_Problem.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to defined a MPI_Problem objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */

#pragma once
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/ProblemBase.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class VAR, class PST>
class MPI_Problem : public ProblemBase<VAR, PST> {
 public:
  MPI_Problem(VAR& variables, PST& pst, const PhysicalConvergence& convergence);

  /////////////////////////////////////////////////////
  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::BlockVector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  /////////////////////////////////////////////////////

  ~MPI_Problem();
};

/**
 * @brief Construct a new Unsteady MPI_Problem object
 *
 * @tparam OPE
 * @tparam SOLVER
 * @param oper
 * @param solver
 * @param convergence
 */
template <class VAR, class PST>
MPI_Problem<VAR, PST>::MPI_Problem(VAR& variables, PST& pst, const PhysicalConvergence& convergence)
    : ProblemBase<VAR, PST>("MPI problem", variables, pst, convergence) {}

/**
 * @brief Do a time-step by calling Step method of the ODE
 *
 * @tparam OPE
 * @tparam VAR
 * @param unk
 * @param current_time
 * @param current_time_step
 */
template <class VAR, class PST>
void MPI_Problem<VAR, PST>::do_time_step(double& next_time, const double& current_time,
                                         double current_time_step, const int iter,
                                         std::vector<std::unique_ptr<mfem::BlockVector>>& vect_unk,
                                         const std::vector<std::vector<std::string>>& unks_info) {
  int rank = mfem::Mpi::WorldRank();

  auto& unk = *(vect_unk[0]);
  unk = static_cast<double>(rank);
  // Store the solution into a temporary mfem::Vector that will be used during updating stage,
  // TODO(cci) : utile? à voir avec le restart de paraview dan mfem et peut etre l'AMR
  this->unknown_.emplace_back(unk);
}

/**
 * @brief Destroy the Problem object
 *
 * @tparam OPE
 * @tparam SOLVER
 */
template <class VAR, class PST>
MPI_Problem<VAR, PST>::~MPI_Problem() {}
