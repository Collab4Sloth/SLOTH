/**
 * @file DSolverBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-08-08
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>
#include <string>

#include "Options/Options.hpp"
#include "Solvers/SolverBase.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class SolverUMFPACK : public SolverBase<mfem::UMFPackSolver, DirectSolverType> {
 public:
  SolverUMFPACK();
  std::shared_ptr<mfem::UMFPackSolver> create_solver(DirectSolverType SOLVER,
                                                     const Parameters& params) override;

  ~SolverUMFPACK();
};

/**
 * @brief Construct a new SolverUMFPACK::SolverUMFPACK object
 *
 */
SolverUMFPACK::SolverUMFPACK() {}

/**
 * @brief Create a direct solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::UMFPackSolver> SolverUMFPACK::create_solver(DirectSolverType SOLVER,
                                                                  const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", UMFPACK_DefaultConstant::print_level);
  auto ss = std::make_shared<mfem::UMFPackSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  return ss;
}

/**
 * @brief Destroy the SolverUMFPACK::SolverUMFPACK object
 *
 */
SolverUMFPACK::~SolverUMFPACK() {}
