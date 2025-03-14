
/**
 * @file HSolverBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Base class for Hypre solvers
 * @version 0.1
 * @date 2024-08-06
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

class SolverHyprePCG : public SolverBase<mfem::HyprePCG, HypreSolverType> {
 public:
  SolverHyprePCG() = default;
  std::shared_ptr<mfem::HyprePCG> create_solver(HypreSolverType SOLVER,
                                                const Parameters& params) override;

  ~SolverHyprePCG() = default;
};

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
DEBILE_INLINE std::shared_ptr<mfem::HyprePCG> SolverHyprePCG::create_solver(HypreSolverType SOLVER,
                                                              const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  auto ss = std::make_shared<mfem::HyprePCG>(MPI_COMM_WORLD);
  return ss;
}


/////////////////////////////

class SolverHypreGMRES : public SolverBase<mfem::HypreGMRES, HypreSolverType> {
 public:
  SolverHypreGMRES() = default;
  std::shared_ptr<mfem::HypreGMRES> create_solver(HypreSolverType SOLVER,
                                                  const Parameters& params) override;

  ~SolverHypreGMRES() = default;
};

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
DEBILE_INLINE std::shared_ptr<mfem::HypreGMRES> SolverHypreGMRES::create_solver(HypreSolverType SOLVER,
                                                                  const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());
  auto ss = std::make_shared<mfem::HypreGMRES>(MPI_COMM_WORLD);
  return ss;
}

/////////////////////////////

class SolverHypreFGMRES : public SolverBase<mfem::HypreFGMRES, HypreSolverType> {
 public:
  SolverHypreFGMRES() = default;
  std::shared_ptr<mfem::HypreFGMRES> create_solver(HypreSolverType SOLVER,
                                                   const Parameters& params) override;

  ~SolverHypreFGMRES() = default;
};

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
DEBILE_INLINE std::shared_ptr<mfem::HypreFGMRES> SolverHypreFGMRES::create_solver(HypreSolverType SOLVER,
                                                                    const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  auto ss = std::make_shared<mfem::HypreFGMRES>(MPI_COMM_WORLD);
  return ss;
}

