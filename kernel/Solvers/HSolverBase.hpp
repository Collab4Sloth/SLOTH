
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
  SolverHyprePCG();
  std::shared_ptr<mfem::HyprePCG> create_solver(HypreSolverType SOLVER,
                                                const Parameters& params) override;

  ~SolverHyprePCG();
};

/**
 * @brief Construct a new HSolverBase::HSolverBase object
 *
 */
SolverHyprePCG::SolverHyprePCG() {}

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
std::shared_ptr<mfem::HyprePCG> SolverHyprePCG::create_solver(HypreSolverType SOLVER,
                                                              const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int iter_max =
      params.get_param_value_or_default<int>("iter_max", HYPRE_PCG_DefaultConstant::iter_max);
  const double tol =
      params.get_param_value_or_default<double>("tol", HYPRE_PCG_DefaultConstant::tol);
  const int print_level = params.get_param_value_or_default<double>(
      "print_level", HYPRE_PCG_DefaultConstant::print_level);

  auto ss = std::make_shared<mfem::HyprePCG>(MPI_COMM_WORLD);

  ss->SetMaxIter(iter_max);
  ss->SetTol(tol);
  ss->SetPrintLevel(print_level);

  return ss;
}

/**
 * @brief Destroy the HSolverBase::HSolverBase object
 *
 */
SolverHyprePCG::~SolverHyprePCG() {}

/////////////////////////////

class SolverHypreGMRES : public SolverBase<mfem::HypreGMRES, HypreSolverType> {
 public:
  SolverHypreGMRES();
  std::shared_ptr<mfem::HypreGMRES> create_solver(HypreSolverType SOLVER,
                                                  const Parameters& params) override;

  ~SolverHypreGMRES();
};

/**
 * @brief Construct a new HSolverBase::HSolverBase object
 *
 */
SolverHypreGMRES::SolverHypreGMRES() {}

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
std::shared_ptr<mfem::HypreGMRES> SolverHypreGMRES::create_solver(HypreSolverType SOLVER,
                                                                  const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int iter_max =
      params.get_param_value_or_default<int>("iter_max", HYPRE_GMRES_DefaultConstant::iter_max);
  const double tol =
      params.get_param_value_or_default<double>("tol", HYPRE_GMRES_DefaultConstant::tol);
  const int kdim =
      params.get_param_value_or_default<double>("kdim", HYPRE_GMRES_DefaultConstant::kdim);
  const int print_level = params.get_param_value_or_default<double>(
      "print_level", HYPRE_GMRES_DefaultConstant::print_level);

  auto ss = std::make_shared<mfem::HypreGMRES>(MPI_COMM_WORLD);

  ss->SetMaxIter(iter_max);
  ss->SetTol(tol);
  ss->SetKDim(kdim);
  ss->SetPrintLevel(print_level);

  return ss;
}

/**
 * @brief Destroy the HSolverBase::HSolverBase object
 *
 */
SolverHypreGMRES::~SolverHypreGMRES() {}

/////////////////////////////

class SolverHypreFGMRES : public SolverBase<mfem::HypreFGMRES, HypreSolverType> {
 public:
  SolverHypreFGMRES();
  std::shared_ptr<mfem::HypreFGMRES> create_solver(HypreSolverType SOLVER,
                                                   const Parameters& params) override;

  ~SolverHypreFGMRES();
};

/**
 * @brief Construct a new HSolverBase::HSolverBase object
 *
 */
SolverHypreFGMRES::SolverHypreFGMRES() {}

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
std::shared_ptr<mfem::HypreFGMRES> SolverHypreFGMRES::create_solver(HypreSolverType SOLVER,
                                                                    const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());
  const int iter_max =
      params.get_param_value_or_default<int>("iter_max", HYPRE_FGMRES_DefaultConstant::iter_max);
  const double tol =
      params.get_param_value_or_default<double>("tol", HYPRE_FGMRES_DefaultConstant::tol);
  const int kdim =
      params.get_param_value_or_default<double>("kdim", HYPRE_FGMRES_DefaultConstant::kdim);
  const int print_level = params.get_param_value_or_default<double>(
      "print_level", HYPRE_FGMRES_DefaultConstant::print_level);

  auto ss = std::make_shared<mfem::HypreFGMRES>(MPI_COMM_WORLD);
  ss->SetMaxIter(iter_max);
  ss->SetTol(tol);
  ss->SetKDim(kdim);
  ss->SetPrintLevel(print_level);

  return ss;
}

/**
 * @brief Destroy the HSolverBase::HSolverBase object
 *
 */
SolverHypreFGMRES::~SolverHypreFGMRES() {}
