/**
 * @file SlothSolver.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class to define a SlothSolver objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>

#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/DSolverBase.hpp"
#include "Solvers/HPrecondBase.hpp"
#include "Solvers/HSolverBase.hpp"
#include "Solvers/IPrecondBase.hpp"
#include "Solvers/ISolverBase.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief List of smart pointers toward MFEM solvers
 *
 */
using VSharedMFEMSolver =
    std::variant<std::shared_ptr<mfem::HyprePCG>, std::shared_ptr<mfem::HypreGMRES>,
                 std::shared_ptr<mfem::HypreSmoother>, std::shared_ptr<mfem::HypreFGMRES>,
                 std::shared_ptr<mfem::HypreBoomerAMG>, std::shared_ptr<mfem::HypreDiagScale>,
                 std::shared_ptr<mfem::HypreILU>, std::shared_ptr<mfem::UMFPackSolver>,
                 std::shared_ptr<mfem::BiCGSTABSolver>, std::shared_ptr<mfem::MINRESSolver>,
                 std::shared_ptr<mfem::CGSolver>, std::shared_ptr<mfem::GMRESSolver>,
                 std::shared_ptr<mfem::DSmoother>, std::shared_ptr<mfem::IterativeSolver>,
                 std::shared_ptr<mfem::Solver>, std::shared_ptr<std::monostate>>;

/**
 * @brief List of Solver and Preconditioner type
 *
 */
using VSolverType = std::variant<HypreSolverType, IterativeSolverType, DirectSolverType,
                                 PreconditionerType, HyprePreconditionerType>;

using VHypreSolver = std::variant<SolverHyprePCG, SolverHypreGMRES, SolverHypreFGMRES>;
using VIterativeSolver = std::variant<SolverBICGSTAB, SolverCG, SolverMINRES, SolverGMRES>;
using VDirectSolver = std::variant<SolverUMFPACK>;

using VHyprePrecond = std::variant<PrecondHypreILU, PrecondHypreBoomerAMG, PrecondHypreDiagScale>;
using VIterativePrecond = std::variant<PrecondDSmoother>;

/**
 * @brief Base class used to manage linear and non linear solvers
 *
 */
class SlothSolver {
 private:
  VSolverType value_;
  const Parameters& params_;

 public:
  SlothSolver(VSolverType value, const Parameters& params);

  VSharedMFEMSolver get_value();

  ~SlothSolver();
};
/**
 * @brief Construct a new SlothSolver object
 *
 * @param value : type of solver (see VSolverType for a list of available solvers)
 * @param params
 */
SlothSolver::SlothSolver(VSolverType value, const Parameters& params)
    : value_(value), params_(params) {}

/**
 * @brief Create the linear solver depending on SolverType
 *
 * @return VSharedMFEMSolver
 */
VSharedMFEMSolver SlothSolver::get_value() {
  return std::visit(
      [this](auto&& arg) -> VSharedMFEMSolver {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, DirectSolverType>) {
          switch (arg) {
            case DirectSolverType::UMFPACK: {
              SolverUMFPACK hh;
              return hh.create_solver(arg, params_);
            }
            default:
              throw std::runtime_error("Unhandled DirectSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, IterativeSolverType>) {
          switch (arg) {
            case IterativeSolverType::BICGSTAB: {
              SolverBICGSTAB hh;
              return hh.create_solver(arg, params_);
            }
            case IterativeSolverType::CG: {
              SolverCG hh;
              return hh.create_solver(arg, params_);
            }
            case IterativeSolverType::GMRES: {
              SolverGMRES hh;
              return hh.create_solver(arg, params_);
            }
            case IterativeSolverType::MINRES: {
              SolverMINRES hh;
              return hh.create_solver(arg, params_);
            }
            default:
              throw std::runtime_error("Unhandled IterativeSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, HypreSolverType>) {
          switch (arg) {
            case HypreSolverType::HYPRE_PCG: {
              SolverHyprePCG hh;
              return hh.create_solver(arg, params_);
            }
            case HypreSolverType::HYPRE_GMRES: {
              SolverHypreGMRES hh;
              return hh.create_solver(arg, params_);
            }
            case HypreSolverType::HYPRE_FGMRES: {
              SolverHypreFGMRES hh;
              return hh.create_solver(arg, params_);
            }
            default:
              throw std::runtime_error("Unhandled HypreSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, HyprePreconditionerType>) {
          switch (arg) {
            case HyprePreconditionerType::HYPRE_ILU: {
              PrecondHypreILU hh;
              return hh.create_solver(arg, params_);
            }
            case HyprePreconditionerType::HYPRE_BOOMER_AMG: {
              PrecondHypreBoomerAMG hh;
              return hh.create_solver(arg, params_);
            }
            case HyprePreconditionerType::HYPRE_DIAG_SCALE: {
              PrecondHypreDiagScale hh;
              return hh.create_solver(arg, params_);
            }
            case HyprePreconditionerType::HYPRE_SMOOTHER: {
              PrecondHypreSmoother hh;
              return hh.create_solver(arg, params_);
            }
            case HyprePreconditionerType::NO: {
              return std::make_shared<std::monostate>();
            }
            default:
              throw std::runtime_error("Unhandled HyprePreconditionerType enum value");
          }
        } else if constexpr (std::is_same_v<T, PreconditionerType>) {
          switch (arg) {
            case PreconditionerType::SMOOTHER: {
              PrecondDSmoother hh;
              return hh.create_solver(arg, params_);
            }
            case PreconditionerType::NO: {
              return std::make_shared<std::monostate>();
            }
            default:
              throw std::runtime_error("Unhandled PreconditionerType enum value");
          }
        } else {
          throw std::runtime_error("Unsupported type");
        }
      },
      value_);
}
/**
 * @brief Destroy the SlothSolver object
 *
 */
SlothSolver::~SlothSolver() {}
