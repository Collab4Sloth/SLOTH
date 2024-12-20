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
#include "Utils/UtilsForVariants.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief List of Solver and Preconditioner type
 *
 */
using VSolverType = std::variant<HypreSolverType, IterativeSolverType, DirectSolverType,
                                 PreconditionerType, HyprePreconditionerType>;

//-------------------
// Solvers
//-------------------
using VHypreSolver =
    std::variant<std::shared_ptr<mfem::HyprePCG>, std::shared_ptr<mfem::HypreGMRES>,
                 std::shared_ptr<mfem::HypreFGMRES>>;
using VIterativeSolver =
    std::variant<std::shared_ptr<mfem::BiCGSTABSolver>, std::shared_ptr<mfem::MINRESSolver>,
                 std::shared_ptr<mfem::CGSolver>, std::shared_ptr<mfem::GMRESSolver>>;
using VDirectSolver = std::variant<std::shared_ptr<mfem::UMFPackSolver>>;

using VSolvers = concat_variant_type<VHypreSolver, VIterativeSolver, VDirectSolver>;

//-------------------
// Preconditionners
//-------------------
using VHyprePrecond =
    std::variant<std::shared_ptr<mfem::HypreILU>, std::shared_ptr<mfem::HypreBoomerAMG>,
                 std::shared_ptr<mfem::HypreDiagScale>>;
using VIterativePrecond =
    std::variant<std::shared_ptr<mfem::DSmoother>, std::shared_ptr<mfem::HypreSmoother>>;
using VPreconds = concat_variant_type<VHyprePrecond, VIterativePrecond>;

/**
 * @brief List of smart pointers toward MFEM solvers
 *
 */
using VSharedMFEMSolver =
    concat_variant_type<VSolvers, VPreconds, std::variant<std::shared_ptr<std::monostate>>>;

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
              return hh.create_solver(arg, params_ + Parameter("description", "UMFPACK"));
            }
            default:
              mfem::mfem_error("Unhandled DirectSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, IterativeSolverType>) {
          switch (arg) {
            case IterativeSolverType::BICGSTAB: {
              SolverBICGSTAB hh;

              return hh.create_solver(arg, params_ + Parameter("description", "BICGSTAB"));
            }
            case IterativeSolverType::CG: {
              SolverCG hh;
              return hh.create_solver(arg, params_ + Parameter("description", "CG"));
            }
            case IterativeSolverType::GMRES: {
              SolverGMRES hh;
              return hh.create_solver(arg, params_ + Parameter("description", "GMRES"));
            }
            case IterativeSolverType::MINRES: {
              SolverMINRES hh;
              return hh.create_solver(arg, params_ + Parameter("description", "MINRES"));
            }
            default:
              mfem::mfem_error("Unhandled IterativeSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, HypreSolverType>) {
          switch (arg) {
            case HypreSolverType::HYPRE_PCG: {
              SolverHyprePCG hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_PCG"));
            }
            case HypreSolverType::HYPRE_GMRES: {
              SolverHypreGMRES hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_GMRES"));
            }
            case HypreSolverType::HYPRE_FGMRES: {
              SolverHypreFGMRES hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_FGMRES"));
            }
            default:
              mfem::mfem_error("Unhandled HypreSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, HyprePreconditionerType>) {
          switch (arg) {
            case HyprePreconditionerType::HYPRE_ILU: {
              PrecondHypreILU hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_ILU"));
            }
            case HyprePreconditionerType::HYPRE_BOOMER_AMG: {
              PrecondHypreBoomerAMG hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_BOOMER_AMG"));
            }
            case HyprePreconditionerType::HYPRE_DIAG_SCALE: {
              PrecondHypreDiagScale hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_DIAG_SCALE"));
            }
            case HyprePreconditionerType::HYPRE_SMOOTHER: {
              PrecondHypreSmoother hh;
              return hh.create_solver(arg, params_ + Parameter("description", "HYPRE_SMOOTHER"));
            }
            case HyprePreconditionerType::NO: {
              return std::make_shared<std::monostate>();
            }
            default:
              mfem::mfem_error("Unhandled HyprePreconditionerType enum value");
          }
        } else if constexpr (std::is_same_v<T, PreconditionerType>) {
          switch (arg) {
            case PreconditionerType::SMOOTHER: {
              PrecondDSmoother hh;
              return hh.create_solver(arg, params_ + Parameter("description", "SMOOTHER"));
            }
            case PreconditionerType::NO: {
              return std::make_shared<std::monostate>();
            }
            default:
              mfem::mfem_error("Unhandled PreconditionerType enum value");
          }
        } else {
          mfem::mfem_error("Unsupported type");
        }
      },
      value_);
}
/**
 * @brief Destroy the SlothSolver object
 *
 */
SlothSolver::~SlothSolver() {}
