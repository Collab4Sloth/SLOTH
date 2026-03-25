/**
 * @file SlothSolver.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class to define a SlothSolver objet
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
#include "Solvers/SlothSolver.hpp"

#include <memory>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/DSolverBase.hpp"
#include "Solvers/HPrecondBase.hpp"
#include "Solvers/HSolverBase.hpp"
#include "Solvers/IPrecondBase.hpp"
#include "Solvers/ISolverBase.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

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
              return hh.create_solver(params_ + Parameter("description", "UMFPACK"));
            }
            default:
              mfem::mfem_error("Unhandled DirectSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, IterativeSolverType>) {
          switch (arg) {
            case IterativeSolverType::BICGSTAB: {
              SolverBICGSTAB hh;

              return hh.create_solver(params_ + Parameter("description", "BICGSTAB"));
            }
            case IterativeSolverType::CG: {
              SolverCG hh;
              return hh.create_solver(params_ + Parameter("description", "CG"));
            }
            case IterativeSolverType::GMRES: {
              SolverGMRES hh;
              return hh.create_solver(params_ + Parameter("description", "GMRES"));
            }
            case IterativeSolverType::MINRES: {
              SolverMINRES hh;
              return hh.create_solver(params_ + Parameter("description", "MINRES"));
            }
            default:
              mfem::mfem_error("Unhandled IterativeSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, HypreSolverType>) {
          switch (arg) {
            case HypreSolverType::HYPRE_PCG: {
              SolverHyprePCG hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_PCG"));
            }
            case HypreSolverType::HYPRE_GMRES: {
              SolverHypreGMRES hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_GMRES"));
            }
            case HypreSolverType::HYPRE_FGMRES: {
              SolverHypreFGMRES hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_FGMRES"));
            }
            default:
              mfem::mfem_error("Unhandled HypreSolverType enum value");
          }
        } else if constexpr (std::is_same_v<T, HyprePreconditionerType>) {
          switch (arg) {
            case HyprePreconditionerType::HYPRE_ILU: {
              PrecondHypreILU hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_ILU"));
            }
            case HyprePreconditionerType::HYPRE_BOOMER_AMG: {
              PrecondHypreBoomerAMG hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_BOOMER_AMG"));
            }
            case HyprePreconditionerType::HYPRE_DIAG_SCALE: {
              PrecondHypreDiagScale hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_DIAG_SCALE"));
            }
            case HyprePreconditionerType::HYPRE_SMOOTHER: {
              PrecondHypreSmoother hh;
              return hh.create_solver(params_ + Parameter("description", "HYPRE_SMOOTHER"));
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
              return hh.create_solver(params_ + Parameter("description", "SMOOTHER"));
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
