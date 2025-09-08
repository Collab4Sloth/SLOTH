/**
 * @file SolversOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Options for solvers and preconditioners
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
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

#include "Utils/Utils.hpp"

#pragma once

///////////////////////////////////////////////////
//////// SOLVERS
///////////////////////////////////////////////////
enum class NLSolverType { NEWTON };
enum class IterativeSolverType { BICGSTAB, GMRES, CG, MINRES };
enum class DirectSolverType { UMFPACK };
enum class HypreSolverType { HYPRE_PCG, HYPRE_GMRES, HYPRE_FGMRES };

///////////////////////////////////////////////////
//////// PRECONDITIONERS
///////////////////////////////////////////////////
enum class PreconditionerType { SMOOTHER, NO };
enum class HyprePreconditionerType {
  HYPRE_ILU,
  HYPRE_BOOMER_AMG,
  HYPRE_DIAG_SCALE,
  HYPRE_SMOOTHER,
  NO
};

//////////////////////////////////////////////////////
//// ALGORITHM
//////////////////////////////////////////////////////
/**
 * @brief Default constant used by Newton algorithm
 *
 */
namespace NewtonDefaultConstant {
const auto iter_max = 100;
const auto abs_tol = 1.e-13;
const auto rel_tol = 1.e-13;
const bool iterative_mode = false;
const auto print_level = 0;
}  // namespace NewtonDefaultConstant

//////////////////////////////////////////////////////
//// SOLVERS
//////////////////////////////////////////////////////
/**
 * @brief Default constant used by BICGSTAB Solver
 *
 */
namespace BICGSTABDefaultConstant {
const auto iter_max = 1000;
const auto abs_tol = 1.e-24;
const auto rel_tol = 1.e-12;
const bool iterative_mode = false;
const auto print_level = 0;
}  // namespace BICGSTABDefaultConstant

/**
 * @brief Default constant used by CG Solver
 *
 */
namespace CGDefaultConstant {
const auto iter_max = 1000;
const auto abs_tol = 1.e-24;
const auto rel_tol = 1.e-12;
const bool iterative_mode = false;
const auto print_level = 0;
}  // namespace CGDefaultConstant

/**
 * @brief Default constant used by MINRES  Solver
 *
 */
namespace MINRESDefaultConstant {
const auto iter_max = 1000;
const auto abs_tol = 1.e-24;
const auto rel_tol = 1.e-12;
const bool iterative_mode = false;
const auto print_level = 0;
}  // namespace MINRESDefaultConstant

/**
 * @brief Default constant used by GMRES Solver
 *
 */
namespace GMRESDefaultConstant {
const auto kdim = 50;
const auto iter_max = 1000;
const auto abs_tol = 1.e-24;
const auto rel_tol = 1.e-12;
const bool iterative_mode = false;
const auto print_level = 0;
}  // namespace GMRESDefaultConstant

/**
 * @brief Default constant used by Mass Solver
 *
 */
namespace UMFPACK_DefaultConstant {
const auto print_level = 0;
}  // namespace UMFPACK_DefaultConstant

//////////////////////////////////////////
//  PRECONDITIONERS
//////////////////////////////////////////
/**
 * @brief Default constant used by DSMOOTHER  Preconditioner
 *
 */
namespace DSMOOTHER_DefaultConstant {
const auto type = 0;  // 0, 1, 2 - scaled Jacobi, scaled l1-Jacobi, scaled lumped-Jacobi
const bool positive_diagonal = false;
const auto print_level = 0;
}  // namespace DSMOOTHER_DefaultConstant

//////////////////////////////////////////
// HYPRE SOLVERS
//////////////////////////////////////////

/**
 * @brief Default constant used by HYPRE_PCG Solver
 *
 */
namespace HYPRE_PCG_DefaultConstant {
const auto iter_max = 100;
const auto tol = 1.e-12;
// const auto abs_tol = 1.e-6;
const auto print_level = 0;
}  // namespace HYPRE_PCG_DefaultConstant

/**
 * @brief Default constant used by HYPRE_GMRES Solver
 *
 */
namespace HYPRE_GMRES_DefaultConstant {
const auto iter_max = 5000;
const auto tol = 1.e-12;
const auto kdim = 100;
const auto print_level = 0;
}  // namespace HYPRE_GMRES_DefaultConstant

/**
 * @brief Default constant used by HYPRE_FGMRES Solver
 *
 */
namespace HYPRE_FGMRES_DefaultConstant {
const auto iter_max = 5000;
const auto tol = 1.e-12;
const auto kdim = 100;
const auto print_level = 0;
}  // namespace HYPRE_FGMRES_DefaultConstant

//////////////////////////////////////////
// HYPRE PRECONDITIONERS
//////////////////////////////////////////
/**
 * @brief Default constant used by HYPRE_ILU  Preconditioner
 *
 */
namespace HYPRE_ILU_DefaultConstant {
const auto type = 0;  // ILU(k) locally and block Jacobi globally
const auto iter_max = 1;
const auto tol = 0.;
const auto print_level = 0;
const auto reorder_type = 0;  // 0 = no reordering, 1 = reverse Cuthill-McKee
}  // namespace HYPRE_ILU_DefaultConstant

/**
 * @brief Default constant used by HYPRE_SMOOTHER  Preconditioner
 *
 */
namespace HYPRE_SMOOTHER_DefaultConstant {
const auto type = 0;
// Taken from MFEM doc
//   Jacobi = 0,       ///< Jacobi
//   l1Jacobi = 1,     ///< l1-scaled Jacobi
//   l1GS = 2,         ///< l1-scaled block Gauss-Seidel/SSOR
//   l1GStr = 4,       ///< truncated l1-scaled block Gauss-Seidel/SSOR
//   lumpedJacobi = 5, ///< lumped Jacobi
//   GS = 6,           ///< Gauss-Seidel
//   OPFS = 10,        /**< On-processor forward solve for matrix w/ triangular
//                          structure */
//   Chebyshev = 16,   ///< Chebyshev
//   Taubin = 1001,    ///< Taubin polynomial smoother
//   FIR = 1002        ///< FIR polynomial smoother

const auto positive_diagonal = true;
}  // namespace HYPRE_SMOOTHER_DefaultConstant

/**
 * @brief Default constant used by HYPRE_BOOMER_AMG  Preconditioner
 *
 */
namespace HYPRE_BOOMER_AMG_DefaultConstant {
const auto iter_max = 100;
const auto tol = 1.e-16;
const auto print_level = 0;
}  // namespace HYPRE_BOOMER_AMG_DefaultConstant

/**
 * @brief Default constant used by HYPRE_DIAG_SCALE  Preconditioner
 *
 */
namespace HYPRE_DIAG_SCALE_DefaultConstant {
// No options
}  // namespace HYPRE_DIAG_SCALE_DefaultConstant
