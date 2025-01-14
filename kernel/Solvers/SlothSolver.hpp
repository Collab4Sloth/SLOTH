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

#pragma once

/**
 * @brief Alias for shared_ptr
 *
 * @tparam T
 */
template <typename T>
using sptr = std::shared_ptr<T>;
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
    std::variant<sptr<mfem::HyprePCG>, sptr<mfem::HypreGMRES>, sptr<mfem::HypreFGMRES>>;
using VIterativeSolver = std::variant<sptr<mfem::BiCGSTABSolver>, sptr<mfem::MINRESSolver>,
                                      sptr<mfem::CGSolver>, sptr<mfem::GMRESSolver>>;
using VDirectSolver = std::variant<sptr<mfem::UMFPackSolver>>;

using VSolvers = concat_variant_type<VHypreSolver, VIterativeSolver, VDirectSolver>;

//-------------------
// Preconditionners
//-------------------
using VHyprePrecond =
    std::variant<sptr<mfem::HypreILU>, sptr<mfem::HypreBoomerAMG>, sptr<mfem::HypreDiagScale>>;
using VIterativePrecond = std::variant<sptr<mfem::DSmoother>, sptr<mfem::HypreSmoother>>;
using VPreconds = concat_variant_type<VHyprePrecond, VIterativePrecond>;

/**
 * @brief List of smart pointers toward MFEM solvers
 *
 */
using VSharedMFEMSolver =
    concat_variant_type<VSolvers, VPreconds, std::variant<sptr<std::monostate>>>;

/**
 * @brief Common methods for Linear and non Linear solvers
 *
 */
struct UtilsSolvers {
  /**
   * @brief Set the mfem Solver object
   *
   * @tparam Solv
   * @tparam Prec
   * @param solv
   * @param prec
   */
  template <typename Solv, typename Prec>
  static void setter_mfem(Solv& solv, Prec& prec) {
    sptr<mfem::Solver> aaPrec = std::dynamic_pointer_cast<mfem::Solver>(prec);
    solv->SetPreconditioner(*aaPrec);
  }

  /**
   * @brief Set the Hypre Solver object
   *
   * @tparam Solv
   * @tparam Prec
   * @param solv
   * @param prec
   */
  template <typename Solv, typename Prec>
  static void setter_hypre(Solv& solv, Prec& prec) {
    sptr<mfem::HypreSolver> aaPrec = std::dynamic_pointer_cast<mfem::HypreSolver>(prec);
    solv->SetPreconditioner(*aaPrec);
  }

  /**
   * @brief Print the name of the Solver
   *
   * @tparam Solv
   */
  template <typename Solv>
  static void print_solver_name(const Solv&) {
    if constexpr (std::is_same<Solv, sptr<mfem::HyprePCG>>::value)
      SlothInfo::debug("Solver used: HyprePCG ");
    else if constexpr (std::is_same<Solv, sptr<mfem::HypreGMRES>>::value)
      SlothInfo::debug("Solver used: HypreGMRES ");
    else if constexpr (std::is_same<Solv, sptr<mfem::HypreFGMRES>>::value)
      SlothInfo::debug("Solver used: HypreFGMRES ");
    else if constexpr (std::is_same<Solv, sptr<mfem::BiCGSTABSolver>>::value)
      SlothInfo::debug("Solver used: BiCGSTABSolver ");
    else if constexpr (std::is_same<Solv, sptr<mfem::MINRESSolver>>::value)
      SlothInfo::debug("Solver used: MINRESSolver ");
    else if constexpr (std::is_same<Solv, sptr<mfem::CGSolver>>::value)
      SlothInfo::debug("Solver used: CGSolver ");
    else if constexpr (std::is_same<Solv, sptr<mfem::GMRESSolver>>::value)
      SlothInfo::debug("Solver used: GMRESSolver ");
    else
      SlothInfo::debug("Solver used: unknown ");
  }

  /**
   * @brief Print the name of the Preconditioner
   *
   * @tparam Prec
   */
  template <typename Prec>
  static void print_prec_name(const Prec&) {
    if constexpr (std::is_same<Prec, sptr<mfem::HypreILU>>::value)
      SlothInfo::debug("Preconditioner used: HypreILU ");
    else if constexpr (std::is_same<Prec, sptr<mfem::HypreBoomerAMG>>::value)
      SlothInfo::debug("Preconditioner used: HypreBoomerAMG ");
    else if constexpr (std::is_same<Prec, sptr<mfem::HypreDiagScale>>::value)
      SlothInfo::debug("Preconditioner used: HypreDiagScale ");
    else if constexpr (std::is_same<Prec, sptr<mfem::HypreSmoother>>::value)
      SlothInfo::debug("Preconditioner used: HypreSmoother ");
    else if constexpr (std::is_same<Prec, sptr<mfem::DSmoother>>::value)
      SlothInfo::debug("Preconditioner used: DSmoother ");
    else if constexpr (std::is_same<Prec, sptr<std::monostate>>::value)
      SlothInfo::debug("No Preconditioner used ");
    else
      SlothInfo::debug("Preconditioner used: unknown ");
  }
};

/**
 * @brief Define a linear solver and its preconditioner (if required)
 *
 */
struct SetPrecondSolver {
  // member
  mfem::Operator& op;

  /**
   * @brief Main function
   *
   * @tparam Solv
   * @tparam Prec
   * @param solv
   * @param prec
   */
  template <typename Solv, typename Prec>
  inline void operator()(Solv&& solv, Prec&& prec) {
    using TT = std::decay_t<decltype(solv)>;
    using PP = std::decay_t<decltype(prec)>;
    UtilsSolvers::print_solver_name(solv);
    UtilsSolvers::print_prec_name(prec);
    if constexpr (!std::is_same_v<PP, sptr<std::monostate>>) {
      if constexpr (is_in_variant_v<TT, VIterativeSolver>) {
        MFEM_VERIFY((is_in_variant_v<PP, VIterativePrecond>),
                    "SetPrecondSolver: IterativeSolver objects  can only be associated with an "
                    "IterativePreconditionner objects");
        if constexpr (is_in_variant_v<PP, VIterativePrecond>) {
          SlothInfo::debug("SetPrecondSolver: setting iterative preconditionner");
          UtilsSolvers::setter_mfem(solv, prec);
        }
      }
      if constexpr (is_in_variant_v<TT, VHypreSolver>) {
        MFEM_VERIFY((is_in_variant_v<PP, VHyprePrecond>),
                    "SetPrecondSolver: HypreSolver objects can only be associated with an "
                    "HyprePreconditionner objects");
        if constexpr (is_in_variant_v<PP, VHyprePrecond>) {
          SlothInfo::debug("SetPrecondSolver: setting hypre preconditionner (not hypre smoother)");
          UtilsSolvers::setter_hypre(solv, prec);
        }
      }
    }

    if constexpr (!std::is_same_v<TT, sptr<std::monostate>>) {
      SlothInfo::debug("SetPrecondSolver: setting operator ");
      solv->SetOperator(this->op);
    }
  }
};

/**
 * @brief Define a non linear solver and its preconditioner (if required)
 *
 */
struct SetPrecondNLSolver {
  // members
  mfem::Operator& op;
  sptr<mfem::NewtonSolver> nl_solver;

  /**
   * @brief Main function
   *
   * @tparam Solv
   * @tparam Prec
   * @param solv
   * @param prec
   */
  template <typename Solv, typename Prec>
  inline void operator()(Solv&& solv, Prec&& prec) {
    using TT = std::decay_t<decltype(solv)>;
    using PP = std::decay_t<decltype(prec)>;
    UtilsSolvers::print_solver_name(solv);
    UtilsSolvers::print_prec_name(prec);
    if constexpr (!std::is_same_v<PP, sptr<std::monostate>>) {
      if constexpr (is_in_variant_v<TT, VIterativeSolver>) {
        MFEM_VERIFY((is_in_variant_v<PP, VIterativePrecond>),
                    "SetPrecondNLSolver: IterativeSolver objects  can only be associated with an "
                    "IterativePreconditionner objects");
        if constexpr (is_in_variant_v<PP, VIterativePrecond>) {
          SlothInfo::debug("SetPrecondNLSolver: setting iterative preconditionner");
          UtilsSolvers::setter_mfem(solv, prec);
        }
      }
      if constexpr (is_in_variant_v<TT, VHypreSolver>) {
        MFEM_VERIFY((is_in_variant_v<PP, VHyprePrecond>),
                    "SetPrecondNLSolver: HypreSolver objects  can only be associated with an "
                    "HyprePreconditionner objects");
        if constexpr (is_in_variant_v<PP, VHyprePrecond>) {
          SlothInfo::debug(
              "SetPrecondNLSolver: setting hypre preconditionner (not hypre smoother)");
          UtilsSolvers::setter_hypre(solv, prec);
        }
      }
    }

    SlothInfo::debug("SetPrecondNLSolver: setting operator ");
    this->nl_solver->SetOperator(this->op);
    if constexpr (!std::is_same_v<TT, sptr<std::monostate>>) {
      SlothInfo::debug("SetPrecondNLSolver: setting solver ");
      this->nl_solver->SetSolver(*solv);
    }
    SlothInfo::debug("SetPrecondNLSolver: end");
  }
};

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
