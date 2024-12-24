/**
 * @file LSolver.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to manage Linear Solver  objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>

#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/SlothSolver.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Generic class for defining linear solvers with or without preconditionners
 *
 */
class LSolver {
 private:
  VSharedMFEMSolver variant_solver_;
  std::shared_ptr<SlothSolver> ss;

  VSharedMFEMSolver variant_precond_;
  std::shared_ptr<SlothSolver> pp;

 public:
  LSolver(VSolverType SOLVER, const Parameters& s_params, mfem::Operator& ope);
  LSolver(VSolverType SOLVER, const Parameters& s_params, VSolverType PRECOND,
          const Parameters& p_params, mfem::Operator& ope);
  VSharedMFEMSolver get_solver();

  ~LSolver();
};

struct SetPreconditioner
{
  // some stuff to reduce the verbosity
/*
  template<typename T>  using sptr = std::shared_ptr<mfem::DSmoother>;
  void warning(std::string solver, std::string precond)
  {
    std::string msg = "LSolver::LSolver Impossible to setup a " + solver + " solver with the preconditionner [" + precond + "]";
    SlothInfo::debug(msg);
  }

  template<typename Solv, typename Prec>
  void setter_mfem(Solv& solv, Prec& prec)
  {
    sptr<mfem::Solver> aaPrec =
    std::dynamic_pointer_cast<mfem::Solver>(prec);
    solv->SetPreconditioner(*aaPrec) ; 
  }
 
  template<typename Solv, typename Prec>
  void setter_hypre(Solv& solv, Prec& prec)
  {
    sptr<mfem::HypreSolver> aaPrec =
    std::dynamic_pointer_cast<mfem::HypreSolver>(prec);
    solv->SetPreconditioner(*aaPrec) ;
  }

  inline void operator()(auto&& solv, sptr<std::monostate>& prec) 
  { 
    SlothInfo::debug("LSolver:: monostate, it does nothing") ; 
  }

  void operator()(sptr<mfem::HypreGMRES>&  solv, sptr<mfem::DSmoother>& prec) { warning("HypreGMRES" , "DSmoother]"); }
  void operator()(sptr<mfem::HyprePCG>&    solv, sptr<mfem::DSmoother>& prec) { warning("HyprePCG"   , "DSmoother]"); }
  void operator()(sptr<mfem::HypreFGMRES>& solv, sptr<mfem::DSmoother>& prec) { warning("HypreFGMRES", "DSmoother]"); }

  template<typename Solv>
  inline void operator()(Solv& solv, sptr<mfem::DSmoother>& prec)
  {
    SlothInfo::debug("LSolver::LSolver setting iterative preconditionner [DSmoother]");
    setter_mfem(solv, prec);
  }

  template<typename Solv>
  inline void operator()(Solv&& solv, sptr<mfem::HypreSmoother>&& prec)
  {
    SlothInfo::debug("LSolver::LSolver setting iterative preconditionner [HypreSmoother]");
    setter_hypre(solv, prec);
  }

  template<typename Solv>
  inline void operator()(Solv&& solv, sptr<mfem::HypreILU>&& prec)
  {
    SlothInfo::debug("LSolver::LSolver setting hypre preconditionner [HypreILU]");
    setter_hypre(solv, prec);
  }

  template<typename Solv>
  inline void operator()(Solv&& solv, sptr<mfem::HypreBoomerAMG>&& prec)
  {
    SlothInfo::debug("LSolver::LSolver setting hypre preconditionner [HypreBoomerAMG]");
    setter_hypre(solv, prec);
  }

  template<typename Solv>
  inline void operator()(Solv&& solv, sptr<mfem::HypreDiagScale>&& prec)
  {
    SlothInfo::debug("LSolver::LSolver setting hypre preconditionner [HypreDiagScale]");
    setter_hypre(solv, prec);
  }
*/
  // member
  mfem::Operator& op;

  /* main function */
  template<typename Solv, typename Prec>
  inline void operator()(Solv&& solv, Prec&& prec) 
  {
    using TT = std::decay_t<decltype(solv)>;
    using PP = std::decay_t<decltype(prec)>;
    if constexpr (is_in_variant_v<TT, VIterativeSolver>) 
    {
      if constexpr (!std::is_same_v<PP, std::shared_ptr<std::monostate>>) 
      {
        MFEM_VERIFY((is_in_variant_v<PP, VIterativePrecond>),
        "LSolver:: IterativeSolver objects  can only be associated with an "
        "IterativePreconditionner objects");
        if constexpr (is_in_variant_v<PP, VIterativePrecond>) 
        {
          SlothInfo::debug("LSolver::LSolver setting iterative preconditionner");
          std::shared_ptr<mfem::Solver> aaPrec = std::dynamic_pointer_cast<mfem::Solver>(prec);
          solv->SetPreconditioner(*aaPrec);
        }
      }
    }
    if constexpr (is_in_variant_v<TT, VHypreSolver>) 
    {
      if constexpr (!std::is_same_v<PP, std::shared_ptr<std::monostate>>) 
      {
        MFEM_VERIFY((is_in_variant_v<PP, VHyprePrecond>),
        "LSolver:: HypreSolver objects  can only be associated with an "
        "HyprePreconditionner objects");
        if constexpr (is_in_variant_v<PP, VHyprePrecond>) 
        {
          SlothInfo::debug("LSolver::LSolver setting hypre preconditionner (not hypre smoother)");
          std::shared_ptr<mfem::HypreSolver> aaPrec = std::dynamic_pointer_cast<mfem::HypreSolver>(prec);
          solv->SetPreconditioner(*aaPrec);
        }
      }
    }
    if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) 
    {
      SlothInfo::debug("LSolver::LSolver setting operator ");
      solv->SetOperator(this->op);
    }
  };
};

/**
 * @brief Construct a new LSolver::LSolver object with a preconditionner
 *
 * @param SOLVER
 * @param s_params
 * @param PRECOND
 * @param p_params
 * @param ope
 */
LSolver::LSolver(VSolverType SOLVER, const Parameters& s_params, VSolverType PRECOND,
                 const Parameters& p_params, mfem::Operator& ope) {
  SlothInfo::debug("LSolver::LSolver start");
  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();

  pp = std::make_shared<SlothSolver>(PRECOND, p_params);
  this->variant_precond_ = pp->get_value();

  SetPreconditioner func_prec = {ope};
  std::visit(func_prec, this->variant_solver_, this->variant_precond_);
  SlothInfo::debug("LSolver::LSolver end");
}

/**
 * @brief Construct a new LSolver::LSolver object without precondnitioner
 *
 * @param SOLVER
 * @param s_params
 * @param ope
 */
LSolver::LSolver(VSolverType SOLVER, const Parameters& s_params, mfem::Operator& ope) {
  SlothInfo::debug("LSolver::LSolver start");

  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();
  std::visit(
      [&](auto&& arg) {
        using TT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
          SlothInfo::debug("LSolver::LSolver setting operator ");
          arg->SetOperator(ope);
        }
      },
      this->variant_solver_);
  SlothInfo::debug("LSolver::LSolver end");
}

/**
 * @brief Return the solver
 *
 * @return std::shared_ptr<T>
 */
VSharedMFEMSolver LSolver::get_solver() { return this->variant_solver_; }

/**
 * @brief Destroy the LSolver::LSolver object
 *
 */
LSolver::~LSolver() {}
