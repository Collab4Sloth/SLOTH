/**
 * @file SlothBDF2Solver.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief BDF2 for SLOTH problems: history from the variables, extrapolation of the lagged
 *        auxiliary variables of a partitioned coupling
 * @version 0.1
 * @date 2026-10-30
 *
 * @copyright CEA (C) 2026
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
#pragma once

#include <utility>
#include <vector>

#include "Solvers/BDF2Solver.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief An auxiliary variable is lagged (seen at t^n by a problem of a partitioned coupling) if
 *        it has been advanced as many times as the variables of the problem: the problem owning
 *        it has not been solved yet in the current step. A variable advanced once more belongs
 *        to a problem already solved in the step; a variable never advanced is a fixed field.
 *
 * @param auxvar auxiliary variable
 * @param vars variables of the problem
 */
template <class T, int DIM>
bool is_lagged_auxiliary_variable(const Variable<T, DIM>& auxvar, Variables<T, DIM>& vars) {
  if (vars.get_variables_number() == 0) return false;
  const int own_updates = vars[0].get_number_of_updates();
  // First step: every variable has 0 update (fixed fields cannot be distinguished)
  return own_updates > 0 && auxvar.get_number_of_updates() == own_updates;
}

/**
 * @brief BDF2 (BDF2Solver) for a SLOTH problem.
 *
 * In addition to BDF2Solver:
 *  - the previous state u^{n-1} is taken from the history of the variables of the problem
 *    (Variable::get_second_to_last), which is transferred by the AMR with the current solution:
 *    the scheme is not restarted after a mesh change (Init keeps the previous time step);
 *  - partitioned coupling: the lagged auxiliary variables (is_lagged_auxiliary_variable) are
 *    temporarily replaced, during the step, by the second-order extrapolation
 *      B~ = (1 + w) B^n - w B^{n-1},  w = dt_n / dt_{n-1},
 *    truncated in their extrapolation bounds if any (Variable::set_extrapolation_bounds), and
 *    restored after the step (also if the step throws). Without it, the lag is a splitting error
 *    O(dt) that limits the order to 1.
 */
template <class T, int DIM>
class SlothBDF2Solver : public BDF2Solver {
 private:
  Variables<T, DIM>* vars_{nullptr};
  std::vector<Variables<T, DIM>*> auxvars_;
  mfem::Vector u_prev_;
  bool info_printed_{false};

  /// Restores the extrapolated auxiliary variables when it goes out of scope
  class ExtrapolationGuard {
    std::vector<Variable<T, DIM>*> extrapolated_;

   public:
    ExtrapolationGuard() = default;
    ExtrapolationGuard(const ExtrapolationGuard&) = delete;
    ExtrapolationGuard& operator=(const ExtrapolationGuard&) = delete;
    void add(Variable<T, DIM>* var) { this->extrapolated_.push_back(var); }
    ~ExtrapolationGuard() {
      for (auto* var : this->extrapolated_) var->restore_state();
    }
  };

  /// u^{n-1}: second-to-last entry of the history of the variables, in the block layout of x
  void set_previous_state_from_variables(const int size) {
    this->u_prev_.SetSize(size);
    int offset = 0;
    for (size_t i = 0; i < this->vars_->get_variables_number(); i++) {
      const mfem::Vector prev = (*this->vars_)[i].get_second_to_last();
      MFEM_VERIFY(offset + prev.Size() <= size,
                  "SlothBDF2Solver: history of the variables inconsistent with the state.");
      for (int j = 0; j < prev.Size(); j++) this->u_prev_(offset + j) = prev(j);
      offset += prev.Size();
    }
    MFEM_VERIFY(offset == size,
                "SlothBDF2Solver: history of the variables inconsistent with the state.");
    auto& state = this->GetState();
    if (state.Size() == 0) {
      state.Append(this->u_prev_);
    } else {
      state.Set(0, this->u_prev_);
    }
  }

  /// Extrapolation of the lagged auxiliary variables (restored by the guard)
  void extrapolate_lagged_auxiliary_variables(const mfem::real_t w, ExtrapolationGuard& guard) {
    const bool print_info = !this->info_printed_ && mfem::Mpi::WorldRank() == 0;
    for (auto* auxvar_vec : this->auxvars_) {
      if (auxvar_vec == nullptr) continue;
      for (size_t i = 0; i < auxvar_vec->get_variables_number(); i++) {
        auto& auxvar = (*auxvar_vec)[i];
        if (!is_lagged_auxiliary_variable(auxvar, *this->vars_)) continue;
        mfem::Vector u_tilde(auxvar.get_last());
        u_tilde *= (1. + w);
        u_tilde.Add(-w, auxvar.get_second_to_last());
        auxvar.set_extrapolated_state(u_tilde);
        guard.add(&auxvar);
        if (print_info) {
          SlothInfo::verbose(" SlothBDF2Solver: auxiliary variable '", auxvar.getVariableName(),
                             "' not yet advanced in the step (partitioned coupling): extrapolated "
                             "to the end of the step");
        }
      }
    }
    this->info_printed_ = true;
  }

 public:
  SlothBDF2Solver() = default;

  /// Variables of the problem and auxiliary variables (TransientOperator::initialize)
  void set_variables(Variables<T, DIM>& vars, std::vector<Variables<T, DIM>*> auxvars) {
    this->vars_ = &vars;
    this->auxvars_ = std::move(auxvars);
  }

  /// Associate the operator (also after a mesh change): the history comes from the variables,
  /// so the scheme is not restarted (the previous time step is kept)
  void Init(mfem::TimeDependentOperator& f) override {
    const mfem::real_t dt_prev = this->GetPreviousTimeStep();
    BDF2Solver::Init(f);
    this->SetPreviousTimeStep(dt_prev);
  }

  void Step(mfem::Vector& x, mfem::real_t& t, mfem::real_t& dt) override {
    MFEM_VERIFY(this->vars_ != nullptr, "SlothBDF2Solver: set_variables must be called first.");
    ExtrapolationGuard guard;
    const mfem::real_t dt_prev = this->GetPreviousTimeStep();
    if (dt_prev > 0.) {
      this->set_previous_state_from_variables(x.Size());
      this->extrapolate_lagged_auxiliary_variables(dt / dt_prev, guard);
    }
    BDF2Solver::Step(x, t, dt);
  }  // guard: auxiliary variables restored
};
