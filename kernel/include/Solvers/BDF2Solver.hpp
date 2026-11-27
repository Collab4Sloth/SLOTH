/**
 * @file BDF2Solver.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Variable-step BDF2 time integrator (mfem::ODESolverWithStates)
 * @version 0.2
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

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Variable-step BDF2 scheme written as a backward Euler step with an effective state and
 *        an effective time step.
 *
 * With w = dt_n / dt_{n-1}, the variable-step BDF2 scheme
 *   x^{n+1} - a x^n + b x^{n-1} = beta dt_n f(x^{n+1}, t^{n+1}),
 *   a = (1+w)^2 / (1+2w),  b = w^2 / (1+2w),  beta = (1+w) / (1+2w),
 * is solved for k = f(x^{n+1}, t^{n+1}) with
 *   x^{n+1} = x* + beta dt_n k,   x* = a x^n - b x^{n-1},
 * that is TimeDependentOperator::ImplicitSolve(beta dt_n, x*, k): any operator implementing
 * ImplicitSolve for backward Euler can be used without modification.
 *
 * Depends on MFEM only. The previous state x^{n-1} is stored as the state of the solver
 * (ODESolverWithStates::GetState(), one vector), as for the multistep solvers of MFEM. It can be
 * set by the user (GetState().Set / Append), for instance after a mesh change, together with
 * the previous time step (SetPreviousTimeStep).
 *
 * - The first step after Init() (no previous state) is a backward Euler step.
 * - Init() restarts the scheme (MFEM convention).
 * - Zero-stability of variable-step BDF2 requires w < 1 + sqrt(2) for arbitrary step sequences
 *   (M. Crouzeix, F. J. Lisbona, The convergence of variable-stepsize, variable-formula,
 *   multistep methods, SIAM J. Numer. Anal. 21(3), 1984): larger ratios are rejected. This bound
 *   is general, hence conservative: larger ratios are known to be stable for parabolic problems.
 */
class BDF2Solver : public mfem::ODESolverWithStates {
 private:
  mfem::Vector k_;
  mfem::Vector x_star_;
  mfem::ODEStateDataVector state_{1};  // x^{n-1}
  mfem::real_t dt_prev_{0.};

 public:
  /// Maximal ratio dt_n / dt_{n-1}: 1 + sqrt(2) (std::sqrt is not constexpr before C++26)
  static constexpr mfem::real_t MaxRatio() { return 1. + 1.4142135623730951; }

  /// Coefficients {a, b, beta} of the variable-step BDF2 scheme for the ratio w = dt_n / dt_{n-1}
  static void Coefficients(const mfem::real_t w, mfem::real_t& a, mfem::real_t& b,
                           mfem::real_t& beta);

  BDF2Solver() = default;

  /// Associate the operator and restart the scheme (next step: backward Euler)
  void Init(mfem::TimeDependentOperator& f) override;

  /// Previous time step dt_{n-1} (0 if there is no previous step)
  mfem::real_t GetPreviousTimeStep() const { return this->dt_prev_; }
  /// Set the previous time step (with a user-provided previous state, e.g. after a mesh change)
  void SetPreviousTimeStep(const mfem::real_t dt_prev) { this->dt_prev_ = dt_prev; }

  mfem::ODEStateData& GetState() override { return this->state_; }
  const mfem::ODEStateData& GetState() const override { return this->state_; }

  void Step(mfem::Vector& x, mfem::real_t& t, mfem::real_t& dt) override;
};
