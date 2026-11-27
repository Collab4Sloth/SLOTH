/**
 * @file BDF2Solver.cpp
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

#include "Solvers/BDF2Solver.hpp"

/**
 * @brief Coefficients {a, b, beta} of the variable-step BDF2 scheme
 *
 * @param w ratio dt_n / dt_{n-1}
 * @param a coefficient of x^n
 * @param b coefficient of x^{n-1}
 * @param beta coefficient of dt_n f(x^{n+1}, t^{n+1})
 */
void BDF2Solver::Coefficients(const mfem::real_t w, mfem::real_t& a, mfem::real_t& b,
                              mfem::real_t& beta) {
  a = (1. + w) * (1. + w) / (1. + 2. * w);
  b = w * w / (1. + 2. * w);
  beta = (1. + w) / (1. + 2. * w);
}

/**
 * @brief Associate the operator and restart the scheme (next step: backward Euler)
 *
 * @param f time-dependent operator implementing ImplicitSolve
 */
void BDF2Solver::Init(mfem::TimeDependentOperator& f) {
  mfem::ODESolver::Init(f);
  this->k_.SetSize(this->f->Width(), this->mem_type);
  this->k_ = 0.;
  this->x_star_.SetSize(this->f->Width(), this->mem_type);
  this->state_.SetSize(this->f->Width(), this->mem_type);
  this->state_.Reset();
  this->dt_prev_ = 0.;
}

/**
 * @brief One BDF2 step (backward Euler if there is no previous state)
 *
 * @param x state: x^n on input, x^{n+1} on output
 * @param t time: t^n on input, t^{n+1} on output
 * @param dt time step dt_n
 */
void BDF2Solver::Step(mfem::Vector& x, mfem::real_t& t, mfem::real_t& dt) {
  mfem::real_t a = 1., b = 0., beta = 1.;  // backward Euler
  const bool bdf2 = this->state_.Size() > 0 && this->dt_prev_ > 0.;
  if (bdf2) {
    MFEM_VERIFY(this->state_[0].Size() == x.Size(),
                "BDF2Solver: previous state and current state of different sizes.");
    const mfem::real_t w = dt / this->dt_prev_;
    MFEM_VERIFY(w < MaxRatio(), "BDF2Solver: time step ratio dt_n/dt_{n-1} = "
                                    << w << " >= 1 + sqrt(2) (zero-stability).");
    Coefficients(w, a, b, beta);
  }

  // x* = a x^n - b x^{n-1}
  this->x_star_.SetSize(x.Size(), this->mem_type);
  if (bdf2) {
    add(a, x, -b, this->state_[0], this->x_star_);
  } else {
    this->x_star_ = x;
  }

  // k = f(x* + beta dt k, t + dt): backward Euler problem with an effective state and step
  this->f->SetTime(t + dt);
  this->f->ImplicitSolve(beta * dt, this->x_star_, this->k_);

  // x^n becomes the previous state, x^{n+1} = x* + beta dt k
  this->state_.Append(x);
  add(this->x_star_, beta * dt, this->k_, x);
  this->dt_prev_ = dt;
  t += dt;
}
