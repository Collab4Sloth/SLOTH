/**
 * @file SlothNewtonSolver.hpp
 * @brief mfem::NewtonSolver that stops as soon as the residual norm diverges.
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
#include <cmath>
#include <sstream>

#include "Time/TimeStepFailure.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Newton solver with divergence detection.
 *
 * mfem::NewtonSolver only stops on convergence or after iter_max iterations. This class checks the
 * residual norm at each iteration (ComputeScalingFactor is called once per iteration, with the
 * residual of the current iterate stored in r) and throws a TimeStepFailure as soon as
 *   ||r_k|| > divergence_factor * ||r_0||   or   ||r_k|| is not finite.
 * With adaptive time-stepping, the time-step is then rejected and reduced.
 *
 * A divergence_factor <= 0 disables the check. It is disabled for steady problems (no time-step to
 * reject, and the residual of a Newton started far from the solution is often non-monotonic).
 */
class SlothNewtonSolver : public mfem::NewtonSolver {
 private:
  double divergence_factor_;
  mutable double initial_norm_ = -1.;
  mutable int iteration_ = 0;

 public:
  SlothNewtonSolver(MPI_Comm comm, const double divergence_factor)
      : mfem::NewtonSolver(comm), divergence_factor_(divergence_factor) {}

  /// Divergence factor (<= 0: no divergence detection)
  void set_divergence_factor(const double factor) { this->divergence_factor_ = factor; }

  void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
    this->initial_norm_ = -1.;
    this->iteration_ = 0;
    mfem::NewtonSolver::Mult(b, x);
  }

  mfem::real_t ComputeScalingFactor(const mfem::Vector& x, const mfem::Vector& b) const override {
    if (this->divergence_factor_ > 0.) {
      // r holds F(x) - b for the current iterate (global norm: same value on every rank)
      const double norm = this->Norm(this->r);
      if (this->initial_norm_ < 0.) {
        this->initial_norm_ = norm;
      } else if (!std::isfinite(norm) || norm > this->divergence_factor_ * this->initial_norm_) {
        std::ostringstream msg;
        msg << "Newton solver diverged at iteration " << this->iteration_ << " (||r|| = " << norm
            << ", ||r_0|| = " << this->initial_norm_ << ", divergence_factor "
            << this->divergence_factor_ << ").";
        throw TimeStepFailure(msg.str());
      }
      ++this->iteration_;
    }
    return mfem::NewtonSolver::ComputeScalingFactor(x, b);
  }
};
