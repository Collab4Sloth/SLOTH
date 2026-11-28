/**
 * @file AverageCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Average of a FunctionCoefficient over the segment [u_n, u_{n+1}] computed by a
 * quadrature rule (midpoint, trapezoid, average vector field...)
 * @version 0.1
 * @date 2026-11-27
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
#include "Coefficients/AverageCoefficient.hpp"

#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr const char* missing_explicit_values =
    "AverageCoefficient: the values at the beginning of the time step are required (the "
    "coefficient must be used with Scheme::SemiImplicit) and must have the same size as the "
    "values.";

/// @brief output = 0 (nothing if the output is not requested)
void reset(mfem::Vector* output) {
  if (output != nullptr) *output = 0.0;
}

/// @brief Work vector of the size of an output, nullptr if the output is not requested
mfem::Vector* work_output(const mfem::Vector* output, mfem::Vector& work) {
  if (output == nullptr) return nullptr;
  work.SetSize(output->Size());
  return &work;
}
}  // namespace

//////////////////////////////////////////////////////////////////////////////
// QuadRule
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief One point s = 1: value at u_{n+1}
 */
QuadRule QuadRule::end_point() { return QuadRule{{1.0}, {1.0}}; }

/**
 * @brief One point s = 1/2: value at (u_n + u_{n+1}) / 2
 */
QuadRule QuadRule::midpoint() { return QuadRule{{0.5}, {1.0}}; }

/**
 * @brief Two points s = 0 and s = 1, weights 1/2
 */
QuadRule QuadRule::trapezoid() { return QuadRule{{0.0, 1.0}, {0.5, 0.5}}; }

/**
 * @brief Gauss-Legendre rule on [0, 1] (MFEM rule of the reference segment), exact for
 * polynomials in s of degree lower than or equal to order
 *
 * @param order Degree of the polynomials integrated exactly (non negative)
 */
QuadRule QuadRule::gauss(const int order) {
  MFEM_VERIFY(order >= 0, "QuadRule::gauss: the order must be non negative.");
  const mfem::IntegrationRule& ir = mfem::IntRules.Get(mfem::Geometry::SEGMENT, order);
  QuadRule rule;
  rule.s.reserve(static_cast<std::size_t>(ir.GetNPoints()));
  rule.w.reserve(static_cast<std::size_t>(ir.GetNPoints()));
  for (int q = 0; q < ir.GetNPoints(); ++q) {
    rule.s.push_back(ir.IntPoint(q).x);
    rule.w.push_back(ir.IntPoint(q).weight);
  }
  return rule;
}

//////////////////////////////////////////////////////////////////////////////
// AverageCoefficient
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Construct the average of a coefficient over the segment [u_n, u_{n+1}]
 *
 * @param coefficient Averaged coefficient g (non-owning, not variational)
 * @param rule Quadrature rule on [0, 1]
 * @param average Quantity averaged: value or gradient of g (see Average)
 */
AverageCoefficient::AverageCoefficient(FunctionCoefficient* coefficient, QuadRule rule,
                                       const Average average)
    : coefficient_(coefficient), rule_(std::move(rule)), average_(average) {
  this->check();
}

/**
 * @brief Check the coefficient and the quadrature rule: points in [0, 1], positive weights of sum
 * equal to 1
 */
void AverageCoefficient::check() const {
  MFEM_VERIFY(this->coefficient_ != nullptr, "AverageCoefficient: null coefficient.");
  MFEM_VERIFY(!this->coefficient_->is_variational(),
              "AverageCoefficient: a variational coefficient is not supported.");
  MFEM_VERIFY(!this->rule_.s.empty() && this->rule_.s.size() == this->rule_.w.size(),
              "AverageCoefficient: the quadrature rule must have as many points as weights "
              "(at least one).");
  for (std::size_t q = 0; q < this->rule_.s.size(); ++q) {
    MFEM_VERIFY(this->rule_.s[q] >= 0.0 && this->rule_.s[q] <= 1.0,
                "AverageCoefficient: the points of the quadrature rule must be in [0, 1].");
    MFEM_VERIFY(this->rule_.w[q] > 0.0,
                "AverageCoefficient: the weights of the quadrature rule must be positive.");
  }
  const double sum = std::accumulate(this->rule_.w.begin(), this->rule_.w.end(), 0.0);
  MFEM_VERIFY(std::abs(sum - 1.0) < 1.e-12,
              "AverageCoefficient: the sum of the weights of the quadrature rule must be "
              "equal to 1.");
}

/**
 * @brief Weight of the first derivative of g at the q-th point: w_q s_q (value), w_q (gradient)
 */
double AverageCoefficient::gradient_weight(const std::size_t q) const {
  const double w = this->rule_.w[q];
  return this->average_ == Average::Value ? w * this->rule_.s[q] : w;
}

/**
 * @brief Weight of the second derivative of g at the q-th point: w_q s_q^2 (value), w_q s_q
 * (gradient)
 */
double AverageCoefficient::hessian_weight(const std::size_t q) const {
  const double ws = this->rule_.w[q] * this->rule_.s[q];
  return this->average_ == Average::Value ? ws * this->rule_.s[q] : ws;
}

/**
 * @brief q-th point of the segment: x = u_n + s_q (u_{n+1} - u_n)
 *
 * @param values u_{n+1}
 * @param exp_values u_n
 * @param q Index of the point of the rule
 * @param x Point of the segment (resized)
 */
void AverageCoefficient::segment_point(const std::span<const double>& values,
                                       const std::span<const double>& exp_values,
                                       const std::size_t q, std::vector<double>& x) const {
  MFEM_VERIFY(exp_values.size() == values.size(), missing_explicit_values);
  const double s = this->rule_.s[q];
  x.resize(values.size());
  for (std::size_t k = 0; k < values.size(); ++k) {
    x[k] = exp_values[k] + s * (values[k] - exp_values[k]);
  }
}

/**
 * @brief Value: average of g (Value), g(u_{n+1}) (Gradient)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
AverageCoefficient::F() {
  return [this](const std::span<const double>& input_vector,
                const std::span<const double>& exp_input_vector,
                const std::span<const double>& auxiliary_vector, const unsigned int dimension) {
    const int dim = static_cast<int>(dimension);
    if (this->average_ == Average::Gradient) {
      MFEM_VERIFY(exp_input_vector.size() == input_vector.size(), missing_explicit_values);
      return this->coefficient_->eval_f(input_vector, exp_input_vector, auxiliary_vector, dim);
    }
    std::vector<double> x;
    double value = 0.0;
    for (std::size_t q = 0; q < this->rule_.s.size(); ++q) {
      this->segment_point(input_vector, exp_input_vector, q, x);
      this->set_point_time(q);
      value +=
          this->rule_.w[q] * this->coefficient_->eval_f(x, exp_input_vector, auxiliary_vector, dim);
    }
    this->restore_time();
    return value;
  };
}

/**
 * @brief First derivatives with respect to u_{n+1} (see Average)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
AverageCoefficient::GradientF() {
  return [this](const std::span<const double>& input_vector,
                const std::span<const double>& exp_input_vector,
                const std::span<const double>& auxiliary_vector, const unsigned int dimension) {
    const int dim = static_cast<int>(dimension);
    const std::size_t n = input_vector.size();
    std::vector<double> gradient(n, 0.0);
    std::vector<double> x;
    for (std::size_t q = 0; q < this->rule_.s.size(); ++q) {
      this->segment_point(input_vector, exp_input_vector, q, x);
      this->set_point_time(q);
      const double a = this->gradient_weight(q);
      for (std::size_t i = 0; i < n; ++i) {
        gradient[i] += a * this->coefficient_->eval_gradient(
                               static_cast<int>(i), x, exp_input_vector, auxiliary_vector, dim);
      }
    }
    this->restore_time();
    return gradient;
  };
}

/**
 * @brief Second derivatives with respect to u_{n+1} (see Average)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
AverageCoefficient::HessianF() {
  return [this](const std::span<const double>& input_vector,
                const std::span<const double>& exp_input_vector,
                const std::span<const double>& auxiliary_vector, const unsigned int dimension) {
    const int dim = static_cast<int>(dimension);
    const std::size_t n = input_vector.size();
    std::vector<double> hessian(n * n, 0.0);
    std::vector<double> x;
    for (std::size_t q = 0; q < this->rule_.s.size(); ++q) {
      this->segment_point(input_vector, exp_input_vector, q, x);
      this->set_point_time(q);
      const double a = this->hessian_weight(q);
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          hessian[i * n + j] +=
              a * this->coefficient_->eval_hessian(static_cast<int>(i), static_cast<int>(j), x,
                                                   exp_input_vector, auxiliary_vector, dim);
        }
      }
    }
    this->restore_time();
    return hessian;
  };
}

/**
 * @brief Set the time of the coefficient and of g
 */
void AverageCoefficient::set_time(double time) {
  this->time_ = time;
  this->coefficient_->set_time(time);
}

/**
 * @brief Set the time step of the current step (average in time) for the coefficient and for g
 *
 * @param dt Time step (0: no average in time)
 */
void AverageCoefficient::set_time_step(double dt) {
  this->dt_ = dt;
  this->coefficient_->set_time_step(dt);
}

/**
 * @brief Time of g at the q-th point of the segment: t - (1 - s_q) dt (nothing if dt = 0)
 */
void AverageCoefficient::set_point_time(const std::size_t q) const {
  if (this->dt_ == 0.) return;
  this->coefficient_->set_time(this->time_ - (1. - this->rule_.s[q]) * this->dt_);
}

/**
 * @brief Time of g set back to the time of the coefficient (nothing if dt = 0)
 */
void AverageCoefficient::restore_time() const {
  if (this->dt_ == 0.) return;
  this->coefficient_->set_time(this->time_);
}

/**
 * @brief True if g depends on the spatial coordinates
 */
bool AverageCoefficient::uses_coordinates() const { return this->coefficient_->uses_coordinates(); }

/**
 * @brief True if g is a function of one variable only
 */
bool AverageCoefficient::is_single_variable() const {
  return this->coefficient_->is_single_variable();
}

/**
 * @brief Set the spatial coordinates of the coefficient and of g
 */
void AverageCoefficient::set_coordinates(std::span<const double> x) {
  this->coordinates_.assign(x.begin(), x.end());
  this->coefficient_->set_coordinates(x);
}

/**
 * @brief Number of variables of the evaluation by batch: the one of g
 */
int AverageCoefficient::batch_nb_variables() const {
  return this->coefficient_->batch_nb_variables();
}

/**
 * @brief Evaluation by batch (see FunctionCoefficient::eval_batch): g evaluated by batch at each
 * point of the rule, on the host or on the device, and the outputs accumulated with the weights
 * of Average
 *
 * @return false if g is not evaluated by batch (point by point evaluation)
 */
bool AverageCoefficient::eval_batch(const int npts, const mfem::Vector& values,
                                    const mfem::Vector& explicit_values,
                                    const mfem::Vector& auxiliary_values,
                                    const mfem::Vector& coordinates, mfem::Vector* F,
                                    mfem::Vector* gradient, mfem::Vector* hessian) const {
  if (this->batch_nb_variables() <= 0) return false;
  reset(F);
  reset(gradient);
  reset(hessian);
  if (npts == 0) return true;
  MFEM_VERIFY(explicit_values.Size() == values.Size(), missing_explicit_values);

  // Gradient: the value is the one of g at u_{n+1}
  const bool average_value = this->average_ == Average::Value;
  if (!average_value && F != nullptr) {
    if (!this->coefficient_->eval_batch(npts, values, explicit_values, auxiliary_values,
                                        coordinates, F, nullptr, nullptr)) {
      return false;
    }
  }

  // u_{n+1} - u_n, then the points u_n + s_q (u_{n+1} - u_n)
  this->delta_.SetSize(values.Size());
  this->points_.SetSize(values.Size());
  subtract(values, explicit_values, this->delta_);
  for (std::size_t q = 0; q < this->rule_.s.size(); ++q) {
    add(explicit_values, this->rule_.s[q], this->delta_, this->points_);
    this->set_point_time(q);
    mfem::Vector* q_F = average_value ? work_output(F, this->q_F_) : nullptr;
    mfem::Vector* q_gradient = work_output(gradient, this->q_gradient_);
    mfem::Vector* q_hessian = work_output(hessian, this->q_hessian_);
    if (!this->coefficient_->eval_batch(npts, this->points_, explicit_values, auxiliary_values,
                                        coordinates, q_F, q_gradient, q_hessian)) {
      this->restore_time();
      return false;
    }
    if (q_F != nullptr) F->Add(this->rule_.w[q], *q_F);
    if (q_gradient != nullptr) gradient->Add(this->gradient_weight(q), *q_gradient);
    if (q_hessian != nullptr) hessian->Add(this->hessian_weight(q), *q_hessian);
  }
  this->restore_time();
  return true;
}
