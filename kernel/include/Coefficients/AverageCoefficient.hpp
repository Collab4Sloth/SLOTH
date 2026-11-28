/**
 * @file AverageCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Average of a FunctionCoefficient over the segment [u_n, u_{n+1}] computed by a
 * quadrature rule (midpoint, trapezoid, average vector field...)
 * @version 0.1
 * @date 2026-11-27
 *
 * Copyright CEA (C) 2026
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

#include <functional>
#include <span>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Quadrature rule on [0, 1] for an average over the segment [u_n, u_{n+1}]:
 *   int_0^1 g(u_n + s (u_{n+1} - u_n)) ds ~ sum_q w_q g(u_n + s_q (u_{n+1} - u_n))
 */
struct QuadRule {
  std::vector<double> s;  ///< points of the rule, in [0, 1]
  std::vector<double> w;  ///< weights of the rule, sum equal to 1

  /// One point s = 1: value at u_{n+1}
  static QuadRule end_point();
  /// One point s = 1/2: value at (u_n + u_{n+1}) / 2
  static QuadRule midpoint();
  /// Two points s = 0 and s = 1, weights 1/2: (g(u_n) + g(u_{n+1})) / 2
  static QuadRule trapezoid();
  /// Gauss-Legendre rule, exact for polynomials in s of degree lower than or equal to order
  static QuadRule gauss(const int order);
};

/**
 * @brief Quantity averaged over the segment by a AverageCoefficient g_avg built on g
 *
 * With x_q = u_n + s_q (u_{n+1} - u_n), the derivatives being taken with respect to u_{n+1}:
 *
 * - Value: average of the value of g (mobility, coefficient depending on the unknowns)
 *     g_avg        = sum_q w_q g(x_q)
 *     dg_avg       = sum_q w_q s_q   dg(x_q)
 *     d2g_avg      = sum_q w_q s_q^2 d2g(x_q)
 *
 * - Gradient: average of the gradient of g (free energy: discrete gradient). With an exact rule,
 *   sum_i dg_avg_i (u_{n+1} - u_n)_i = g(u_{n+1}) - g(u_n) (average vector field).
 *     g_avg        = g(u_{n+1})            (free energy at the end of the time step)
 *     dg_avg       = sum_q w_q dg(x_q)
 *     d2g_avg      = sum_q w_q s_q d2g(x_q)  (derivative of dg_avg with respect to u_{n+1})
 */
enum class Average { Value, Gradient };

/**
 * @brief Average of a FunctionCoefficient g over the segment [u_n, u_{n+1}]
 *
 * u_{n+1} are the values given to the coefficient and u_n its explicit values: it must be used in
 * a Coefficient with Scheme::SemiImplicit, which gives both. The explicit values must have the
 * same size as the values. The explicit values are also given unchanged to g.
 *
 * Examples:
 * - mobility at the midpoint: AverageCoefficient(&M, QuadRule::midpoint(),
 *   Average::Value);
 * - average vector field of a free energy W of degree p: AverageCoefficient(&W,
 *   QuadRule::gauss(p - 1), Average::Gradient).
 *
 * @remark Non-owning: g must outlive the coefficient. A variational coefficient is not supported.
 * @remark Average in time: the point s of the segment is evaluated at the time t - (1 - s) dt,
 * t being the time of the coefficient (end of the step) and dt the time step given by
 * set_time_step (no shift while no time step is given). g is then set back to the time t.
 */
class AverageCoefficient : public FunctionCoefficient {
 private:
  FunctionCoefficient* coefficient_;  ///< averaged coefficient g
  QuadRule rule_;
  Average average_;
  double dt_{0.};  ///< time step of the current step (see set_time_step)

  // Work arrays of eval_batch
  mutable mfem::Vector delta_, points_;                ///< u_{n+1} - u_n, points of the segment
  mutable mfem::Vector q_F_, q_gradient_, q_hessian_;  ///< outputs of g at one point of the rule

  void check() const;
  double gradient_weight(const std::size_t q) const;
  double hessian_weight(const std::size_t q) const;
  void set_point_time(const std::size_t q) const;
  void restore_time() const;
  void segment_point(const std::span<const double>& values,
                     const std::span<const double>& exp_values, const std::size_t q,
                     std::vector<double>& x) const;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  AverageCoefficient(FunctionCoefficient* coefficient, QuadRule rule, const Average average);
  virtual ~AverageCoefficient() = default;

  const QuadRule& rule() const { return this->rule_; }
  Average average() const { return this->average_; }

  void set_time(double time) override;
  void set_time_step(double dt) override;
  bool uses_coordinates() const override;
  bool is_single_variable() const override;
  void set_coordinates(std::span<const double> x) override;

  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const override;
  int batch_nb_variables() const override;
};
