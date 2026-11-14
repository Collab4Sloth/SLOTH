/**
 * @file CommonCoefficients.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief List of common coefficients
 * @version 0.1
 * @date 2026-03-21
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

#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

#include "Coefficients/EvalAtPoints.hpp"
#include "Coefficients/FunctionCoefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Evaluation of W at one point, on the host or on the device
 *
 * Same expressions as W (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct WEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 0;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    [[maybe_unused]] const double x = values[0];
    if (F != nullptr) {
      *F = prefactor_ * (pow(x, 2) * pow(1 - x, 2));
    }
    if (gradient != nullptr) {
      gradient[0] = prefactor_ * (pow(x, 2) * (2 * x - 2) + 2 * x * pow(1 - x, 2));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ * (2 * pow(x, 2) + 4 * x * (2 * x - 2) + 2 * pow(1 - x, 2));
    }
  }
};

/**
 *
 * @brief Coefficient based on expression: x*x*(1.0-x) * (1.0-x)
 *
 */
class W : public FunctionCoefficient {
 private:
  WEval eval_;  ///< F, gradient and Hessian (host and device)

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
  bool is_single_variable() const override { return true; }
  W() = default;
  explicit W(const double prefactor) { eval_.prefactor_ = prefactor; }
  virtual ~W() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const WEval& get_eval() const { return eval_; }
  /// Evaluation at a set of points with the evaluation structure (see eval_at_points)
  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const final {
    eval_at_points(eval_, npts, this->time_, values, explicit_values, auxiliary_values, coordinates,
                   F, gradient, hessian);
    return true;
  }
  int batch_nb_variables() const final { return WEval::nb_variables; }
};

/**
 * @brief Evaluation of Fw at one point, on the host or on the device
 *
 * Same expressions as Fw (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct FwEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 0;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    [[maybe_unused]] const double x = values[0];
    if (F != nullptr) {
      *F = prefactor_ * (((1.0 / 4.0) * pow(x, 2) - 1.0 / 4.0) * (pow(x, 2) - 1));
    }
    if (gradient != nullptr) {
      gradient[0] = prefactor_ * (2 * x * ((1.0 / 4.0) * pow(x, 2) - 1.0 / 4.0) +
                                  (1.0 / 2.0) * x * (pow(x, 2) - 1));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ * (3 * pow(x, 2) - 1);
    }
  }
};

/**
 *
 * @brief Coefficient based on expression: 0.25 * (x * x - 1.0) * (x * x - 1.0)
 *
 */
class Fw : public FunctionCoefficient {
 private:
  FwEval eval_;  ///< F, gradient and Hessian (host and device)

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
  bool is_single_variable() const override { return true; }
  Fw() = default;
  explicit Fw(const double prefactor) { eval_.prefactor_ = prefactor; }

  virtual ~Fw() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const FwEval& get_eval() const { return eval_; }
  /// Evaluation at a set of points with the evaluation structure (see eval_at_points)
  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const final {
    eval_at_points(eval_, npts, this->time_, values, explicit_values, auxiliary_values, coordinates,
                   F, gradient, hessian);
    return true;
  }
  int batch_nb_variables() const final { return FwEval::nb_variables; }
};

/**
 * @brief Evaluation of H at one point, on the host or on the device
 *
 * Same expressions as H (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct HEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 0;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    [[maybe_unused]] const double x = values[0];
    if (F != nullptr) {
      *F = prefactor_ * (pow(x, 3) * (6 * pow(x, 2) - 15 * x + 10));
    }
    if (gradient != nullptr) {
      gradient[0] =
          prefactor_ * (pow(x, 3) * (12 * x - 15) + 3 * pow(x, 2) * (6 * pow(x, 2) - 15 * x + 10));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ * (12 * pow(x, 3) + 6 * pow(x, 2) * (12 * x - 15) +
                                 6 * x * (6 * pow(x, 2) - 15 * x + 10));
    }
  }
};

/**
 *
 * @brief Coefficient based on expression: x * x * x * (6.0 * x * x - 15.0 * x + 10.0)
 *
 */
class H : public FunctionCoefficient {
 private:
  HEval eval_;  ///< F, gradient and Hessian (host and device)

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
  bool is_single_variable() const override { return true; }
  H() = default;
  explicit H(const double prefactor) { eval_.prefactor_ = prefactor; }

  virtual ~H() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const HEval& get_eval() const { return eval_; }
  /// Evaluation at a set of points with the evaluation structure (see eval_at_points)
  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const final {
    eval_at_points(eval_, npts, this->time_, values, explicit_values, auxiliary_values, coordinates,
                   F, gradient, hessian);
    return true;
  }
  int batch_nb_variables() const final { return HEval::nb_variables; }
};

/**
 * @brief Evaluation of Log at one point, on the host or on the device
 *
 * Same expressions as Log (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 *
 * Bounded variables (F, gradient and Hessian evaluated at the clamped values):
 *   x in [eps, (1.0 - eps)]
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct LogEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 0;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;
  double eps = 1.e-10;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    // Bounded variables: clamped to their bounds
    [[maybe_unused]] const double x = fmin((1.0 - eps), fmax(eps, values[0]));
    if (F != nullptr) {
      *F = prefactor_ * (x * log(x) + (1 - x) * log(1 - x));
    }
    if (gradient != nullptr) {
      gradient[0] = prefactor_ * (log(x) - log(1 - x));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ * (1.0 / (1 - x) + 1.0 / x);
    }
  }
};

/**
 *
 * @brief Coefficient based on expression: x * log(x) + (1.0 - x)*log(1.0 - x)
 *
 */
class Log : public FunctionCoefficient {
 private:
  LogEval eval_;  ///< F, gradient and Hessian (host and device)

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
  bool is_single_variable() const override { return true; }
  Log() = default;
  explicit Log(const double prefactor) { eval_.prefactor_ = prefactor; }

  virtual ~Log() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const LogEval& get_eval() const { return eval_; }
  /// Evaluation at a set of points with the evaluation structure (see eval_at_points)
  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const final {
    eval_at_points(eval_, npts, this->time_, values, explicit_values, auxiliary_values, coordinates,
                   F, gradient, hessian);
    return true;
  }
  int batch_nb_variables() const final { return LogEval::nb_variables; }
};

/**
 *
 * @brief Coefficient based on expression: 0.5*dot(x,x)
 *
 */
class GradientEnergy : public FunctionCoefficient {
 private:
  double prefactor_;

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
  GradientEnergy() : prefactor_(1.0) {}
  explicit GradientEnergy(const double prefactor) : prefactor_(prefactor) {}

  virtual ~GradientEnergy() = default;
};
