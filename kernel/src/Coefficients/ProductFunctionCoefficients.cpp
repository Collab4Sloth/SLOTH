/**
 * @file ProductFunctionCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class defining analytical expression, first and second derivatives of a Product of
 * FunctionCoefficient
 * @version 0.1
 * @date 2025-02-03
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
#include "Coefficients/ProductFunctionCoefficients.hpp"

#include <any>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"

/**
 * @brief Build a FunctionCoefficient as the product of several FunctionCoefficients
 *
 * @tparam Args
 *
 */
ProductCoefficient::ProductCoefficient(std::initializer_list<FunctionCoefficient*> coeffs)
    : vect_coefficients_(coeffs) {}

/**
 * @brief  Compute the value of the coefficient as the product of the values of the coefficients
 *
 * @return std::function<double(const std::span<const double>&, const std::span<const double>&,
 * const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
ProductCoefficient::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double prod_F = 1.0;
    for (const auto& coef : this->vect_coefficients_) {
      prod_F *= coef->eval_f(input_vector, exp_input_vector, auxiliary_vector, dimension);
    }

    return prod_F;
  };
  return func;
}

/**
 * @brief Gradient, by the product rule (no division by the factors, so that the gradient is
 * defined where a factor vanishes):
 *   dF/dx_i = sum_k (df_k/dx_i) prod_{l != k} f_l
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
ProductCoefficient::GradientF() {
  auto func = [this](const std::span<const double>& input_vector,
                     const std::span<const double>& exp_input_vector,
                     const std::span<const double>& auxiliary_vector,
                     const unsigned int dimension) {
    const std::size_t n = input_vector.size();
    const std::size_t nc = this->vect_coefficients_.size();

    // Values of the factors
    std::vector<double> f(nc);
    for (std::size_t k = 0; k < nc; ++k) {
      f[k] = this->vect_coefficients_[k]->eval_f(input_vector, exp_input_vector, auxiliary_vector,
                                                 dimension);
    }
    // Product of all the factors but k
    auto product_without = [&f, nc](std::size_t k) {
      double p = 1.0;
      for (std::size_t l = 0; l < nc; ++l) {
        if (l != k) p *= f[l];
      }
      return p;
    };

    std::vector<double> gradient(n, 0.0);
    for (std::size_t k = 0; k < nc; ++k) {
      const double others = product_without(k);
      for (std::size_t i = 0; i < n; ++i) {
        gradient[i] +=
            this->vect_coefficients_[k]->eval_gradient(
                static_cast<int>(i), input_vector, exp_input_vector, auxiliary_vector, dimension) *
            others;
      }
    }
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian, by the product rule (no division by the factors):
 *   d2F/dx_i dx_j = sum_k (d2f_k/dx_i dx_j) prod_{l != k} f_l
 *                 + sum_{k != m} (df_k/dx_i) (df_m/dx_j) prod_{l != k, m} f_l
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
ProductCoefficient::HessianF() {
  auto func = [this](const std::span<const double>& input_vector,
                     const std::span<const double>& exp_input_vector,
                     const std::span<const double>& auxiliary_vector,
                     const unsigned int dimension) {
    const std::size_t n = input_vector.size();
    const std::size_t nc = this->vect_coefficients_.size();

    // Values and gradients of the factors: f[k], g[k][i]
    std::vector<double> f(nc);
    std::vector<std::vector<double>> g(nc, std::vector<double>(n));
    for (std::size_t k = 0; k < nc; ++k) {
      const auto& coef = this->vect_coefficients_[k];
      f[k] = coef->eval_f(input_vector, exp_input_vector, auxiliary_vector, dimension);
      for (std::size_t i = 0; i < n; ++i) {
        g[k][i] = coef->eval_gradient(static_cast<int>(i), input_vector, exp_input_vector,
                                      auxiliary_vector, dimension);
      }
    }
    // Product of all the factors but k and m (m = k: all but k)
    auto product_without = [&f, nc](std::size_t k, std::size_t m) {
      double p = 1.0;
      for (std::size_t l = 0; l < nc; ++l) {
        if (l != k && l != m) p *= f[l];
      }
      return p;
    };

    std::vector<double> hessian(n * n, 0.0);
    for (std::size_t k = 0; k < nc; ++k) {
      const double others = product_without(k, k);
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          hessian[i * n + j] += this->vect_coefficients_[k]->eval_hessian(
                                    static_cast<int>(i), static_cast<int>(j), input_vector,
                                    exp_input_vector, auxiliary_vector, dimension) *
                                others;
        }
      }
      for (std::size_t m = 0; m < nc; ++m) {
        if (m == k) continue;
        const double others_km = product_without(k, m);
        for (std::size_t i = 0; i < n; ++i) {
          for (std::size_t j = 0; j < n; ++j) {
            hessian[i * n + j] += g[k][i] * g[m][j] * others_km;
          }
        }
      }
    }
    return hessian;
  };
  return func;
}

/**
 * @brief Set time for all FunctionCoefficients
 *
 * @param time
 */
void ProductCoefficient::set_time(double time) {
  this->time_ = time;

  for (auto* coef : vect_coefficients_) {
    coef->set_time(time);
  }
}

/**
 * @brief True if at least one FunctionCoefficient depends on the spatial coordinates
 */
bool ProductCoefficient::uses_coordinates() const {
  for (const auto* coef : vect_coefficients_) {
    if (coef->uses_coordinates()) return true;
  }
  return false;
}

/**
 * @brief Set the spatial coordinates for all FunctionCoefficients
 *
 * @param x Physical coordinates
 */
void ProductCoefficient::set_coordinates(std::span<const double> x) {
  this->coordinates_.assign(x.begin(), x.end());
  for (auto* coef : vect_coefficients_) {
    coef->set_coordinates(x);
  }
}

/**
 * @brief True if every term is a function of one variable only (then so is the product)
 */
bool ProductCoefficient::is_single_variable() const {
  if (vect_coefficients_.empty()) return false;
  for (const auto* coef : vect_coefficients_) {
    if (!coef->is_single_variable()) return false;
  }
  return true;
}

/**
 * @brief Number of variables of the evaluation by batch: the one of the terms, if all of them are
 * evaluated by batch with the same number of variables, 0 otherwise (point by point evaluation).
 */
int ProductCoefficient::batch_nb_variables() const {
  int n = 0;
  for (const auto* coef : vect_coefficients_) {
    const int nc = coef->batch_nb_variables();
    if (nc <= 0 || (n > 0 && nc != n)) return 0;
    n = nc;
  }
  return n;
}

/**
 * @brief Evaluation by batch (see FunctionCoefficient::eval_batch): product of the factors, each
 * one evaluated by batch, on the host or on the device.
 *
 * The product P of the first factors is multiplied by each new factor f (product rule):
 *   d2(P f) = d2P f + dP (x) df + df (x) dP + P d2f,   d(P f) = dP f + P df,   P f
 *
 * @return false if a factor is not evaluated by batch (point by point evaluation)
 */
bool ProductCoefficient::eval_batch(const int npts, const mfem::Vector& values,
                                    const mfem::Vector& explicit_values,
                                    const mfem::Vector& auxiliary_values,
                                    const mfem::Vector& coordinates, mfem::Vector* F,
                                    mfem::Vector* gradient, mfem::Vector* hessian) const {
  const int n = this->batch_nb_variables();
  if (n <= 0) return false;
  const bool want_hessian = hessian != nullptr;
  const bool need_gradient = gradient != nullptr || want_hessian;

  // Product of the first factors and its gradient: in the outputs if requested
  mfem::Vector& P = F != nullptr ? *F : this->product_F_;
  mfem::Vector& dP = gradient != nullptr ? *gradient : this->product_gradient_;
  P.SetSize(npts);
  if (need_gradient) dP.SetSize(npts * n);

  bool first = true;
  for (const auto* coef : vect_coefficients_) {
    if (first) {
      if (!coef->eval_batch(npts, values, explicit_values, auxiliary_values, coordinates, &P,
                            need_gradient ? &dP : nullptr, hessian)) {
        return false;
      }
      first = false;
      continue;
    }
    this->term_F_.SetSize(npts);
    if (need_gradient) this->term_gradient_.SetSize(npts * n);
    if (want_hessian) this->term_hessian_.SetSize(npts * n * n);
    if (!coef->eval_batch(npts, values, explicit_values, auxiliary_values, coordinates,
                          &this->term_F_, need_gradient ? &this->term_gradient_ : nullptr,
                          want_hessian ? &this->term_hessian_ : nullptr)) {
      return false;
    }

    // Data read and written where the loop runs (host or device)
    P.UseDevice(true);
    this->term_F_.UseDevice(true);
    double* p = P.ReadWrite();
    const double* f = this->term_F_.Read();
    double* dp = nullptr;
    const double* df = nullptr;
    if (need_gradient) {
      dP.UseDevice(true);
      this->term_gradient_.UseDevice(true);
      dp = dP.ReadWrite();
      df = this->term_gradient_.Read();
    }
    double* d2p = nullptr;
    const double* d2f = nullptr;
    if (want_hessian) {
      hessian->UseDevice(true);
      this->term_hessian_.UseDevice(true);
      d2p = hessian->ReadWrite();
      d2f = this->term_hessian_.Read();
    }
    mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
      const double pq = p[q];
      const double fq = f[q];
      // Second derivatives first: they use the gradient and the value of the previous product
      if (d2p != nullptr) {
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            const int ij = (q * n + i) * n + j;
            d2p[ij] = d2p[ij] * fq + dp[q * n + i] * df[q * n + j] + df[q * n + i] * dp[q * n + j] +
                      pq * d2f[ij];
          }
        }
      }
      if (dp != nullptr) {
        for (int i = 0; i < n; ++i) {
          dp[q * n + i] = dp[q * n + i] * fq + pq * df[q * n + i];
        }
      }
      p[q] = pq * fq;
    });
  }
  return true;
}
