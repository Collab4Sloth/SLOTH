/**
 * @file SumFunctionCoefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class defining analytical expression, first and second derivatives of a sum of
 * FunctionCoefficient objects
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
#include "Coefficients/SumFunctionCoefficients.hpp"

#include <any>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"

/**
 * @brief Build a FunctionCoefficient as the summation of several FunctionCoefficients
 *
 * @tparam Args
 *
 */
SumCoefficient::SumCoefficient(std::initializer_list<FunctionCoefficient*> coeffs)
    : vect_coefficients_(coeffs) {}

/**
 * @brief  Compute the value of the coefficient as a summation of the value of coefficients
 *
 * @return std::function<double(const std::span<const double>&, const std::span<const double>&,
 * const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
SumCoefficient::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double sum_F = 0.0;
    for (const auto& coef : this->vect_coefficients_) {
      sum_F += coef->eval_f(input_vector, exp_input_vector, auxiliary_vector, dimension);
    }

    return sum_F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
SumCoefficient::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(input_vector.size(), 0.0);

    for (unsigned int i = 0; i < input_vector.size(); ++i) {
      for (const auto& coef : this->vect_coefficients_) {
        gradient[i] +=
            coef->eval_gradient(i, input_vector, exp_input_vector, auxiliary_vector, dimension);
      }
    }

    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
SumCoefficient::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    const int size = input_vector.size();
    std::vector<double> hessian(size * size, 0.0);

    for (unsigned int i = 0; i < input_vector.size(); ++i) {
      for (unsigned int j = 0; j < input_vector.size(); ++j) {
        for (const auto& coef : this->vect_coefficients_) {
          hessian[i * size + j] +=
              coef->eval_hessian(i, j, input_vector, exp_input_vector, auxiliary_vector, dimension);
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
void SumCoefficient::set_time(double time) {
  this->time_ = time;

  for (auto* coef : vect_coefficients_) {
    coef->set_time(time);
  }
}

/**
 * @brief True if at least one FunctionCoefficient depends on the spatial coordinates
 */
bool SumCoefficient::uses_coordinates() const {
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
void SumCoefficient::set_coordinates(std::span<const double> x) {
  this->coordinates_.assign(x.begin(), x.end());
  for (auto* coef : vect_coefficients_) {
    coef->set_coordinates(x);
  }
}

/**
 * @brief True if every term is a function of one variable only (then so is the sum)
 */
bool SumCoefficient::is_single_variable() const {
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
int SumCoefficient::batch_nb_variables() const {
  int n = 0;
  for (const auto* coef : vect_coefficients_) {
    const int nc = coef->batch_nb_variables();
    if (nc <= 0 || (n > 0 && nc != n)) return 0;
    n = nc;
  }
  return n;
}

/**
 * @brief Evaluation by batch (see FunctionCoefficient::eval_batch): sum of the outputs of the
 * terms, each one evaluated by batch, on the host or on the device.
 *
 * @return false if a term is not evaluated by batch (point by point evaluation)
 */
bool SumCoefficient::eval_batch(const int npts, const mfem::Vector& values,
                                const mfem::Vector& explicit_values,
                                const mfem::Vector& auxiliary_values,
                                const mfem::Vector& coordinates, mfem::Vector* F,
                                mfem::Vector* gradient, mfem::Vector* hessian) const {
  if (this->batch_nb_variables() <= 0) return false;
  bool first = true;
  for (const auto* coef : vect_coefficients_) {
    if (first) {
      // First term: written in the outputs
      if (!coef->eval_batch(npts, values, explicit_values, auxiliary_values, coordinates, F,
                            gradient, hessian)) {
        return false;
      }
      first = false;
      continue;
    }
    if (F != nullptr) this->term_F_.SetSize(F->Size());
    if (gradient != nullptr) this->term_gradient_.SetSize(gradient->Size());
    if (hessian != nullptr) this->term_hessian_.SetSize(hessian->Size());
    if (!coef->eval_batch(npts, values, explicit_values, auxiliary_values, coordinates,
                          F != nullptr ? &this->term_F_ : nullptr,
                          gradient != nullptr ? &this->term_gradient_ : nullptr,
                          hessian != nullptr ? &this->term_hessian_ : nullptr)) {
      return false;
    }
    if (F != nullptr) *F += this->term_F_;
    if (gradient != nullptr) *gradient += this->term_gradient_;
    if (hessian != nullptr) *hessian += this->term_hessian_;
  }
  return true;
}
