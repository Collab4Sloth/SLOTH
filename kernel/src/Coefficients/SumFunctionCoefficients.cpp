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

#include <algorithm>
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

namespace {
constexpr const char* explicit_term_error =
    "SumCoefficient: a term with Scheme::Explicit is evaluated with the values at the beginning of "
    "the time step: the sum must be wrapped in a Coefficient with Scheme::SemiImplicit.";

/// @brief Values given to a term: the explicit ones for an explicit term
const std::span<const double>& term_input(const bool explicit_term,
                                          const std::span<const double>& input_vector,
                                          const std::span<const double>& exp_input_vector) {
  if (!explicit_term) return input_vector;
  MFEM_VERIFY(!exp_input_vector.empty(), explicit_term_error);
  return exp_input_vector;
}

/// @brief Work vector of the size of an output, nullptr if the output is not requested
mfem::Vector* term_output(const mfem::Vector* output, mfem::Vector& work) {
  if (output == nullptr) return nullptr;
  work.SetSize(output->Size());
  return &work;
}

/// @brief output += term (nothing if the output is not requested)
void accumulate(mfem::Vector* output, const mfem::Vector* term) {
  if (output != nullptr && term != nullptr) *output += *term;
}

/// @brief output = 0 (nothing if the output is not requested)
void reset(mfem::Vector* output) {
  if (output != nullptr) *output = 0.0;
}
}  // namespace

/**
 * @brief Build a FunctionCoefficient as the summation of several FunctionCoefficients
 *
 * @tparam Args
 *
 */
SumCoefficient::SumCoefficient(std::initializer_list<FunctionCoefficient*> coeffs)
    : vect_coefficients_(coeffs), vect_schemes_(coeffs.size(), Scheme::Implicit) {
  this->check_terms();
}

/**
 * @brief Build a FunctionCoefficient as the summation of several FunctionCoefficients, each one
 * with its own scheme (see SumTerm)
 *
 * @param terms Terms of the sum with their scheme
 */
SumCoefficient::SumCoefficient(std::initializer_list<SumTerm> terms) {
  for (const auto& term : terms) {
    this->vect_coefficients_.push_back(term.coefficient);
    this->vect_schemes_.push_back(term.scheme);
  }
  this->check_terms();
}

/**
 * @brief Check the terms: not null, scheme Implicit, Explicit or SemiImplicit
 */
void SumCoefficient::check_terms() const {
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    MFEM_VERIFY(this->vect_coefficients_[k] != nullptr, "SumCoefficient: null term.");
    const Scheme scheme = this->vect_schemes_[k];
    MFEM_VERIFY(
        scheme == Scheme::Implicit || scheme == Scheme::Explicit || scheme == Scheme::SemiImplicit,
        "SumCoefficient: the scheme of a term must be Implicit, Explicit or SemiImplicit.");
  }
}

/**
 * @brief True if the k-th term is explicit
 */
bool SumCoefficient::is_explicit_term(const std::size_t k) const {
  return this->vect_schemes_[k] == Scheme::Explicit;
}

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
    for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
      const auto& input = term_input(this->is_explicit_term(k), input_vector, exp_input_vector);
      sum_F +=
          this->vect_coefficients_[k]->eval_f(input, exp_input_vector, auxiliary_vector, dimension);
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

    for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
      const auto& input = term_input(this->is_explicit_term(k), input_vector, exp_input_vector);
      for (unsigned int i = 0; i < input_vector.size(); ++i) {
        gradient[i] += this->vect_coefficients_[k]->eval_gradient(i, input, exp_input_vector,
                                                                  auxiliary_vector, dimension);
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

    // Explicit terms: derivatives at the beginning of the time step, no second derivative
    for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
      if (this->is_explicit_term(k)) continue;
      for (unsigned int i = 0; i < input_vector.size(); ++i) {
        for (unsigned int j = 0; j < input_vector.size(); ++j) {
          hessian[i * size + j] += this->vect_coefficients_[k]->eval_hessian(
              i, j, input_vector, exp_input_vector, auxiliary_vector, dimension);
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
 * @brief Set the time step for all FunctionCoefficients
 *
 * @param dt Time step of the current step
 */
void SumCoefficient::set_time_step(double dt) {
  for (auto* coef : vect_coefficients_) {
    coef->set_time_step(dt);
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
 * terms, each one evaluated by batch, on the host or on the device. An explicit term is evaluated
 * with the explicit values and has no second derivative.
 *
 * @return false if a term is not evaluated by batch (point by point evaluation)
 */
bool SumCoefficient::eval_batch(const int npts, const mfem::Vector& values,
                                const mfem::Vector& explicit_values,
                                const mfem::Vector& auxiliary_values,
                                const mfem::Vector& coordinates, mfem::Vector* F,
                                mfem::Vector* gradient, mfem::Vector* hessian) const {
  if (this->batch_nb_variables() <= 0) return false;
  reset(F);
  reset(gradient);
  reset(hessian);
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    const bool explicit_term = this->is_explicit_term(k);
    MFEM_VERIFY(!explicit_term || npts == 0 || explicit_values.Size() > 0, explicit_term_error);
    mfem::Vector* term_F = term_output(F, this->term_F_);
    mfem::Vector* term_gradient = term_output(gradient, this->term_gradient_);
    mfem::Vector* term_hessian =
        explicit_term ? nullptr : term_output(hessian, this->term_hessian_);
    if (!this->vect_coefficients_[k]->eval_batch(npts, explicit_term ? explicit_values : values,
                                                 explicit_values, auxiliary_values, coordinates,
                                                 term_F, term_gradient, term_hessian)) {
      return false;
    }
    accumulate(F, term_F);
    accumulate(gradient, term_gradient);
    accumulate(hessian, term_hessian);
  }
  return true;
}

/**
 * @brief True if at least one term is variational: the sum is then evaluated by
 * eval_variational_batch only
 */
bool SumCoefficient::is_variational() const {
  return std::ranges::any_of(this->vect_coefficients_,
                             [](const FunctionCoefficient* c) { return c->is_variational(); });
}

/**
 * @brief Dimension of the gradients of the variational terms, 0 if there is none
 */
int SumCoefficient::variational_dimension() const {
  for (const auto* coef : this->vect_coefficients_) {
    if (coef->is_variational()) return coef->variational_dimension();
  }
  return 0;
}

/**
 * @brief Number of variables of the variational terms, 0 if there is none
 */
int SumCoefficient::variational_nb_variables() const {
  for (const auto* coef : this->vect_coefficients_) {
    if (coef->is_variational()) return coef->variational_nb_variables();
  }
  return 0;
}

/**
 * @brief True if a variational term reads the gradients of its explicit variables or is explicit
 * (evaluated with the gradients at the beginning of the time step)
 */
bool SumCoefficient::uses_explicit_gradients() const {
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    const auto* coef = this->vect_coefficients_[k];
    if (!coef->is_variational()) continue;
    if (coef->uses_explicit_gradients() || this->is_explicit_term(k)) return true;
  }
  return false;
}

/**
 * @brief True if a variational term reads the gradients of its auxiliary variables
 */
bool SumCoefficient::uses_auxiliary_gradients() const {
  return std::ranges::any_of(this->vect_coefficients_, [](const FunctionCoefficient* c) {
    return c->is_variational() && c->uses_auxiliary_gradients();
  });
}

/**
 * @brief True if the sum is variational and all its variational terms are isotropic and not
 * explicit (non-variational terms: no gradient term)
 */
bool SumCoefficient::is_isotropic() const {
  if (!this->is_variational()) return false;
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    const auto* coef = this->vect_coefficients_[k];
    if (!coef->is_variational()) continue;
    if (this->is_explicit_term(k) || !coef->is_isotropic()) return false;
  }
  return true;
}

/**
 * @brief Evaluation of a variational sum (see FunctionCoefficient::eval_variational_batch)
 *
 * - variational term: all the outputs;
 * - non-variational term, function of the variables of the sum: F, df/du and d2f/du2 only;
 * - explicit term: evaluated with the explicit values (and the explicit gradients for a
 *   variational term), no second derivative.
 *
 * @return false if the sum is not variational or if a term is not evaluated by batch
 */
bool SumCoefficient::eval_variational_batch(
    const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
    const mfem::Vector& auxiliary_values, const mfem::Vector& gradients,
    const mfem::Vector& explicit_gradients, const mfem::Vector& auxiliary_gradients,
    const mfem::Vector& coordinates, const VariationalBatch& outputs) const {
  if (!this->is_variational()) return false;
  const int n = this->variational_nb_variables();
  const int d = this->variational_dimension();
  for (const auto* coef : this->vect_coefficients_) {
    if (coef->is_variational()) {
      MFEM_VERIFY(coef->variational_nb_variables() == n && coef->variational_dimension() == d,
                  "SumCoefficient: the variational terms must have the same variables and the "
                  "same dimension.");
    } else {
      MFEM_VERIFY(coef->batch_nb_variables() == n,
                  "SumCoefficient: a non-variational term of a variational sum must be a function "
                  "of the " +
                      std::to_string(n) + " variables of the sum, evaluated by batch.");
    }
  }

  reset(outputs.F);
  reset(outputs.df_dvar);
  reset(outputs.df_dgrad);
  reset(outputs.d2f_dvar2);
  reset(outputs.d2f_dvar_dgrad);
  reset(outputs.d2f_dgrad2);
  reset(outputs.d2f_dgrad2_isotropic);
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    const auto* coef = this->vect_coefficients_[k];
    const bool explicit_term = this->is_explicit_term(k);
    MFEM_VERIFY(!explicit_term || npts == 0 || explicit_values.Size() > 0, explicit_term_error);
    const mfem::Vector& term_values = explicit_term ? explicit_values : values;
    VariationalBatch term;
    term.F = term_output(outputs.F, this->term_F_);
    term.df_dvar = term_output(outputs.df_dvar, this->term_gradient_);
    if (!explicit_term) term.d2f_dvar2 = term_output(outputs.d2f_dvar2, this->term_hessian_);

    if (coef->is_variational()) {
      MFEM_VERIFY(!explicit_term || npts == 0 || explicit_gradients.Size() == gradients.Size(),
                  "SumCoefficient: a variational term with Scheme::Explicit needs the gradients "
                  "at the beginning of the time step.");
      term.df_dgrad = term_output(outputs.df_dgrad, this->term_df_dgrad_);
      if (!explicit_term) {
        term.d2f_dvar_dgrad = term_output(outputs.d2f_dvar_dgrad, this->term_d2f_dvar_dgrad_);
        term.d2f_dgrad2 = term_output(outputs.d2f_dgrad2, this->term_d2f_dgrad2_);
        term.d2f_dgrad2_isotropic =
            term_output(outputs.d2f_dgrad2_isotropic, this->term_d2f_dgrad2_isotropic_);
      }
      if (!coef->eval_variational_batch(npts, term_values, explicit_values, auxiliary_values,
                                        explicit_term ? explicit_gradients : gradients,
                                        explicit_gradients, auxiliary_gradients, coordinates,
                                        term)) {
        return false;
      }
    } else if (!coef->eval_batch(npts, term_values, explicit_values, auxiliary_values, coordinates,
                                 term.F, term.df_dvar, term.d2f_dvar2)) {
      return false;
    }
    accumulate(outputs.F, term.F);
    accumulate(outputs.df_dvar, term.df_dvar);
    accumulate(outputs.df_dgrad, term.df_dgrad);
    accumulate(outputs.d2f_dvar2, term.d2f_dvar2);
    accumulate(outputs.d2f_dvar_dgrad, term.d2f_dvar_dgrad);
    accumulate(outputs.d2f_dgrad2, term.d2f_dgrad2);
    accumulate(outputs.d2f_dgrad2_isotropic, term.d2f_dgrad2_isotropic);
  }
  return true;
}
