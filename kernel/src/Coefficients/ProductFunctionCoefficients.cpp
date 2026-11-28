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

#include <algorithm>
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

namespace {
constexpr const char* explicit_factor_error =
    "ProductCoefficient: a factor with Scheme::Explicit is evaluated with the values at the "
    "beginning of the time step: the product must be wrapped in a Coefficient with "
    "Scheme::SemiImplicit.";

/// @brief Values given to a factor: the explicit ones for an explicit factor
const std::span<const double>& factor_input(const bool explicit_factor,
                                            const std::span<const double>& input_vector,
                                            const std::span<const double>& exp_input_vector) {
  if (!explicit_factor) return input_vector;
  MFEM_VERIFY(!exp_input_vector.empty(), explicit_factor_error);
  return exp_input_vector;
}

/// @brief Work vector of the given size
mfem::Vector* sized(mfem::Vector& work, const int size) {
  work.SetSize(size);
  return &work;
}

/// @brief output = value (nothing if the output is not requested)
void copy_output(mfem::Vector* output, const mfem::Vector& value) {
  if (output != nullptr) *output = value;
}
}  // namespace

/**
 * @brief Build a FunctionCoefficient as the product of several FunctionCoefficients
 *
 * @param coeffs Factors (all implicit)
 */
ProductCoefficient::ProductCoefficient(std::initializer_list<FunctionCoefficient*> coeffs)
    : vect_coefficients_(coeffs), vect_schemes_(coeffs.size(), Scheme::Implicit) {
  this->check_terms();
}

/**
 * @brief Build a FunctionCoefficient as the product of several FunctionCoefficients, each one
 * with its own scheme (see ProductTerm)
 *
 * @param terms Factors of the product with their scheme
 */
ProductCoefficient::ProductCoefficient(std::initializer_list<ProductTerm> terms) {
  for (const auto& term : terms) {
    this->vect_coefficients_.push_back(term.coefficient);
    this->vect_schemes_.push_back(term.scheme);
  }
  this->check_terms();
}

/**
 * @brief Check the factors: not null, scheme Implicit, Explicit or SemiImplicit, no explicit
 * factor in a variational product
 */
void ProductCoefficient::check_terms() const {
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    MFEM_VERIFY(this->vect_coefficients_[k] != nullptr, "ProductCoefficient: null factor.");
    const Scheme scheme = this->vect_schemes_[k];
    MFEM_VERIFY(
        scheme == Scheme::Implicit || scheme == Scheme::Explicit || scheme == Scheme::SemiImplicit,
        "ProductCoefficient: the scheme of a factor must be Implicit, Explicit or SemiImplicit.");
  }
  MFEM_VERIFY(!(this->is_variational() && this->has_explicit_term()),
              "ProductCoefficient: a factor with Scheme::Explicit is not supported in a "
              "variational product (non-symmetric second derivatives with respect to the "
              "variables and their gradients).");
}

/**
 * @brief True if the k-th factor is explicit
 */
bool ProductCoefficient::is_explicit_term(const std::size_t k) const {
  return this->vect_schemes_[k] == Scheme::Explicit;
}

/**
 * @brief True if at least one factor is explicit
 */
bool ProductCoefficient::has_explicit_term() const {
  return std::ranges::any_of(this->vect_schemes_,
                             [](const Scheme s) { return s == Scheme::Explicit; });
}

/**
 * @brief  Compute the value of the coefficient as a product of the value of coefficients
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
    for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
      const auto& input = factor_input(this->is_explicit_term(k), input_vector, exp_input_vector);
      prod_F *=
          this->vect_coefficients_[k]->eval_f(input, exp_input_vector, auxiliary_vector, dimension);
    }

    return prod_F;
  };
  return func;
}

/**
 * @brief Gradient: sum over the factors k of grad f_k times the product of the other factors
 * (explicit factor: value and gradient at the beginning of the time step)
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
      const auto& input = factor_input(this->is_explicit_term(k), input_vector, exp_input_vector);
      f[k] =
          this->vect_coefficients_[k]->eval_f(input, exp_input_vector, auxiliary_vector, dimension);
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
      const auto& input = factor_input(this->is_explicit_term(k), input_vector, exp_input_vector);
      const double others = product_without(k);
      for (std::size_t i = 0; i < n; ++i) {
        gradient[i] +=
            this->vect_coefficients_[k]->eval_gradient(static_cast<int>(i), input, exp_input_vector,
                                                       auxiliary_vector, dimension) *
            others;
      }
    }
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian: derivatives of the gradient, the value and the gradient of the explicit factors
 * being held constant
 * @remark Hessian matrix stored in vector : H(i,j) = d(gradient[i])/du_j -> H(i*n+j)
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
      const auto& input = factor_input(this->is_explicit_term(k), input_vector, exp_input_vector);
      f[k] = coef->eval_f(input, exp_input_vector, auxiliary_vector, dimension);
      for (std::size_t i = 0; i < n; ++i) {
        g[k][i] = coef->eval_gradient(static_cast<int>(i), input, exp_input_vector,
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
      // Second derivatives of an implicit factor
      if (!this->is_explicit_term(k)) {
        const double others = product_without(k, k);
        for (std::size_t i = 0; i < n; ++i) {
          for (std::size_t j = 0; j < n; ++j) {
            hessian[i * n + j] += this->vect_coefficients_[k]->eval_hessian(
                                      static_cast<int>(i), static_cast<int>(j), input_vector,
                                      exp_input_vector, auxiliary_vector, dimension) *
                                  others;
          }
        }
      }
      // Gradient of factor k times the derivative of an implicit factor m
      for (std::size_t m = 0; m < nc; ++m) {
        if (m == k || this->is_explicit_term(m)) continue;
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
 * @brief Set the time step for all FunctionCoefficients
 *
 * @param dt Time step of the current step
 */
void ProductCoefficient::set_time_step(double dt) {
  for (auto* coef : vect_coefficients_) {
    coef->set_time_step(dt);
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
 * @brief True if every factor is a function of one variable only (then so is the product)
 */
bool ProductCoefficient::is_single_variable() const {
  if (vect_coefficients_.empty()) return false;
  for (const auto* coef : vect_coefficients_) {
    if (!coef->is_single_variable()) return false;
  }
  return true;
}

/**
 * @brief Number of variables of the evaluation by batch: the one of the factors, if all of them
 * are evaluated by batch with the same number of variables, 0 otherwise (point by point
 * evaluation).
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
 * @brief Evaluation by batch (see FunctionCoefficient::eval_batch): product rule applied factor
 * after factor, on the host or on the device. An explicit factor is evaluated with the explicit
 * values; its value and its gradient are held constant in the second derivatives:
 *   P' = P f,  dP' = dP f + P df
 *   implicit f: d2P' = d2P f + dP df^T + df dPi^T + P d2f
 *   explicit f: d2P' = d2P f + df dPi^T
 * with dPi the gradient of P, explicit factors held constant (dPi = dP without explicit factor).
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
  // Gradient with the explicit factors held constant: only for the second derivatives
  const bool separate_implicit_gradient = want_hessian && this->has_explicit_term();

  // Product of the first factors and its gradient: in the outputs if requested
  mfem::Vector& P = F != nullptr ? *F : this->product_F_;
  mfem::Vector& dP = gradient != nullptr ? *gradient : this->product_gradient_;
  mfem::Vector& dPi = separate_implicit_gradient ? this->product_implicit_gradient_ : dP;
  P.SetSize(npts);
  if (need_gradient) dP.SetSize(npts * n);

  bool first = true;
  for (std::size_t k = 0; k < this->vect_coefficients_.size(); ++k) {
    const auto* coef = this->vect_coefficients_[k];
    const bool explicit_factor = this->is_explicit_term(k);
    MFEM_VERIFY(!explicit_factor || npts == 0 || explicit_values.Size() > 0, explicit_factor_error);
    const mfem::Vector& factor_values = explicit_factor ? explicit_values : values;
    if (first) {
      if (!coef->eval_batch(npts, factor_values, explicit_values, auxiliary_values, coordinates, &P,
                            need_gradient ? &dP : nullptr, explicit_factor ? nullptr : hessian)) {
        return false;
      }
      if (want_hessian && explicit_factor) *hessian = 0.0;
      if (separate_implicit_gradient) {
        dPi.SetSize(npts * n);
        if (explicit_factor) {
          dPi = 0.0;
        } else {
          dPi = dP;
        }
      }
      first = false;
      continue;
    }
    this->term_F_.SetSize(npts);
    if (need_gradient) this->term_gradient_.SetSize(npts * n);
    const bool term_hessian = want_hessian && !explicit_factor;
    if (term_hessian) this->term_hessian_.SetSize(npts * n * n);
    if (!coef->eval_batch(npts, factor_values, explicit_values, auxiliary_values, coordinates,
                          &this->term_F_, need_gradient ? &this->term_gradient_ : nullptr,
                          term_hessian ? &this->term_hessian_ : nullptr)) {
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
    double* dpi = nullptr;
    if (separate_implicit_gradient) {
      dPi.UseDevice(true);
      dpi = dPi.ReadWrite();
    }
    const double* dpi_read = separate_implicit_gradient ? dpi : dp;
    double* d2p = nullptr;
    const double* d2f = nullptr;
    if (want_hessian) {
      hessian->UseDevice(true);
      d2p = hessian->ReadWrite();
    }
    if (term_hessian) {
      this->term_hessian_.UseDevice(true);
      d2f = this->term_hessian_.Read();
    }
    mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
      const double pq = p[q];
      const double fq = f[q];
      // Second derivatives first: they use the gradients and the value of the previous product
      if (d2p != nullptr) {
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            const int ij = (q * n + i) * n + j;
            if (explicit_factor) {
              d2p[ij] = d2p[ij] * fq + df[q * n + i] * dpi_read[q * n + j];
            } else {
              d2p[ij] = d2p[ij] * fq + dp[q * n + i] * df[q * n + j] +
                        df[q * n + i] * dpi_read[q * n + j] + pq * d2f[ij];
            }
          }
        }
      }
      if (dpi != nullptr) {
        for (int i = 0; i < n; ++i) {
          dpi[q * n + i] =
              explicit_factor ? dpi[q * n + i] * fq : dpi[q * n + i] * fq + pq * df[q * n + i];
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

/**
 * @brief True if at least one factor is variational: the product is then evaluated by
 * eval_variational_batch only
 */
bool ProductCoefficient::is_variational() const {
  return std::ranges::any_of(this->vect_coefficients_,
                             [](const FunctionCoefficient* c) { return c->is_variational(); });
}

/**
 * @brief Dimension of the gradients of the variational factors, 0 if there is none
 */
int ProductCoefficient::variational_dimension() const {
  for (const auto* coef : this->vect_coefficients_) {
    if (coef->is_variational()) return coef->variational_dimension();
  }
  return 0;
}

/**
 * @brief Number of variables of the variational factors, 0 if there is none
 */
int ProductCoefficient::variational_nb_variables() const {
  for (const auto* coef : this->vect_coefficients_) {
    if (coef->is_variational()) return coef->variational_nb_variables();
  }
  return 0;
}

/**
 * @brief True if a variational factor reads the gradients of its explicit variables
 */
bool ProductCoefficient::uses_explicit_gradients() const {
  return std::ranges::any_of(this->vect_coefficients_, [](const FunctionCoefficient* c) {
    return c->is_variational() && c->uses_explicit_gradients();
  });
}

/**
 * @brief True if a variational factor reads the gradients of its auxiliary variables
 */
bool ProductCoefficient::uses_auxiliary_gradients() const {
  return std::ranges::any_of(this->vect_coefficients_, [](const FunctionCoefficient* c) {
    return c->is_variational() && c->uses_auxiliary_gradients();
  });
}

/**
 * @brief Evaluation of a variational product (see FunctionCoefficient::eval_variational_batch):
 * product rule with respect to z = (u, grad u), factor after factor, on the host or on the
 * device. A non-variational factor, function of the variables of the product, has no derivative
 * with respect to the gradients. All the factors are implicit (see check_terms).
 *
 * @return false if the product is not variational or if a factor is not evaluated by batch
 */
bool ProductCoefficient::eval_variational_batch(
    const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
    const mfem::Vector& auxiliary_values, const mfem::Vector& gradients,
    const mfem::Vector& explicit_gradients, const mfem::Vector& auxiliary_gradients,
    const mfem::Vector& coordinates, const VariationalBatch& outputs) const {
  if (!this->is_variational()) return false;
  const int n = this->variational_nb_variables();
  const int d = this->variational_dimension();
  const int nd = n * d;
  for (const auto* coef : this->vect_coefficients_) {
    if (coef->is_variational()) {
      MFEM_VERIFY(coef->variational_nb_variables() == n && coef->variational_dimension() == d,
                  "ProductCoefficient: the variational factors must have the same variables and "
                  "the same dimension.");
    } else {
      MFEM_VERIFY(coef->batch_nb_variables() == n,
                  "ProductCoefficient: a non-variational factor of a variational product must be "
                  "a function of the " +
                      std::to_string(n) + " variables of the product, evaluated by batch.");
    }
  }
  const bool second = outputs.d2f_dvar2 != nullptr || outputs.d2f_dvar_dgrad != nullptr ||
                      outputs.d2f_dgrad2 != nullptr;

  // Product of the first factors and its derivatives (work arrays)
  mfem::Vector& P = this->product_F_;
  mfem::Vector& Pv = this->product_gradient_;
  mfem::Vector& Pg = this->product_df_dgrad_;
  mfem::Vector& Pvv = this->product_d2f_dvar2_;
  mfem::Vector& Pvg = this->product_d2f_dvar_dgrad_;
  mfem::Vector& Pgg = this->product_d2f_dgrad2_;
  P.SetSize(npts);
  Pv.SetSize(npts * n);
  Pg.SetSize(npts * nd);
  if (second) {
    Pvv.SetSize(npts * n * n);
    Pvg.SetSize(npts * n * nd);
    Pgg.SetSize(npts * nd * nd);
  }

  bool first = true;
  for (const auto* coef : this->vect_coefficients_) {
    const bool variational_factor = coef->is_variational();
    // Outputs of the factor: written in the product for the first one
    VariationalBatch term;
    term.F = first ? &P : sized(this->term_F_, npts);
    term.df_dvar = first ? &Pv : sized(this->term_gradient_, npts * n);
    if (second) term.d2f_dvar2 = first ? &Pvv : sized(this->term_hessian_, npts * n * n);
    if (variational_factor) {
      term.df_dgrad = first ? &Pg : sized(this->term_df_dgrad_, npts * nd);
      if (second) {
        term.d2f_dvar_dgrad = first ? &Pvg : sized(this->term_d2f_dvar_dgrad_, npts * n * nd);
        term.d2f_dgrad2 = first ? &Pgg : sized(this->term_d2f_dgrad2_, npts * nd * nd);
      }
      if (!coef->eval_variational_batch(npts, values, explicit_values, auxiliary_values, gradients,
                                        explicit_gradients, auxiliary_gradients, coordinates,
                                        term)) {
        return false;
      }
    } else if (!coef->eval_batch(npts, values, explicit_values, auxiliary_values, coordinates,
                                 term.F, term.df_dvar, term.d2f_dvar2)) {
      return false;
    }
    if (first) {
      // Non-variational first factor: no derivative with respect to the gradients
      if (!variational_factor) {
        Pg = 0.0;
        if (second) {
          Pvg = 0.0;
          Pgg = 0.0;
        }
      }
      first = false;
      continue;
    }

    // Data read and written where the loop runs (host or device)
    auto rw = [](mfem::Vector& v) {
      v.UseDevice(true);
      return v.ReadWrite();
    };
    auto rd = [](const mfem::Vector* v) -> const double* {
      if (v == nullptr) return nullptr;
      v->UseDevice(true);
      return v->Read();
    };
    double* p = rw(P);
    double* pv = rw(Pv);
    double* pg = rw(Pg);
    double* pvv = second ? rw(Pvv) : nullptr;
    double* pvg = second ? rw(Pvg) : nullptr;
    double* pgg = second ? rw(Pgg) : nullptr;
    const double* f = rd(term.F);
    const double* fv = rd(term.df_dvar);
    const double* fg = rd(term.df_dgrad);  // nullptr: non-variational factor
    const double* fvv = rd(term.d2f_dvar2);
    const double* fvg = rd(term.d2f_dvar_dgrad);
    const double* fgg = rd(term.d2f_dgrad2);
    mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
      const double pq = p[q];
      const double fq = f[q];
      const double* Gv = pv + q * n;
      const double* Gg = pg + q * nd;
      const double* gv = fv + q * n;
      const double* gg = fg != nullptr ? fg + q * nd : nullptr;
      // Second derivatives first: they use the gradients and the value of the previous product
      if (pvv != nullptr) {
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            const int ij = (q * n + i) * n + j;
            pvv[ij] = pvv[ij] * fq + Gv[i] * gv[j] + gv[i] * Gv[j] + pq * fvv[ij];
          }
          for (int a = 0; a < nd; ++a) {
            const int ia = (q * n + i) * nd + a;
            const double cross = gg != nullptr ? Gv[i] * gg[a] + pq * fvg[ia] : 0.;
            pvg[ia] = pvg[ia] * fq + cross + gv[i] * Gg[a];
          }
        }
        for (int a = 0; a < nd; ++a) {
          for (int b = 0; b < nd; ++b) {
            const int ab = (q * nd + a) * nd + b;
            const double own = gg != nullptr ? Gg[a] * gg[b] + gg[a] * Gg[b] + pq * fgg[ab] : 0.;
            pgg[ab] = pgg[ab] * fq + own;
          }
        }
      }
      for (int i = 0; i < n; ++i) pv[q * n + i] = Gv[i] * fq + pq * gv[i];
      for (int a = 0; a < nd; ++a) {
        pg[q * nd + a] = Gg[a] * fq + (gg != nullptr ? pq * gg[a] : 0.);
      }
      p[q] = pq * fq;
    });
  }

  copy_output(outputs.F, P);
  copy_output(outputs.df_dvar, Pv);
  copy_output(outputs.df_dgrad, Pg);
  if (second) {
    copy_output(outputs.d2f_dvar2, Pvv);
    copy_output(outputs.d2f_dvar_dgrad, Pvg);
    copy_output(outputs.d2f_dgrad2, Pgg);
  }
  return true;
}
