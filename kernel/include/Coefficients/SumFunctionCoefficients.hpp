/**
 * @file SumFunctionCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class defining analytical expression, first and second derivatives of a sum of
 * FunctionCoefficient objects
 * @version 0.1
 * @date 2025-02-03
 *
 * Copyright CEA (C) 2025
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

#include "Coefficients/CoefficientTerm.hpp"
#include "Coefficients/FunctionCoefficient.hpp"

/// Term of a SumCoefficient with its own scheme (see CoefficientTerm)
using SumTerm = CoefficientTerm;

/**
 * @brief Define a summation of FunctionCoefficient objects as a FunctionCoefficient
 * @remark All the terms share the same variables. Each term has its own scheme (see SumTerm).
 * Variational and non-variational terms can be summed: the sum is then variational, the
 * non-variational terms only contribute to f, df/du and d2f/du2.
 *
 */
class SumCoefficient : public FunctionCoefficient {
 private:
  std::vector<FunctionCoefficient*> vect_coefficients_;
  std::vector<Scheme> vect_schemes_;  ///< scheme of each term
  // Work arrays of eval_batch and eval_variational_batch
  mutable mfem::Vector term_F_, term_gradient_, term_hessian_;  ///< outputs of one term
  mutable mfem::Vector term_df_dgrad_, term_d2f_dvar_dgrad_, term_d2f_dgrad2_;
  mutable mfem::Vector term_d2f_dgrad2_isotropic_;

  bool is_explicit_term(std::size_t k) const;
  void check_terms() const;

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
  explicit SumCoefficient(std::initializer_list<FunctionCoefficient*> coeffs);
  explicit SumCoefficient(std::initializer_list<SumTerm> terms);
  SumCoefficient() = default;

  virtual ~SumCoefficient() = default;

  void set_time(double time) override;
  bool uses_coordinates() const override;
  bool is_single_variable() const override;
  void set_coordinates(std::span<const double> x) override;

  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const override;
  int batch_nb_variables() const override;

  bool is_variational() const override;
  int variational_dimension() const override;
  int variational_nb_variables() const override;
  bool uses_explicit_gradients() const override;
  bool uses_auxiliary_gradients() const override;
  bool is_isotropic() const override;
  bool eval_variational_batch(const int npts, const mfem::Vector& values,
                              const mfem::Vector& explicit_values,
                              const mfem::Vector& auxiliary_values, const mfem::Vector& gradients,
                              const mfem::Vector& explicit_gradients,
                              const mfem::Vector& auxiliary_gradients,
                              const mfem::Vector& coordinates,
                              const VariationalBatch& outputs) const override;
};
