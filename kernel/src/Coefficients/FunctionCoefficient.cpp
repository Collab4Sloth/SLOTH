/**
 * @file FunctionCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class defining analytical expression, first and second derivatives of a Coefficient
 * @version 0.1
 * @date 2025-11-06

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

#include "Coefficients/FunctionCoefficient.hpp"

#include <any>
#include <cmath>
#include <concepts>
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

namespace {
/**
 * @brief Dimension passed to F(), GradientF() and HessianF(): the one given if it is
 * non-negative, the number of input values otherwise (e.g. GradientEnergy then uses all the
 * components of the gradient given as input).
 */
unsigned int resolve_dimension(std::optional<int> dimension,
                               const std::span<const double>& input_vector) {
  if (dimension.has_value() && *dimension >= 0) {
    return static_cast<unsigned int>(*dimension);
  }
  return static_cast<unsigned int>(input_vector.size());
}
}  // namespace

/**
 * @brief Construct a new FunctionCoefficient::FunctionCoefficient object
 *
 */
FunctionCoefficient::FunctionCoefficient() : time_(0.0) {}

/**
 * @brief Evaluates the coefficient.
 *
 * The coefficient is computed from the current variables.
 *
 * @param input_vector Values of the current variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the coefficient.
 */
double FunctionCoefficient::eval_f(const std::span<const double>& input_vector,
                                   std::optional<int> dimension) {
  static const std::span<const double> empty_aux{};
  static const std::span<const double> empty_exp{};
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  return F()(input_vector, empty_exp, empty_aux, dim);
}

/**
 * @brief Evaluates the coefficient.
 *
 * @remark Work for Implicit/Explict coefficients with auxiliary variables and Semi-implicit
 * coefficients without auxiliary variables
 *
 * The coefficient is computed from the current and auxiliary (explicit) variables.
 *
 * @param input_vector Values of the current variables.
 * @param auxiliary_vector Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the coefficient.
 */
double FunctionCoefficient::eval_f(const std::span<const double>& input_vector,
                                   const std::span<const double>& auxiliary_vector,
                                   std::optional<int> dimension) {
  static const std::span<const double> empty_exp{};
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  return F()(input_vector, empty_exp, auxiliary_vector, dim);
}

/**
 * @brief Evaluates the coefficient.
 *
 * The coefficient is computed from the current and auxiliary variables.
 *
 * @param input_vector Values of the current variables.
 * @param exp_input_vector Values of the previous variables.
 * @param auxiliary_vector Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the coefficient.
 */
double FunctionCoefficient::eval_f(const std::span<const double>& input_vector,
                                   const std::span<const double>& exp_input_vector,
                                   const std::span<const double>& auxiliary_vector,
                                   std::optional<int> dimension) {
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  return F()(input_vector, exp_input_vector, auxiliary_vector, dim);
}

/**
 * @brief Evaluates the component i of the gradient.
 *
 * The gradient is computed from the current variables.
 *
 * @param i  index of the gradient .
 * @param input_vector Values of the current variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the gradient component i.
 */
double FunctionCoefficient::eval_gradient(const int i, const std::span<const double>& input_vector,
                                          std::optional<int> dimension) {
  static const std::span<const double> empty_aux{};
  static const std::span<const double> empty_exp{};
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  return GradientF()(input_vector, empty_exp, empty_aux, dim)[i];
}

/**
 * @brief Evaluates the component i of the gradient.
 *
 * @remark Work for Implicit/Explict coefficients with auxiliary variables and Semi-implicit
 * coefficients without auxiliary variables
 *
 * The gradient is computed from the current and auxiliary (explicit) variables.
 *
 * @param i  index of the gradient .
 * @param input_vector Values of the current variables.
 * @param auxiliary_vector Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the gradient component i.
 */
double FunctionCoefficient::eval_gradient(const int i, const std::span<const double>& input_vector,
                                          const std::span<const double>& auxiliary_vector,
                                          std::optional<int> dimension) {
  static const std::span<const double> empty_exp{};
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  return GradientF()(input_vector, empty_exp, auxiliary_vector, dim)[i];
}

/**
 * @brief Evaluates the component i of the gradient.
 *
 * The gradient is computed from the current and auxiliary variables.
 *
 * @param i  index of the gradient .
 * @param input_vector Values of the current variables.
 * @param exp_input_vector Values of the previous variables.
 * @param auxiliary_vector Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the gradient component i.
 */
double FunctionCoefficient::eval_gradient(const int i, const std::span<const double>& input_vector,
                                          const std::span<const double>& exp_input_vector,
                                          const std::span<const double>& auxiliary_vector,
                                          std::optional<int> dimension) {
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  return GradientF()(input_vector, exp_input_vector, auxiliary_vector, dim)[i];
}

/**
 * @brief Evaluates the (i,j) component of the Hessian matrix.
 *
 * The Hessian is computed from the current variables and returned
 * as a flattened array. The component is accessed using row-major ordering.
 *
 * @param i Row index of the Hessian matrix.
 * @param j Column index of the Hessian matrix.
 * @param input_vector Values of the current variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the Hessian component (i, j).
 */
double FunctionCoefficient::eval_hessian(const int i, const int j,
                                         const std::span<const double>& input_vector,
                                         std::optional<int> dimension) {
  static const std::span<const double> empty_aux{};
  static const std::span<const double> empty_exp{};
  const unsigned int dim = resolve_dimension(dimension, input_vector);
  const int size = input_vector.size();
  return HessianF()(input_vector, empty_exp, empty_aux, dim)[i * size + j];
}

/**
 * @brief Evaluates the (i,j) component of the Hessian matrix.
 *
 * @remark Work for Implicit/Explict coefficients with auxiliary variables and Semi-implicit
 * coefficients without auxiliary variables
 *
 * The Hessian is computed from the current and auxiliary (explicit) variables and returned
 * as a flattened array. The component is accessed using row-major ordering.
 *
 * @param i Row index of the Hessian matrix.
 * @param j Column index of the Hessian matrix.
 * @param input_vector Values of the current variables.
 * @param auxiliary_vector Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the Hessian component (i, j).
 */
double FunctionCoefficient::eval_hessian(const int i, const int j,
                                         const std::span<const double>& input_vector,
                                         const std::span<const double>& auxiliary_vector,
                                         std::optional<int> dimension) {
  static const std::span<const double> empty_exp{};
  const int size = input_vector.size();
  const unsigned int dim = resolve_dimension(dimension, input_vector);

  return HessianF()(input_vector, empty_exp, auxiliary_vector, dim)[i * size + j];
}

/**
 * @brief Evaluates the (i,j) component of the Hessian matrix.
 *
 * The Hessian is computed from the current and auxiliary variables and returned
 * as a flattened array. The component is accessed using row-major ordering.
 *
 * @param i Row index of the Hessian matrix.
 * @param j Column index of the Hessian matrix.
 * @param input_vector Values of the current variables.
 * @param exp_input_vector Values of the previous variables.
 * @param auxiliary_vector Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (number of input values if not given or negative).
 *
 * @return Value of the Hessian component (i, j).
 */
double FunctionCoefficient::eval_hessian(const int i, const int j,
                                         const std::span<const double>& input_vector,
                                         const std::span<const double>& exp_input_vector,
                                         const std::span<const double>& auxiliary_vector,
                                         std::optional<int> dimension) {
  const int size = input_vector.size();
  const unsigned int dim = resolve_dimension(dimension, input_vector);

  return HessianF()(input_vector, exp_input_vector, auxiliary_vector, dim)[i * size + j];
}
/**
 * @brief Update the time associated with the coefficient
 *
 */
void FunctionCoefficient::set_time(double time) { this->time_ = time; }

/**
 * @brief Return the time associated with the coefficient
 *
 */
double FunctionCoefficient::get_time() const { return this->time_; }

/**
 * @brief Set the physical coordinates of the point where the coefficient is evaluated
 *
 * @remark Only used by coefficients depending on the spatial coordinates
 *         (see uses_coordinates()). No reallocation after the first call.
 *
 * @param x Physical coordinates
 */
void FunctionCoefficient::set_coordinates(std::span<const double> x) {
  this->coordinates_.assign(x.begin(), x.end());
}