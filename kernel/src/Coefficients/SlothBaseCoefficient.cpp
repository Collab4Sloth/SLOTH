/**
 * @file SlothBaseCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class definining a coefficient from an analytical formula
 * @remark Depends on ExprTk library (MIT license)
 * @version 0.1
 * @date 2025-11-06
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

#include "Coefficients/SlothBaseCoefficient.hpp"

#include <iostream>
#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "Coefficients/ConstantCoefficient.hpp"
#include "Coefficients/FunctionCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Constructs a SlothBaseCoefficient with a user-defined coefficient.
 *
 * @param qty Quantity associated with the coefficient.
 * @param scheme Discretization scheme used by the coefficient.
 * @param coef Coefficient implementation (function-based).
 */
SlothBaseCoefficient::SlothBaseCoefficient(GlossaryQuantity qty, Scheme scheme,
                                           std::shared_ptr<FunctionCoefficient> coef)
    : coefficient_quantity_(qty), coefficient_(coef), scheme_(scheme) {}

/**
 * @brief Constructs a SlothBaseCoefficient with a constant value.
 *
 * @param qty Quantity associated with the coefficient.
 * @param coef Constant coefficient value.
 */
SlothBaseCoefficient::SlothBaseCoefficient(GlossaryQuantity qty, double coef)
    : coefficient_quantity_(qty),
      coefficient_(std::make_shared<ConstantCoefficient>(coef)),
      scheme_(Scheme::Constant) {}

/**
 * @brief Computes the value of the Coefficient.
 *
 * @return Value of the coefficient.
 */
double SlothBaseCoefficient::compute() {
  static const std::vector<double> values{};
  return this->coefficient_->eval_f(values);
}

/**
 * @brief Computes the value of the Coefficient.
 *
 * The coefficient is evaluated as a function of the current variable values.
 *
 * @param values Values of the current variables.
 * @param dimension Optional spatial dimension (if applicable).
 *
 * @return Value of the coefficient.
 */
double SlothBaseCoefficient::compute(const std::span<const double>& values,
                                     std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_f(values, dimension);
  } else {
    return this->coefficient_->eval_f(values);
  }
}

/**
 * @brief Computes the value of the Coefficient.
 *
 * The coefficient is evaluated as a function of the current variable values
 * and auxiliary variables.
 *
 * @param values Values of the current variables.
 * @param auxiliary_values Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (if applicable).
 *
 * @return Value of the coefficient.
 */
double SlothBaseCoefficient::compute(const std::span<const double>& values,
                                     const std::span<const double>& auxiliary_values,
                                     std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_f(values, auxiliary_values, dimension);
  } else {
    return this->coefficient_->eval_f(values, auxiliary_values);
  }
}

/**
 * @brief Computes the component id of the gradient.
 *
 * The component is evaluated as a function of the current variable values.
 *
 * @param id  index of the gradient matrix.
 * @param values Values of the current variables.
 * @param dimension Optional spatial dimension (if applicable).
 *
 * @return Value of the gradient component id.
 */
double SlothBaseCoefficient::compute_gradient(const int id, const std::span<const double>& values,
                                              std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_gradient(id, values, dimension);
  } else {
    return this->coefficient_->eval_gradient(id, values);
  }
}

/**
 * @brief Computes the component id of the gradient.
 *
 * The component is evaluated as a function of the current variable values
 * and auxiliary variables.
 *
 * @param id  index of the gradient matrix.
 * @param values Values of the current variables.
 * @param auxiliary_values Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (if applicable).
 *
 * @return Value of the gradient component id.
 */
double SlothBaseCoefficient::compute_gradient(const int id, const std::span<const double>& values,
                                              const std::span<const double>& auxiliary_values,
                                              std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_gradient(id, values, auxiliary_values, dimension);
  } else {
    return this->coefficient_->eval_gradient(id, values, auxiliary_values);
  }
}

/**
 * @brief Computes the (i,j) component of the Hessian matrix.
 *
 * The component is evaluated as a function of the current variable values.
 *
 * @param id Row index of the Hessian matrix.
 * @param jd Column index of the Hessian matrix.
 * @param values Values of the current variables.
 * @param dimension Optional spatial dimension (if applicable).
 *
 * @return Value of the Hessian component (id, jd).
 */
double SlothBaseCoefficient::compute_hessian(const int id, const int jd,
                                             const std::span<const double>& values,
                                             std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_hessian(id, jd, values, dimension);
  } else {
    return this->coefficient_->eval_hessian(id, jd, values);
  }
}

/**
 * @brief Computes the (i,j) component of the Hessian matrix.
 *
 * The component is evaluated as a function of the current variable values
 * and auxiliary variables.
 *
 * @param id Row index of the Hessian matrix.
 * @param jd Column index of the Hessian matrix.
 * @param values Values of the current variables.
 * @param auxiliary_values Values of the auxiliary variables.
 * @param dimension Optional spatial dimension (if applicable).
 *
 * @return Value of the Hessian component (id, jd).
 */
double SlothBaseCoefficient::compute_hessian(const int id, const int jd,
                                             const std::span<const double>& values,
                                             const std::span<const double>& auxiliary_values,
                                             std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_hessian(id, jd, values, auxiliary_values, dimension);
  } else {
    return this->coefficient_->eval_hessian(id, jd, values, auxiliary_values);
  }
}

/**
 * @brief  Return the type of the quantity associated with the coefficient
 *
 * @return the type of the quantity.
 */
GlossaryType SlothBaseCoefficient::get_type() const { return this->coefficient_quantity_.type; }
/**
 * @brief Returns the quantity associated with the coefficient.
 *
 * @return the quantity.
 */
const GlossaryQuantity SlothBaseCoefficient::get_quantity() { return this->coefficient_quantity_; }

/**
 * @brief Returns the ID of the quantity associated with the coefficient.
 *
 * @return ID of the associated quantity.
 */
unsigned int SlothBaseCoefficient::get_id() const { return this->coefficient_quantity_.id; }

/**
 * @brief Indicates whether the coefficient is discretized with an implicit scheme.
 *
 * @return true if the associated scheme is Scheme::Implicit,
 *         false otherwise.
 */
bool SlothBaseCoefficient::is_implicit() const { return this->scheme_ == Scheme::Implicit; }

/**
 * @brief Indicates whether the coefficient is discretized with an explicit scheme.
 *
 * @return true if the associated scheme is Scheme::Explicit,
 *         false otherwise.
 */
bool SlothBaseCoefficient::is_explicit() const { return this->scheme_ == Scheme::Explicit; }

/**
 * @brief Indicates whether the coefficient is discretized with a semi-implicit scheme.
 *
 * @return true if the associated scheme is Scheme::SemiImplicit,
 *         false otherwise.
 */
bool SlothBaseCoefficient::is_semi_implicit() const {
  return this->scheme_ == Scheme::SemiImplicit;
}

/**
 * @brief Indicates whether the coefficient is scalar.
 *
 * @return true if the associated scheme is Scheme::Constant,
 *         false otherwise.
 */
bool SlothBaseCoefficient::is_scalar() const { return this->scheme_ == Scheme::Constant; }
