/**
 * @file SlothBaseCoefficient.hpp
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

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "FunctionCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

/**
 * @brief Base class definining a coefficient from an analytical formula
 *
 */
class SlothBaseCoefficient {
 private:
  GlossaryQuantity coefficient_type_;
  std::shared_ptr<FunctionCoefficient> coefficient_;

 public:
  SlothBaseCoefficient(GlossaryQuantity type, std::shared_ptr<FunctionCoefficient> coef);
  ~SlothBaseCoefficient();

  template <typename... Args>
  double compute(Args... args);

  template <typename... Args>
  double compute_gradient(const int i, Args... args);

  template <typename... Args>
  double compute_hessian(const int i, const int j, Args... args);

  GlossaryType get_type() const;
};
/**
 * @brief Construct a new SlothBaseCoefficient::SlothBaseCoefficient object
 *
 * @tparam Args
 * @param type
 * @param coef
 */
SlothBaseCoefficient::SlothBaseCoefficient(GlossaryQuantity type,
                                           std::shared_ptr<FunctionCoefficient> coef)
    : coefficient_type_(type), coefficient_(coef) {}

/**
 * @brief Compute  analytical expression
 *
 * @tparam Args
 * @param args
 * @return double
 */
template <typename... Args>
double SlothBaseCoefficient::compute(Args... args) {
  std::vector<double> values = {static_cast<double>(args)...};
  // if (values.size() != this->variable_names_.size()) {
  //   throw std::runtime_error("Number of variables not consistent with analytical expression");
  // }
  return this->coefficient_->eval_f(values);
}

/**
 * @brief Compute the i-th element of the gradient of analytical expression
 *
 * @tparam Args
 * @param id
 * @param args
 * @return double
 */
template <typename... Args>
double SlothBaseCoefficient::compute_gradient(const int id, Args... args) {
  std::vector<double> values = {static_cast<double>(args)...};
  // if (values.size() != this->variable_names_.size()) {
  //   throw std::runtime_error("Number of variables not consistent with analytical expression");
  // }
  return this->coefficient_->eval_gradient(id, values);
}

/**
 * @brief Compute the (i,j)-th element of the hessian of analytical expression
 *
 * @tparam Args
 * @param id
 * @param jd
 * @param args
 * @return double
 */
template <typename... Args>
double SlothBaseCoefficient::compute_hessian(const int id, const int jd, Args... args) {
  std::vector<double> values = {static_cast<double>(args)...};
  // if (values.size() != this->variable_names_.size()) {
  //   throw std::runtime_error("Number of variables not consistent with analytical expression");
  // }
  return this->coefficient_->eval_hessian(id, jd, values);
}

/**
 * @brief  Return the GlobalQuantity associated with the coefficient
 *
 * @return GlossaryType
 */
GlossaryType SlothBaseCoefficient::get_type() const { return this->coefficient_type_.type; }

/**
 * @brief Destroy the SlothBaseCoefficient::SlothBaseCoefficient object
 *
 */
SlothBaseCoefficient::~SlothBaseCoefficient() {}
