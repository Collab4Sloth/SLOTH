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
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "ConstantCoefficient.hpp"
#include "FunctionCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

/**
 * @brief Base class definining a coefficient from an analytical formula
 *
 */
class SlothBaseCoefficient {
 private:
  GlossaryQuantity coefficient_quantity_;
  std::shared_ptr<FunctionCoefficient> coefficient_;

  Scheme scheme_;

 public:
  SlothBaseCoefficient(GlossaryQuantity type, Scheme scheme,
                       std::shared_ptr<FunctionCoefficient> coef);
  SlothBaseCoefficient(GlossaryQuantity type, double coef);

  ~SlothBaseCoefficient();

  double compute();
  double compute(const std::vector<double>& values, std::optional<int> dimension = std::nullopt);
  double compute(const std::vector<double>& values, const std::vector<double>& auxiliary_values,
                 std::optional<int> dimension = std::nullopt);

  double compute_gradient(const int i, const std::vector<double>& values,
                          std::optional<int> dimension = std::nullopt);
  double compute_gradient(const int i, const std::vector<double>& values,
                          const std::vector<double>& auxiliary_values,
                          std::optional<int> dimension = std::nullopt);

  double compute_hessian(const int i, const int j, const std::vector<double>& values,
                         std::optional<int> dimension = std::nullopt);
  double compute_hessian(const int i, const int j, const std::vector<double>& values,
                         const std::vector<double>& auxiliary_values,
                         std::optional<int> dimension = std::nullopt);

  GlossaryType get_type() const;
  const GlossaryQuantity get_quantity();
  unsigned int get_id();
  bool is_implicit() const;
  bool is_explicit() const;
  bool is_semi_implicit() const;
  bool is_scalar() const;
};
/**
 * @brief Construct a new SlothBaseCoefficient::SlothBaseCoefficient object
 *
 * @tparam Args
 * @param type
 * @param coef
 */
SlothBaseCoefficient::SlothBaseCoefficient(GlossaryQuantity qty, Scheme scheme,
                                           std::shared_ptr<FunctionCoefficient> coef)
    : coefficient_quantity_(qty), coefficient_(coef), scheme_(scheme) {}

SlothBaseCoefficient::SlothBaseCoefficient(GlossaryQuantity qty, double coef)
    : coefficient_quantity_(qty),
      coefficient_(std::make_shared<ConstantCoefficient>(coef)),
      scheme_(Scheme::Constant) {}

/**
 * @brief Compute  analytical expression
 *
 * @tparam Args
 * @param args
 * @return double
 */
double SlothBaseCoefficient::compute() {
  static const std::vector<double> values{};
  return this->coefficient_->eval_f(values);
}
double SlothBaseCoefficient::compute(const std::vector<double>& values,
                                     std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_f(values, dimension);
  } else {
    return this->coefficient_->eval_f(values);
  }
}

double SlothBaseCoefficient::compute(const std::vector<double>& values,
                                     const std::vector<double>& auxiliary_values,
                                     std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_f(values, auxiliary_values, dimension);
  } else {
    return this->coefficient_->eval_f(values, auxiliary_values);
  }
}

/**
 * @brief Compute the i-th element of the gradient of analytical expression
 *
 * @tparam Args
 * @param id
 * @param args
 * @return double
 */

double SlothBaseCoefficient::compute_gradient(const int id, const std::vector<double>& values,
                                              std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_gradient(id, values, dimension);
  } else {
    return this->coefficient_->eval_gradient(id, values);
  }
}
double SlothBaseCoefficient::compute_gradient(const int id, const std::vector<double>& values,
                                              const std::vector<double>& auxiliary_values,
                                              std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_gradient(id, values, auxiliary_values, dimension);
  } else {
    return this->coefficient_->eval_gradient(id, values, auxiliary_values);
  }
}

/**
 * @brief Compute the (i,j)-th element of the hessian of analytical expression
 *
 * @param id
 * @param jd
 * @param args
 * @return double
 */

double SlothBaseCoefficient::compute_hessian(const int id, const int jd,
                                             const std::vector<double>& values,
                                             std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_hessian(id, jd, values, dimension);
  } else {
    return this->coefficient_->eval_hessian(id, jd, values);
  }
}

double SlothBaseCoefficient::compute_hessian(const int id, const int jd,
                                             const std::vector<double>& values,
                                             const std::vector<double>& auxiliary_values,
                                             std::optional<int> dimension) {
  if (dimension.has_value()) {
    return this->coefficient_->eval_hessian(id, jd, values, auxiliary_values, dimension);
  } else {
    return this->coefficient_->eval_hessian(id, jd, values, auxiliary_values);
  }
}

/**
 * @brief  Return the GlobalQuantity associated with the coefficient
 *
 * @return GlossaryType
 */
inline GlossaryType SlothBaseCoefficient::get_type() const {
  return this->coefficient_quantity_.type;
}

inline const GlossaryQuantity SlothBaseCoefficient::get_quantity() {
  return this->coefficient_quantity_;
}

inline unsigned int SlothBaseCoefficient::get_id() { return this->coefficient_quantity_.id; }

inline bool SlothBaseCoefficient::is_implicit() const { return this->scheme_ == Scheme::Implicit; }

inline bool SlothBaseCoefficient::is_explicit() const { return this->scheme_ == Scheme::Explicit; }

inline bool SlothBaseCoefficient::is_semi_implicit() const {
  return this->scheme_ == Scheme::SemiImplicit;
}
inline bool SlothBaseCoefficient::is_scalar() const { return this->scheme_ == Scheme::Constant; }

/**
 * @brief Destroy the SlothBaseCoefficient::SlothBaseCoefficient object
 *
 */
SlothBaseCoefficient::~SlothBaseCoefficient() {}
