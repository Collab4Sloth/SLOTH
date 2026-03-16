/**
 * @file Coefficient.cpp
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

#include "Coefficients/Coefficient.hpp"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"

#ifdef SLOTH_USE_EXPRTK
#include "Coefficients/ExprTkSloth.hpp"  // NOLINT [no include the directory when naming exprtk include file]
#endif
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Constructs a Coefficient object from a user-provided function coefficient.
 *
 * This constructor allows creating a coefficient with a custom FunctionCoefficient
 * implementation.
 *
 * @param type Quantity associated with this coefficient.
 * @param scheme Discretization scheme used by this coefficient.
 * @param coef Shared pointer to the function-based coefficient.
 */
Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                         std::shared_ptr<FunctionCoefficient> coef)
    : SlothBaseCoefficient(type, scheme, std::move(coef)) {}

/**
 * @brief Constructs a Coefficient object with a constant value.
 *
 * @param type Quantity associated with this coefficient.
 * @param coef Constant coefficient value.
 */
Coefficient::Coefficient(GlossaryQuantity type, double coef)
    : SlothBaseCoefficient(type, std::move(coef)) {}
