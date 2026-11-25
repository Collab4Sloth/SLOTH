/**
 * @file CoefficientTerm.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Term of a sum or of a product of FunctionCoefficient objects, with its own scheme
 * @version 0.1
 * @date 2026-11-23
 *
 * Copyright CEA (C) 2026
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

#include "Coefficients/FunctionCoefficient.hpp"
#include "Options/TimeOptions.hpp"

/**
 * @brief Term of a SumCoefficient or factor of a ProductCoefficient, with its own scheme
 *
 * - Scheme::Implicit (default): the term is evaluated with the values given to the sum or the
 *   product.
 * - Scheme::Explicit: the term is evaluated with the explicit values given to the sum or the
 *   product (values at the beginning of the time step). The scheme applies to its derivatives:
 *   its value and its derivatives are data, not differentiated in the second derivatives. The
 *   sum or the product must then be wrapped in a semi-implicit Coefficient.
 * - Scheme::SemiImplicit: same as Scheme::Implicit, the term reads the explicit values itself.
 */
struct CoefficientTerm {
  FunctionCoefficient* coefficient = nullptr;
  Scheme scheme = Scheme::Implicit;
};
