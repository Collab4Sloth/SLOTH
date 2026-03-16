/**
 * @file Coefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Coefficients class used to build and manage a collection of Coefficient
 * @version 0.1
 * @date 2025-11-26
 *
 * @anchor coefficients
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

#include "Coefficients/Coefficients.hpp"

#include <any>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "Glossary/Glossary.hpp"

/**
 * @brief Adds a new coefficient to the container.
 *
 * @param coef Coefficient to be added.
 */
void Coefficients::add(Coefficient coef) { this->vect_coefficients_.push_back(std::move(coef)); }

/**
 * @brief Returns the vector of coefficients.
 *
 * @return Vector containing all coefficients.
 */
std::vector<Coefficient> Coefficients::getCoefficients() const { return this->vect_coefficients_; }

/**
 * @brief Returns the number of coefficients.
 *
 * @return Number of stored coefficients.
 */
size_t Coefficients::size() noexcept { return this->vect_coefficients_.size(); }

/**
 * @brief Returns the i-th coefficient.
 *
 * @param i Index of the coefficient.
 * @return Reference to the i-th coefficient.
 *
 * @throws std::out_of_range if i is out of bounds.
 */
Coefficient Coefficients::operator[](size_t i) {
  if (i >= vect_coefficients_.size()) throw std::out_of_range("Index out of range");
  return this->vect_coefficients_[i];
}

/**
 * @brief Returns a copy of the glossary types associated with the coefficients.
 *
 * @return Vector of glossary types.
 */
std::vector<GlossaryType> Coefficients::get_types() const { return this->vect_coefficient_types_; }
