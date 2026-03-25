/**
 * @file Coefficients.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Coefficients class used to build and manage a collection of Coefficient
 * @version 0.1
 * @date 2025-11-26
 *
 * @anchor coefficients
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

#pragma once
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

template <class T, class Coefficient>
concept CoeffVar = std::same_as<std::remove_cvref_t<T>, Coefficient>;

/**
 * @brief Class used to manage a collection of Coefficient
 *
 */
class Coefficients {
 private:
  std::vector<Coefficient> vect_coefficients_;
  std::vector<GlossaryType> vect_coefficient_types_;

 public:
  template <CoeffVar<Coefficient>... Args>
    requires((sizeof...(Args) > 0))
  explicit Coefficients(const Args&... args);

  Coefficients() = default;
  virtual ~Coefficients() = default;

  void add(Coefficient coef);
  std::vector<Coefficient> getCoefficients() const;
  size_t size() noexcept;
  Coefficient operator[](size_t i);

  std::vector<GlossaryType> get_types() const;
};

/**
 * @brief Constructs a Coefficients object from a list of Coefficient objects.
 *
 * @tparam Args Variadic template parameters constrained to Coefficient-compatible types.
 * @param args List of coefficients used to initialize the container.
 */
template <CoeffVar<Coefficient>... Args>
  requires((sizeof...(Args) > 0))
Coefficients::Coefficients(const Args&... args)
    : vect_coefficients_{args...}, vect_coefficient_types_{args.get_type()...} {}
