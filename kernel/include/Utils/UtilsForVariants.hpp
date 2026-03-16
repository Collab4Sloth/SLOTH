/**
 * @file UtilsForVariants.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief manage std::Variants
 * @version 0.1
 * @date 2025-09-05
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
#include <memory>
#include <type_traits>
#include <variant>

#pragma once

template <typename T, typename Variant>
inline constexpr bool is_in_variant_v = false;

template <typename T, typename... Types>
inline constexpr bool is_in_variant_v<T, std::variant<Types...>> =
    (std::is_same_v<T, Types> || ...);

// Concatenation of a list of std::variant
template <class... Args>
struct concat_variants;

/**
 * @brief Concatenation of std::variants
 *
 * @tparam T
 * @tparam Args1
 * @tparam Args2
 * @tparam Remaining
 */
template <template <class...> class T, class... Args1, class... Args2,
          typename... RemainingVariants>
struct concat_variants<T<Args1...>, T<Args2...>, RemainingVariants...> {
  using type = typename concat_variants<T<Args1..., Args2...>, RemainingVariants...>::type;
};

/**
 * @brief Only one variant (useless but to avoid error?)
 *
 * @tparam T
 * @tparam Args1
 */
template <template <class...> class T, class... Args1>
struct concat_variants<T<Args1...>> {
  using type = T<Args1...>;
};

/**
 * @brief Alias for the template by providing variadic types and directly retrieving the
 * concatenated type.
 *
 * @tparam Args
 */
template <class... Args>
using concat_variant_type = typename concat_variants<Args...>::type;
