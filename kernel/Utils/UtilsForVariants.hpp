/**
 * @file UtilsForVariants.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Usefull methods to manage std::Variants
 * @version 0.1
 * @date 2024-12-13
 *
 * Copyright CEA (c) 2024
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
