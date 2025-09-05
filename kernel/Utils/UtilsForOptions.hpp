/**
 * @file UtilsForOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Usefull methods for Options
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
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#pragma once

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
namespace PhaseFieldPrivate {
template <typename T>
struct mmap : private std::vector<std::pair<const char* const, T>> {
  using mpair = std::pair<const char* const, T>;
  mmap(const std::initializer_list<mpair>&);
  T find(const char* const, const std::string&);
};

/**
 * @brief Construct a new mmap<E Type>::mmap object
 *
 * @tparam T
 * @param values
 */
template <typename T>
mmap<T>::mmap(const std::initializer_list<typename mmap::mpair>& values)
    : std::vector<typename mmap::mpair>(values) {}

/**
 * @brief
 *
 * @tparam T
 * @param n
 * @param v
 * @return T
 */
template <typename T>
T mmap<T>::find(const char* const n, const std::string& v) {
  const auto pe = this->end();
  const auto p = std::find_if(this->begin(), pe, [&v](const mpair& e) {
    // Convert in uppercase for sake of generality
    std::string s1 = v;
    std::string s2 = e.first;
    transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
    transform(s2.begin(), s2.end(), s2.begin(), ::toupper);

    return s1 == s2;
  });
  if (p == pe) {
    std::string msg =
        "EnumNotFound::EnumNotFound : "
        "string '" +
        std::string(n) +
        "' "
        "is not a valid value for '" +
        v + "' keyword. See documentation";
    mfem::mfem_error(msg.c_str());
  }
  return p->second;
}
}  // namespace PhaseFieldPrivate
