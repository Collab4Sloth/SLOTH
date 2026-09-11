/**
 * @file CouplingFactory.hpp
 * @brief Helper building Coupling<PB, PB, ..., PB> (N times) from
 *        a std::vector<PB>, with N given explicitly as a template parameter.
 * @version 0.1
 *
 * @anchor couplings
 *
 * @copyright CEA (C) 2026
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
#include <cassert>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "Couplings/Coupling.hpp"

/**
 * @brief Maps an index to PB, whatever the index -> used to repeat PB N times
 *        in a template parameter pack (Coupling<RepeatedProblemType<PB,0>,
 * RepeatedProblemType<PB,1>, ...> = Coupling<PB, PB, ...>).
 *
 * @tparam PB Problem type to repeat.
 */
template <class PB, std::size_t>
using RepeatedProblemType = PB;

/**
 * @brief Implementation detail of HomogeneousCoupling: expands the index sequence
 *        to build the Coupling<PB, PB, ..., PB> return type and forward each
 *        vector element by index.
 */
template <class PB, std::size_t N, std::size_t... I>
Coupling<RepeatedProblemType<PB, I>...> HomogeneousCoupling(const std::string& name,
                                                            std::vector<PB> vect,
                                                            std::index_sequence<I...>) {
  return Coupling<RepeatedProblemType<PB, I>...>(name, std::move(vect[I])...);
}

/**
 * @brief Build a Coupling<PB, PB, ..., PB> (N times) from a std::vector<PB>
 *        whose size is only known at runtime.
 *
 *
 * @tparam N Number of problems in `vect`. Must match `vect.size()` exactly
 *           (checked with assert in debug builds).
 * @tparam PB Problem type shared by every element of `vect`.
 * @param name Name of the coupling.
 * @param vect Problems to couple, all of type PB, size exactly N.
 * @return Coupling<PB, PB, ..., PB> (N times), ready to pass to
 *         TimeDiscretization.
 */
template <std::size_t N, class PB>
auto setCoupling(const std::string& name, std::vector<PB> vect) {
  MFEM_VERIFY(vect.size() == N, "setCoupling<N>: vect must contain exactly N elements");
  return HomogeneousCoupling<PB, N>(name, std::move(vect), std::make_index_sequence<N>{});
}