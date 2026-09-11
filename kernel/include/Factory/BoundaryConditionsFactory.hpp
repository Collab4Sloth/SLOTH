/**
 * @file BoundaryConditionsFactory.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief BoundaryConditions factory to manage boundary conditions
 * @version 2.1
 * @date 2026-09-05
 *
 * @anchor bcs
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

#include <cstddef>
#include <vector>

#include "BCs/BoundaryConditions.hpp"
/**
 * @brief Build N BoundaryConditions<T, DIM>, one per element of `spatials`,
 *        all sharing the same `boundaries...` list.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @tparam BoundaryArgs Deduced argument types forwarded to every
 *                      BoundaryConditions<T, DIM> constructor call.
 * @param N Number of BoundaryConditions to build. Must match spatials.size().
 * @param spatials One SpatialDiscretization<T, DIM>* per grain (SPAS).
 * @param boundaries Forwarded as-is to each BoundaryConditions constructor.
 * @return std::vector<BoundaryConditions<T, DIM>> of size N.
 */
template <class T, int DIM, class... BoundaryArgs>
std::vector<BoundaryConditions<T, DIM>> setBoundaryConditions(
    std::size_t N, const std::vector<SpatialDiscretization<T, DIM>*>& spatials,
    const BoundaryArgs&... boundaries) {
  MFEM_VERIFY(N >= 1, "setBoundaryConditions: N must be >= 1");
  MFEM_VERIFY(spatials.size() == N,
              "setBoundaryConditions: spatials must contain exactly N elements");

  std::vector<BoundaryConditions<T, DIM>> bcs;
  bcs.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    bcs.emplace_back(spatials[i], boundaries...);
  }
  return bcs;
}