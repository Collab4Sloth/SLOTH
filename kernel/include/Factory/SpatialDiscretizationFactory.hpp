/**
 * @file SpatialDiscretizationFactory.hpp
 * @brief Helper building N SpatialDiscretization<T, DIM> objects sharing the
 *        same underlying mesh, from any of SpatialDiscretization's
 *        constructors.
 * @version 0.1
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
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "Spatial/Spatial.hpp"

/**
 * @brief Build N SpatialDiscretization<T, DIM> objects sharing the same mesh,
 *        returned directly as a `SPAS`-compatible `std::vector<SpatialDiscretization<T, DIM>*>`:
 *        the first one is built by forwarding `args...` verbatim to whichever
 *        `SpatialDiscretization<T, DIM>` constructor they match - any of the
 *        four (mesh-file, inline-mesh, periodic inline-mesh, or
 *        existing-mesh overloads) - and the remaining N-1 share that first
 *        one's mesh, exactly like building `spatial_phi1`, `spatial_phi2`,
 *        ... from `spatial_phi0.get_mesh()` by hand.
 *
 * @details A single templated, perfectly-forwarding function is used instead
 *          of one factory per constructor, so this stays in sync with
 *          `SpatialDiscretization`'s constructors automatically. `args...`
 *          must be exactly the argument list of one of the existing
 *          constructors (everything except `N` and `shared_is_periodic`).
 *
 *          `T` and `DIM` cannot be deduced from `args...`, so they must
 *          always be given explicitly at the call site:
 *          `setSpatialDiscretization<mfem::H1_FECollection, 2>(30, false, ...)`.
 *
 *          In all constructors, `fe_order` is always the 2nd parameter;
 *          `std::get<1>` on the forwarded pack recovers it generically to
 *          build the shared-mesh elements via the existing-mesh constructor.
 *
 * @warning OWNERSHIP: each element is allocated with `new` and returned as a
 *          raw pointer - the caller now owns them and MUST eventually
 *          release them with `deleteSpatialDiscretization` (below), not with
 *          a manual loop of `delete`. `result[0]` owns the mesh (its
 *          destructor deletes it); `result[1..N-1]` reference it without
 *          owning it. Deleting `result[0]` before `result[1..N-1]` would
 *          free the mesh while they still use it (use-after-free) -
 *          `deleteSpatialDiscretization` deletes in the correct (reverse)
 *          order for you.
 *
 * @tparam T Finite element collection type. Must be given explicitly.
 * @tparam DIM Spatial dimension. Must be given explicitly.
 * @tparam Args Deduced argument types, forwarded as-is to the
 *             `SpatialDiscretization<T, DIM>` constructor used to build the
 *             first element.
 * @param N Total number of SpatialDiscretization objects to build (must be
 *          >= 1). The first is built from `args...`, the remaining N-1
 *          share its mesh.
 * @param shared_is_periodic Value passed as `is_periodic_mesh` to the N-1
 *                           mesh-sharing elements. Unused when N == 1.
 * @param args Arguments forwarded verbatim to the `SpatialDiscretization<T, DIM>`
 *            constructor used for the first element.
 * @return std::vector<SpatialDiscretization<T, DIM>*> of size N, owned by
 *         the caller - release with `deleteSpatialDiscretization`.
 */
template <class T, int DIM, class... Args>
std::vector<SpatialDiscretization<T, DIM>*> setSpatialDiscretization(std::size_t N,
                                                                     bool shared_is_periodic,
                                                                     Args&&... args) {
  static_assert(sizeof...(Args) >= 2,
                "setSpatialDiscretization: args must match one of "
                "SpatialDiscretization's constructors, e.g. (mesh_type, fe_order, ...) "
                "or (existing_mesh, fe_order, ...).");
  MFEM_VERIFY(N >= 1, "setSpatialDiscretization: N must be >= 1");

  // Phase 1: build safely. Each element is owned by a std::unique_ptr in
  // `storage`, so std::vector reallocation (if capacity were ever exceeded)
  // would only move unique_ptr's themselves, never touch the pointee -
  // sidesteps SpatialDiscretization having no move constructor entirely.
  std::vector<std::unique_ptr<SpatialDiscretization<T, DIM>>> storage;
  storage.reserve(N);

  // fe_order is always the 2nd parameter of every SpatialDiscretization
  // constructor - grab it before args... is consumed by forwarding.
  const int fe_order = std::get<1>(std::forward_as_tuple(args...));

  // storage[0]: built by forwarding args... to whichever constructor matches.
  storage.push_back(std::make_unique<SpatialDiscretization<T, DIM>>(std::forward<Args>(args)...));

  // storage[1..N-1]: share storage[0]'s mesh, own their own FE space.
  mfem::ParMesh* shared_mesh = storage.front()->get_mesh();
  for (std::size_t i = 1; i < N; ++i) {
    storage.push_back(
        std::make_unique<SpatialDiscretization<T, DIM>>(shared_mesh, fe_order, shared_is_periodic));
  }

  // Phase 2: transfer ownership to the caller as raw pointers (SPAS). Each
  // unique_ptr::release() hands off ownership without deleting anything;
  // `storage` itself is emptied of ownership responsibility as it goes out
  // of scope right after (its destructor runs on now-null unique_ptr's).
  std::vector<SpatialDiscretization<T, DIM>*> result;
  result.reserve(N);
  for (auto& ptr : storage) {
    result.push_back(ptr.release());
  }

  return result;
}

/**
 * @brief Release a `SPAS`-shaped vector produced by `setSpatialDiscretization`.
 *
 * @details Deletes elements in reverse order (`spatials.back()` down to
 *          `spatials.front()`), so the mesh-owning element (`spatials[0]`,
 *          built from the mesh-creating constructor) is always deleted
 *          last, after every element referencing its mesh has already been
 *          destroyed. Clears `spatials` afterwards so it cannot be used
 *          again by mistake.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param spatials Vector produced by `setSpatialDiscretization`. Left empty
 *                after the call.
 */
template <class T, int DIM>
void deleteSpatialDiscretization(std::vector<SpatialDiscretization<T, DIM>*>& spatials) {
  for (auto it = spatials.rbegin(); it != spatials.rend(); ++it) {
    delete *it;
  }
  spatials.clear();
}