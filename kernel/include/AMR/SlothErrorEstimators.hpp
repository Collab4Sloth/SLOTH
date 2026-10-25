/**
 * @file SlothErrorEstimators.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class used to build MFEM (or custom) error estimators at
 *        runtime, selected via an ErrorEstimatorType
 * @version 0.1
 * @date 2026-08-05
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
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Base class used to build MFEM (or custom) error estimators at
 *        runtime, selected via an ErrorEstimatorType
 */
class SlothErrorEstimators {
 private:
  /**
   * @brief Finite element space used by an estimator (flux space...), with its
   *        collection. It does not depend on the solution field: it is shared
   *        by all the fields of the same mesh and polynomial order, and reused
   *        along the whole computation (updated when the mesh has changed).
   *
   * Members are declared so that the default destruction order is valid:
   * space, then collection.
   */
  struct CachedSpace {
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    std::unique_ptr<mfem::ParFiniteElementSpace> fes;
  };
  /// Key of a cached space: mesh, polynomial order of the solution field, name of the space
  using SpaceKey = std::tuple<const mfem::ParMesh*, int, std::string>;

  ErrorEstimatorType value_;
  mfem::BilinearFormIntegrator* integ_{nullptr};
  std::map<SpaceKey, CachedSpace> spaces_;

  mfem::ParFiniteElementSpace& space(
      const std::string& name, mfem::ParMesh& mesh, const int order, const int vdim,
      const std::function<std::unique_ptr<mfem::FiniteElementCollection>()>& make_fec);

 public:
  SlothErrorEstimators(ErrorEstimatorType value, mfem::BilinearFormIntegrator* integ);
  explicit SlothErrorEstimators(ErrorEstimatorType value);

  SlothErrorEstimators(const SlothErrorEstimators&) = delete;
  SlothErrorEstimators& operator=(const SlothErrorEstimators&) = delete;

  std::shared_ptr<mfem::ErrorEstimator> get_value(mfem::ParGridFunction& x,
                                                  mfem::ParFiniteElementSpace& fespace,
                                                  mfem::ParMesh& mesh);

  /// Number of finite element spaces kept (all meshes, orders and names)
  std::size_t number_of_cached_spaces() const { return spaces_.size(); }

  ~SlothErrorEstimators() = default;
};
