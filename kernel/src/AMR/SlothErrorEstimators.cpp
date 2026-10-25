/**
 * @file SlothErrorEstimators.cpp
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
#include "AMR/SlothErrorEstimators.hpp"

#include <memory>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new SlothErrorEstimators object for an estimator
 *        type that requires a BilinearFormIntegrator (e.g.
 *        L2_ZIENKIEWICZ_ZHU).
 *
 * @param value Type of estimator to build (see ErrorEstimatorType).
 * @param integ Integrator used by the underlying estimator to compute
 *             per-element flux/error. The caller retains ownership and
 *             must keep it alive for as long as this object is used.
 */
SlothErrorEstimators::SlothErrorEstimators(ErrorEstimatorType value,
                                           mfem::BilinearFormIntegrator* integ)
    : value_(value), integ_(integ) {}

/**
 * @brief Construct a new SlothErrorEstimators object for an estimator
 *        type that requires no auxiliary integrator (e.g. a custom
 *        gradient-based estimator).
 *
 * @param value Type of estimator to build (see ErrorEstimatorType).
 */
SlothErrorEstimators::SlothErrorEstimators(ErrorEstimatorType value)
    : value_(value), integ_(nullptr) {}

/**
 * @brief Release the estimator, then the flux spaces, then their
 *        collections (the estimator references the spaces, which
 *        reference the collections).
 */
void SlothErrorEstimators::CachedEstimator::release() {
  this->estimator.reset();
  this->smooth_flux_fes.reset();
  this->flux_fes.reset();
  this->smooth_flux_fec.reset();
  this->flux_fec.reset();
  this->mesh = nullptr;
  this->order = -1;
}

/**
 * @brief Build the flux spaces and the mfem estimator of one solution
 *        field. The spaces are passed by reference to the mfem estimator,
 *        which therefore does not own them.
 *
 * @param entry Cache entry to fill (released first).
 * @param x Grid function (solution field) the error is estimated on.
 * @param order Polynomial order of the finite element space of `x`.
 * @param mesh Mesh the finite element space of `x` is built on.
 */
void SlothErrorEstimators::build(CachedEstimator& entry, mfem::ParGridFunction& x, const int order,
                                 mfem::ParMesh& mesh) {
  entry.release();
  entry.mesh = &mesh;
  entry.order = order;
  const int dim = mesh.Dimension();
  const int sdim = mesh.SpaceDimension();

  switch (this->value_) {
    case ErrorEstimatorType::L2_ZIENKIEWICZ_ZHU: {
      MFEM_VERIFY(this->integ_ != nullptr,
                  "SlothErrorEstimators: L2_ZIENKIEWICZ_ZHU requires an integrator.");
      entry.flux_fec = std::make_unique<mfem::L2_FECollection>(order, dim);
      entry.smooth_flux_fec = std::make_unique<mfem::RT_FECollection>(order - 1, dim);
      entry.flux_fes =
          std::make_unique<mfem::ParFiniteElementSpace>(&mesh, entry.flux_fec.get(), sdim);
      entry.smooth_flux_fes =
          std::make_unique<mfem::ParFiniteElementSpace>(&mesh, entry.smooth_flux_fec.get());
      entry.estimator = std::make_shared<mfem::L2ZienkiewiczZhuEstimator>(
          *this->integ_, x, *entry.flux_fes, *entry.smooth_flux_fes);
      break;
    }
    case ErrorEstimatorType::KELLY: {
      MFEM_VERIFY(this->integ_ != nullptr, "SlothErrorEstimators: KELLY requires an integrator.");
      entry.flux_fec = std::make_unique<mfem::L2_FECollection>(order, dim);
      entry.flux_fes =
          std::make_unique<mfem::ParFiniteElementSpace>(&mesh, entry.flux_fec.get(), sdim);
      entry.estimator =
          std::make_shared<mfem::KellyErrorEstimator>(*this->integ_, x, *entry.flux_fes);
      break;
    }
    default:
      mfem::mfem_error("SlothErrorEstimators::get_value: unknown estimator type.");
  }
}

/**
 * @brief Return the error estimator of `x`, ready to compute the local
 *        errors of the current solution on the current mesh.
 *
 * @details The estimator and its flux spaces are built once per solution
 *          field, and reused at each AMR pass:
 *          - no allocation of collections/spaces at each pass;
 *          - after a mesh modification, the flux spaces are updated by the
 *            mfem estimator itself (`ComputeEstimates()` calls
 *            `Update(false)`, which is valid whatever the number of mesh
 *            modifications since the previous estimate);
 *          - `Reset()` is called so that the estimates are recomputed even
 *            when the mesh is unchanged (the solution has changed).
 *          The estimator is rebuilt only if the mesh or the polynomial
 *          order associated with `x` has changed.
 *
 *          The returned estimator is owned by this object: it must not be
 *          used after the next call to get_value() for the same `x`, nor
 *          after the destruction of this object.
 *
 * @param x Grid function (solution field) the error is estimated on.
 * @param fespace Finite element space `x` lives on.
 * @param mesh Mesh `fespace` is built on.
 * @return std::shared_ptr<mfem::ErrorEstimator> Estimator of `x`.
 */
std::shared_ptr<mfem::ErrorEstimator> SlothErrorEstimators::get_value(
    mfem::ParGridFunction& x, mfem::ParFiniteElementSpace& fespace, mfem::ParMesh& mesh) {
  // Order from the collection: valid even on a rank without any element
  const int order = fespace.FEColl()->GetOrder();

  auto& entry = this->cache_[&x];
  if (entry.estimator == nullptr || entry.mesh != &mesh || entry.order != order) {
    this->build(entry, x, order, mesh);
  } else {
    entry.estimator->Reset();
  }
  return entry.estimator;
}

/**
 * @brief Destroy the SlothErrorEstimators object
 *
 */
SlothErrorEstimators::~SlothErrorEstimators() {
  for (auto& [gf, entry] : this->cache_) {
    entry.release();
  }
}
