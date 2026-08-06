/**
 * @file SlothErrorEstimators.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class used to build MFEM (or custom) error estimators at
 *        runtime, selected via an ErrorEstimatorType
 * @version 0.1
 * @date 2026-08-05
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
 * @brief Build a fresh error estimator instance for the current state
 *        of `x`/`fespace`/`mesh`.
 *
 * @param x Grid function (solution field) the error is estimated on.
 * @param fespace Finite element space `x` lives on.
 * @param mesh Mesh `fespace` is built on.
 * @return std::shared_ptr<mfem::ErrorEstimator> Newly built estimator,
 *        owned by the caller.
 */
std::shared_ptr<mfem::ErrorEstimator> SlothErrorEstimators::get_value(
    mfem::ParGridFunction& x, mfem::ParFiniteElementSpace& fespace, mfem::ParMesh& mesh) {
  switch (this->value_) {
    case ErrorEstimatorType::L2_ZIENKIEWICZ_ZHU: {
      MFEM_VERIFY(this->integ_ != nullptr,
                  "SlothErrorEstimators: L2_ZIENKIEWICZ_ZHU requires a BilinearFormIntegrator, "
                  "pass one via the (type, integ) constructor.");

      const int order = fespace.GetOrder(0);
      this->flux_fes_ = new mfem::ParFiniteElementSpace(
          &mesh, new mfem::L2_FECollection(order, mesh.Dimension()), mesh.SpaceDimension());
      this->smooth_flux_fes_ = new mfem::ParFiniteElementSpace(
          &mesh, new mfem::RT_FECollection(order - 1, mesh.Dimension()));
      return std::make_shared<mfem::L2ZienkiewiczZhuEstimator>(*this->integ_, x, this->flux_fes_,
                                                               this->smooth_flux_fes_);
      break;
    }
    case ErrorEstimatorType::KELLY: {
      MFEM_VERIFY(this->integ_ != nullptr,
                  "SlothErrorEstimators: KELLY requires a BilinearFormIntegrator.");
      const int order = fespace.GetOrder(0);
      this->flux_fes_ = new mfem::ParFiniteElementSpace(
          &mesh, new mfem::L2_FECollection(order, mesh.Dimension()), mesh.SpaceDimension());
      return std::make_shared<mfem::KellyErrorEstimator>(*this->integ_, x, this->flux_fes_);
      break;
    }
    default:
      mfem::mfem_error("SlothEstimator::get_value: unknown estimator type.");
  }
  return nullptr;
}

/**
 * @brief Resynchronize the internal flux finite element spaces of the
 *        last estimator built by get_value(), after a mesh mutation
 *        that happened while that estimator is still alive.
 *
 * @details No-op if the current estimator type does not use
 *          separate flux spaces (e.g. not yet relevant for Kelly).
 */
void SlothErrorEstimators::UpdateFluxSpaces() {
  if (this->flux_fes_ != nullptr) {
    this->flux_fes_->Update();
  }
  if (this->smooth_flux_fes_ != nullptr) {
    this->smooth_flux_fes_->Update();
  }
}

/**
 * @brief Destroy the SlothErrorEstimators object
 *
 */
SlothErrorEstimators::~SlothErrorEstimators() {}
