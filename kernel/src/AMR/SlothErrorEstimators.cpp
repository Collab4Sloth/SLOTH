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
 * @brief Finite element space `name` for the given mesh and polynomial order,
 *        built at the first request (collection given by `make_fec`, `vdim`
 *        components) and updated (Update(false)) if the mesh has changed since
 *        its last use (refinement, derefinement or rebalancing).
 *
 * @param name Name of the space (spaces of the same name are shared).
 * @param mesh Mesh of the solution field.
 * @param order Polynomial order of the finite element space of the solution.
 * @param vdim Number of components of the space.
 * @param make_fec Builds the finite element collection (first request only).
 */
mfem::ParFiniteElementSpace& SlothErrorEstimators::space(
    const std::string& name, mfem::ParMesh& mesh, const int order, const int vdim,
    const std::function<std::unique_ptr<mfem::FiniteElementCollection>()>& make_fec) {
  auto& cached = this->spaces_[{&mesh, order, name}];
  if (cached.fes == nullptr) {
    cached.fec = make_fec();
    cached.fes = std::make_unique<mfem::ParFiniteElementSpace>(&mesh, cached.fec.get(), vdim);
  } else if (cached.fes->GetSequence() != mesh.GetSequence()) {
    cached.fes->Update(false);
  }
  return *cached.fes;
}

/**
 * @brief Return an error estimator of `x`, ready to compute the local errors
 *        of the current solution on the current mesh.
 *
 * @details The spaces used by the estimators (flux spaces) are built once per
 *          (mesh, polynomial order) and shared by all the solution fields (see
 *          space()); the mfem estimator itself, which only refers to `x` and to
 *          these spaces, is built at each call. Nothing is kept for `x`: a field
 *          may be destroyed at any time.
 *
 *          A new estimator is added as a new case below, requesting the spaces
 *          it needs with space().
 *
 *          The returned estimator refers to spaces owned by this object: it must
 *          not be used after the destruction of this object.
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
  const int dim = mesh.Dimension();
  const int sdim = mesh.SpaceDimension();
  // Discontinuous flux (gradient of the solution), shared by the estimators below
  auto l2_flux = [&]() -> mfem::ParFiniteElementSpace& {
    return this->space("L2 flux", mesh, order, sdim,
                       [&] { return std::make_unique<mfem::L2_FECollection>(order, dim); });
  };

  switch (this->value_) {
    case ErrorEstimatorType::L2_ZIENKIEWICZ_ZHU: {
      MFEM_VERIFY(this->integ_ != nullptr,
                  "SlothErrorEstimators: L2_ZIENKIEWICZ_ZHU requires an integrator.");
      auto& smooth_flux = this->space("RT smoothed flux", mesh, order, 1, [&] {
        return std::make_unique<mfem::RT_FECollection>(order - 1, dim);
      });
      return std::make_shared<mfem::L2ZienkiewiczZhuEstimator>(*this->integ_, x, l2_flux(),
                                                               smooth_flux);
    }
    case ErrorEstimatorType::KELLY: {
      MFEM_VERIFY(this->integ_ != nullptr, "SlothErrorEstimators: KELLY requires an integrator.");
      return std::make_shared<mfem::KellyErrorEstimator>(*this->integ_, x, l2_flux());
    }
    default:
      mfem::mfem_error("SlothErrorEstimators::get_value: unknown estimator type.");
  }
  return nullptr;
}
