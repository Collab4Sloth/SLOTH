/**
 * @file SlothErrorEstimators.hpp
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
#pragma once
#include <memory>

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
  ErrorEstimatorType value_;
  mfem::BilinearFormIntegrator* integ_{nullptr};
  mfem::ParFiniteElementSpace* flux_fes_{nullptr};
  mfem::ParFiniteElementSpace* smooth_flux_fes_{nullptr};

 public:
  SlothErrorEstimators(ErrorEstimatorType value, mfem::BilinearFormIntegrator* integ);
  explicit SlothErrorEstimators(ErrorEstimatorType value);

  std::shared_ptr<mfem::ErrorEstimator> get_value(mfem::ParGridFunction& x,
                                                  mfem::ParFiniteElementSpace& fespace,
                                                  mfem::ParMesh& mesh);
  void UpdateFluxSpaces();
  ~SlothErrorEstimators();
};
