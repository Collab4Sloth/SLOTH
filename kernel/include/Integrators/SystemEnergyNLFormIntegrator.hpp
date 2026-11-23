/**
 * @file SystemEnergyNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Integrator computing the energies of the system (no residual, no Jacobian)
 * @version 0.1
 * @date 2026-11-23
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
#include <vector>

#include "Coefficients/Coefficients.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Integrator of the energies of the system: SlothNLFormIntegrator::system_energies with
 * the coefficients of the operator (FreeEnergy and GradientEnergy of block 0). It is not added to
 * the nonlinear forms: it has no residual and no Jacobian.
 *
 * @tparam VARS Variables
 */
template <class VARS>
class SystemEnergyNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 public:
  /**
   * @param geometry Geometry of the problem
   * @param u Current values of the variables
   * @param params Parameters of the operator
   * @param auxvars Auxiliary variables (their current values are used)
   * @param coefficients Coefficients of the operator
   */
  SystemEnergyNLFormIntegrator(Geometry geometry, const std::vector<mfem::ParGridFunction>& u,
                               const Parameters& params, std::vector<VARS*> auxvars,
                               const std::vector<Coefficients>& coefficients)
      : SlothNLFormIntegrator<VARS>(geometry, 0., u, {}, params, auxvars, coefficients) {
    this->integrator_name_ = "SystemEnergyNLFormIntegrator";
  }

  void init() override {}

 protected:
  void get_coefficients() override {}

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& /*el*/,
                             mfem::ElementTransformation& /*Tr*/,
                             const mfem::Array<const mfem::Vector*>& /*elfun*/,
                             const mfem::Array<mfem::Vector*>& /*elvect*/) override {
    MFEM_ABORT("SystemEnergyNLFormIntegrator: no residual, energies of the system only.");
  }

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& /*el*/,
                           mfem::ElementTransformation& /*Tr*/,
                           const mfem::Array<const mfem::Vector*>& /*elfun*/,
                           const mfem::Array2D<mfem::DenseMatrix*>& /*elmat*/) override {
    MFEM_ABORT("SystemEnergyNLFormIntegrator: no Jacobian, energies of the system only.");
  }
};
