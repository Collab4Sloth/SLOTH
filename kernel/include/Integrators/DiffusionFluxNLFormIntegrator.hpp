/**
 * @file DiffusionFluxNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief inter-diffusion integrator
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

#pragma once

#include <list>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief  Class dedicated to the VF of an inter-diffusion equation
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class DiffusionFluxNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  void add_diffusion_flux(mfem::ElementTransformation& Tr, const int nElement,
                          const mfem::IntegrationPoint& ip, const int dim);

 protected:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, Flux_;

  std::list<GlossaryType> expected_list_;
  Coefficients stab_diffusion;

  double compute_coefficient(Coefficient coef, const std::span<const double>& values);

  std::vector<SlothGridFunction> sloth_u_old_;
  virtual void get_parameters() {}

  virtual std::vector<mfem::Vector> get_flux_gradient(mfem::ElementTransformation& Tr,
                                                      const int nElement,
                                                      const mfem::IntegrationPoint& ip,
                                                      const int dim) = 0;
  virtual std::vector<double> get_flux_coefficient(const int nElement,
                                                   const mfem::IntegrationPoint& ip) = 0;
  void get_coefficients() override;
  void init() override;

 public:
  DiffusionFluxNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                const std::vector<mfem::ParGridFunction>& aux_old,
                                const Parameters& params, std::vector<VARS*> auxvars,
                                const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvect) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmat) override;
};

#include "Integrators/DiffusionFluxNLFormIntegrator.tpp"
