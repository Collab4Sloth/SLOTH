/**
 * @file RobinNLFormIntegrator.hpp
 * @author Marine Harel (marine.harel@cea.fr)
 * @brief Robin integrator
 * @version 0.1
 * @date 2026-02-23
 *
 * @copyright  CEA (C) 2026
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
#include <string>
#include <utility>
#include <vector>

#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief  Class dedicated to the Robin boundary condition integrator
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class RobinNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi;

 protected:
  std::vector<mfem::ParGridFunction> vaux_gf_;
  const unsigned int bdr_id_ = 0;  // Boundary id
  const unsigned int blk_ = 0;     // Block to which apply Robin
  Coefficient robin_a;
  Coefficient robin_b;

  double compute_coefficient(Coefficient coef, const std::span<const double>& values,
                             const std::span<const double>& aux_values);
  double compute_gradient_coefficient(Coefficient coef, const int blk,
                                      const std::span<const double>& values,
                                      const std::span<const double>& aux_values);

  void get_coefficients() override;
  void init() override;

 public:
  RobinNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                        const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
                        std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients,
                        const unsigned int block, const unsigned int bdr_id);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
};

#include "Integrators/RobinNLFormIntegrator.tpp"
