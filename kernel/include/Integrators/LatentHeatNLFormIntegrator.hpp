/**
 * @file LatentHeatNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Latent heat contribution for heat transfer equation
 * @version 0.1
 * @date 2025-09-05
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

#include <algorithm>
#include <list>
#include <memory>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Latent heat contribution for heat transfer equations
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class LatentHeatNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  std::vector<mfem::ParGridFunction> phi_gf_;
  std::vector<mfem::ParGridFunction> phi_old_gf_;
  double latent_time_step_;

  std::list<GlossaryType> expected_list_{GlossaryType::Mobility};
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

 protected:
  Coefficients mobility;

  double get_latent_heat_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir,
                               unsigned int blk);

  double compute_coefficient(Coefficient coef, const std::span<const double>& values);

  void get_coefficients() override;
  void check_variables_consistency();

 public:
  LatentHeatNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                             const std::vector<mfem::ParGridFunction>& aux_old,
                             const Parameters& params, std::vector<VARS*> auxvars,
                             const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
  void init() override;
};

#include "Integrators/LatentHeatNLFormIntegrator.tpp"
