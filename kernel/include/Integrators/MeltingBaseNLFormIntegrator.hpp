/**
 * @file MeltingBaseNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for melting integrators
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
 * @brief Base class for melting integrators
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class MeltingBaseNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  std::list<GlossaryType> expected_list_{GlossaryType::Mobility, GlossaryType::PhaseFieldPotential};
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

 protected:
  Coefficients mobility;
  Coefficients interpolation_potential;
  std::vector<mfem::ParGridFunction> vaux_gf_;

  virtual double get_phase_change_at_ip(unsigned int blk, const std::span<const double>& values,
                                        const std::span<const double>& aux_values) = 0;

  virtual double get_seed_at_ip(unsigned int blk, const std::span<const double>& values,
                                const std::span<const double>& aux_values);

  void get_coefficients() override;

 public:
  MeltingBaseNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
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

#include "Integrators/MeltingBaseNLFormIntegrator.tpp"
