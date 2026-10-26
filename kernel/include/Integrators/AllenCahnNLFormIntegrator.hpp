/**
 * @file AllenCahnNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the Allen Cahn equation
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
 * @brief Class dedicated to the VF of the RHS of the AllenCahn equation (see TimeDerivative
 * integrator for the LHS)
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class AllenCahnNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  std::list<GlossaryType> expected_list_{GlossaryType::Mobility, GlossaryType::Capillary,
                                         GlossaryType::FreeEnergy};
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  std::vector<double> mob_at_ip_;   ///< mobility of each block at the current point
  std::vector<double> dmob_at_ip_;  ///< dM/du_j at the current point

  // General case: blocks with different finite elements (e.g. different orders). Same
  // formulas as AssembleElementVector/AssembleElementGrad, with the shape functions of each
  // block and rectangular off-diagonal element matrices.
  std::vector<mfem::Vector> shapes_;  ///< shape functions of each block at the current point
  std::vector<double> u_generic_;     ///< {u, u_old} at the current point (general case)
  std::vector<double> aux_generic_;   ///< auxiliary values at the current point (general case)

  const mfem::IntegrationRule& generic_integration_rule(
      const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr) const;
  void generic_values_at_ip(const mfem::Array<const mfem::FiniteElement*>& el,
                            mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                            const mfem::Array<const mfem::Vector*>& elfun);

  void AssembleElementVectorGeneric(const mfem::Array<const mfem::FiniteElement*>& el,
                                    mfem::ElementTransformation& Tr,
                                    const mfem::Array<const mfem::Vector*>& elfun,
                                    const mfem::Array<mfem::Vector*>& elvec);
  void AssembleElementGradGeneric(const mfem::Array<const mfem::FiniteElement*>& el,
                                  mfem::ElementTransformation& Tr,
                                  const mfem::Array<const mfem::Vector*>& elfun,
                                  const mfem::Array2D<mfem::DenseMatrix*>& elmats);

 protected:
  Coefficients lambda;
  Coefficients mobility;
  Coefficients double_well_energy;

  std::vector<mfem::ParGridFunction> vaux_gf_;

  void get_coefficients() override;

  virtual void check_variables_consistency();

  void add_mobility_derivative(const int blk, const double lamb, const mfem::DenseMatrix& gradPsi,
                               const mfem::Vector& Psi, const mfem::Vector& elfun_blk,
                               const std::span<const double>& u, const std::span<const double>& aux,
                               const double weight_coef,
                               const mfem::Array2D<mfem::DenseMatrix*>& elmats,
                               const std::vector<mfem::Vector>* col_shapes = nullptr);

  virtual double get_mob_at_ip(unsigned int blk, const std::span<const double>& values,
                               const std::span<const double>& aux_values);

 public:
  void init() override;
  AllenCahnNLFormIntegrator(Geometry geometry, const double time_step,
                            const std::vector<mfem::ParGridFunction>& u_old,
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
};

#include "Integrators/AllenCahnNLFormIntegrator.tpp"
