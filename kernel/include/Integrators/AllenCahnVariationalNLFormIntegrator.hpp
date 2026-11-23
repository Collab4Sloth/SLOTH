/**
 * @file AllenCahnVariationalNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Allen-Cahn integrator with a variational free energy f(u, grad u)
 * @version 0.1
 * @date 2026-11-22
 *
 * Copyright CEA (C) 2026
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
#include <optional>
#include <vector>

#include "Coefficients/FieldsAtPoints.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Coefficients/VariationalEval.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Right-hand side of the Allen-Cahn equations with a variational free energy (see the
 * TimeDerivative integrator for the left-hand side).
 *
 * The free energy F = int f(u, grad u) is a variational coefficient ("variational": true in the
 * JSON file) of the n blocks u_0, ..., u_{n-1}, gradient terms included. Equation of the block b:
 *   du_b/dt = - M_b dF/du_b,   dF/du_b = df/du_b - div(df/d(grad u_b)),
 * weak form of the right-hand side, for a test function psi:
 *   R_b(psi) = int M_b (df/du_b psi + df/d(grad u_b) . grad psi).
 * Its Jacobian involves the second derivatives of f with respect to the variables and to the
 * components of their gradients, and the derivatives of the mobility.
 *
 * Coefficients (each block):
 *   - Mobility: mobility M_b of the block,
 *   - FreeEnergy: variational coefficient f. The energy of the first block is used for all the
 *     blocks: it is a function of the n blocks (in their order) and of their gradients.
 * The capillary coefficient is part of f: no Capillary coefficient.
 *
 * Scheme of the free energy: implicit (values and gradients of u), explicit (values and gradients
 * of u_old, no contribution to the Jacobian) or semi-implicit (variables u, explicit variables
 * u_old).
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class AllenCahnVariationalNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  /// Base class of the nonlinear integrators of SLOTH
  using SlothNLBase = SlothNLFormIntegrator<VARS>;
  using CoefficientAtPoints = typename SlothNLBase::CoefficientAtPoints;

  std::list<GlossaryType> expected_list_{GlossaryType::Mobility, GlossaryType::FreeEnergy};

  // Shape functions and their physical gradients at one point
  mfem::Vector Psi_;
  mfem::DenseMatrix gradPsi_;

  // Mobility of each block at the integration points of the current element
  std::vector<CoefficientAtPoints> mobility_at_points_;

  // Inputs of the evaluation of the free energy at the integration points of the current element
  std::vector<mfem::Vector> dofs_u_, dofs_u_old_;
  mfem::Vector values_, explicit_values_, auxiliary_values_;
  mfem::Vector gradients_, explicit_gradients_, auxiliary_gradients_, coordinates_;
  mfem::DenseMatrix points_;

  // Derivatives of the free energy at the integration points of the current element (layout of
  // VariationalBatch: values of one point contiguous)
  mfem::Vector df_dvar_, df_dgrad_, d2f_dvar2_, d2f_dvar_dgrad_, d2f_dgrad2_;
  bool with_second_derivatives_ = false;  ///< second derivatives computed (implicit energy)

  // Work data of the assembly
  mfem::Vector flux_, row_, column_;
  mfem::DenseMatrix block_hessian_, gradPsi_hessian_;
  std::vector<mfem::Vector> residual_integrand_;  ///< integrand of R_b at one point (each block)

  void energy_at_points(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                        const mfem::IntegrationRule& ir,
                        const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian);
  void mobility_at_points(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                          const mfem::IntegrationRule& ir,
                          const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian);
  double integration_weight(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip);

 protected:
  Coefficients mobility;
  std::optional<Coefficient> free_energy;

  std::vector<mfem::ParGridFunction> vaux_gf_;

  void get_coefficients() override;

 public:
  void init() override;
  AllenCahnVariationalNLFormIntegrator(Geometry geometry, const double time_step,
                                       const std::vector<mfem::ParGridFunction>& u_old,
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
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
};

#include "Integrators/AllenCahnVariationalNLFormIntegrator.tpp"
