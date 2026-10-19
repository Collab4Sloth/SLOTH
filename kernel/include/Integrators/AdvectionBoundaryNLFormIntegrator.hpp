/**
 * @file AdvectionBoundaryNLFormIntegrator.hpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Inflow/outflow boundary term of the conservative advection operator div(v u)
 * @version 0.5
 * @date 2026-10-17
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

#include <string>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Integrators/VelocityField.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Upwinded boundary term of the advection operator on boundary bdr_id:
 *
 *   R_r  = ( (v.n)^+ u_c + (v.n)^- u_in , psi_r )_Gamma
 *   J_rc = ( (v.n)^+ phi_c, psi_r )_Gamma
 *
 * with (a)^+ = max(a, 0), (a)^- = min(a, 0):
 *  - outflow (v.n > 0): the interior value u_c leaves the domain,
 *  - inflow  (v.n < 0): the incoming flux (v.n) u_in is imposed.
 *
 * v is given by the (domain) velocity coefficients of the advected block (see VelocityField:
 * Velocity_X/Y/Z or v = alpha grad(a)), u_in by the Inflow coefficient attached to the
 * boundary. Both can depend on the coordinates, the time and the auxiliary variables, not on the
 * unknowns.
 *
 * The integrator is added automatically (see OperatorBase::build_rhs_nonlinear_form) on each
 * boundary holding an Inflow coefficient. Default: r = c = block.
 * See AdvectionBoundaryCH for the split Cahn-Hilliard case.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class AdvectionBoundaryNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  // Work arrays reused across elements and integration points
  std::vector<mfem::Vector> shape_;
  std::vector<double> u_values_;
  std::vector<double> aux_values_;
  mfem::Vector nor_;
  mfem::Vector v_;

  void prepare_element(const mfem::Array<const mfem::FiniteElement*>& el,
                       mfem::ElementTransformation& Tr);
  double evaluate_at_ip(const mfem::Array<const mfem::FiniteElement*>& el,
                        const mfem::Array<const mfem::Vector*>& elfun,
                        mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip);
  static const mfem::IntegrationRule& get_integration_rule(
      const mfem::Array<const mfem::FiniteElement*>& el, const mfem::ElementTransformation& Tr);

 protected:
  std::vector<mfem::ParGridFunction> vaux_gf_;
  const unsigned int vel_blk_;     ///< block holding the Velocity coefficients
  const unsigned int row_blk_;     ///< residual row receiving the boundary term
  const unsigned int var_blk_;     ///< advected variable
  const unsigned int inflow_blk_;  ///< block holding the Inflow coefficient
  const unsigned int bdr_id_;      ///< boundary id

  int sdim_{0};
  VelocityField velocity_;
  Coefficient inflow_;
  bool needs_aux_{false};

  void get_coefficients() override;

  AdvectionBoundaryNLFormIntegrator(Geometry geometry, const double time_step,
                                    const std::vector<mfem::ParGridFunction>& u_old,
                                    const std::vector<mfem::ParGridFunction>& aux_old,
                                    const Parameters& params, std::vector<VARS*> auxvars,
                                    const std::vector<Coefficients>& coefficients,
                                    const unsigned int vel_blk, const unsigned int row_blk,
                                    const unsigned int var_blk, const unsigned int inflow_blk,
                                    const unsigned int bdr_id);

 public:
  AdvectionBoundaryNLFormIntegrator(Geometry geometry, const double time_step,
                                    const std::vector<mfem::ParGridFunction>& u_old,
                                    const std::vector<mfem::ParGridFunction>& aux_old,
                                    const Parameters& params, std::vector<VARS*> auxvars,
                                    const std::vector<Coefficients>& coefficients,
                                    const unsigned int block, const unsigned int bdr_id);

  void init() override;

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
};

#include "Integrators/AdvectionBoundaryNLFormIntegrator.tpp"
