/**
 * @file AdvectionBoundaryNLFormIntegrator.tpp
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

#include <algorithm>
#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/AdvectionBoundaryNLFormIntegrator.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AdvectionBoundaryNLFormIntegrator object (all roles = block).
 *
 * @param geometry     Geometry of the problem (Cartesian, axisymmetric, ...).
 * @param time_step    Current time step.
 * @param u_old        Vector of previous-time-step solution fields.
 * @param aux_old      Vector of previous-time-step auxiliary fields.
 * @param params       Parameters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the integrator.
 * @param coefficients List of coefficients (Velocity components, Inflow).
 * @param block        Block to which the boundary term is applied.
 * @param bdr_id       Boundary id.
 */
template <class VARS>
AdvectionBoundaryNLFormIntegrator<VARS>::AdvectionBoundaryNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients,
    const unsigned int block, const unsigned int bdr_id)
    : AdvectionBoundaryNLFormIntegrator(geometry, time_step, u_old, aux_old, params,
                                        std::move(auxvars), coefficients, block, block, block,
                                        block, bdr_id) {}

/**
 * @brief Protected constructor with explicit velocity block, residual row, advected variable
 * and inflow block.
 */
template <class VARS>
AdvectionBoundaryNLFormIntegrator<VARS>::AdvectionBoundaryNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients,
    const unsigned int vel_blk, const unsigned int row_blk, const unsigned int var_blk,
    const unsigned int inflow_blk, const unsigned int bdr_id)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, std::move(auxvars),
                                  coefficients),
      vel_blk_(vel_blk),
      row_blk_(row_blk),
      var_blk_(var_blk),
      inflow_blk_(inflow_blk),
      bdr_id_(bdr_id),
      inflow_(Coefficient(Glossary::Default, 0.0)) {
  this->integrator_name_ = "AdvectionBoundary";
}

/**
 * @brief Initialize the integrator: coefficients and auxiliary fields.
 */
template <class VARS>
void AdvectionBoundaryNLFormIntegrator<VARS>::init() {
  const unsigned int max_blk =
      std::max({this->vel_blk_, this->row_blk_, this->var_blk_, this->inflow_blk_});
  MFEM_VERIFY(max_blk < this->nb_blk_,
              this->integrator_name_ + "NLFormIntegrator: block index out of range.");

  this->get_coefficients();

  this->vaux_gf_.reserve(this->vaux_gf_.size() + this->aux_infos_.size());
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief Configure the velocity of the advected block (see VelocityField) and retrieve the
 * Inflow coefficient of the boundary, or raise an error.
 */
template <class VARS>
void AdvectionBoundaryNLFormIntegrator<VARS>::get_coefficients() {
  this->sdim_ = this->u_old_[0].ParFESpace()->GetParMesh()->SpaceDimension();

  std::string grad_symbol;
  if (this->params_.has_parameter("velocity_gradient")) {
    grad_symbol = this->params_.template get_param_value<std::string>("velocity_gradient");
  }
  const auto get_coef = [this](GlossaryType type, int id) {
    return this->get_coefficient(this->vel_blk_, type, id);
  };
  const bool has_velocity = this->velocity_.configure(
      this->integrator_name_, this->vel_blk_, this->sdim_, get_coef, grad_symbol, this->aux_infos_);
  MFEM_VERIFY(has_velocity, this->integrator_name_ +
                                "NLFormIntegrator: no velocity (Velocity or VelocityPrefactor "
                                "coefficients) for block " +
                                std::to_string(this->vel_blk_) + ".");
  this->needs_aux_ = this->velocity_.needs_aux_values();

  auto inflow = this->get_coefficient(this->inflow_blk_, GlossaryType::Inflow, 0, this->bdr_id_);
  MFEM_VERIFY(inflow.has_value(), this->integrator_name_ +
                                      "NLFormIntegrator: Inflow coefficient not found for "
                                      "boundary " +
                                      std::to_string(this->bdr_id_) + ".");
  this->inflow_ = *inflow;
  this->needs_aux_ = this->needs_aux_ || (!inflow->is_scalar() && !inflow->is_file());
}

/**
 * @brief Integration rule shared by all blocks (exact for the highest-order block).
 */
template <class VARS>
const mfem::IntegrationRule& AdvectionBoundaryNLFormIntegrator<VARS>::get_integration_rule(
    const mfem::Array<const mfem::FiniteElement*>& el, const mfem::ElementTransformation& Tr) {
  int order = 0;
  for (int b = 0; b < el.Size(); ++b) {
    order = std::max(order, el[b]->GetOrder());
  }
  return mfem::IntRules.Get(el[0]->GetGeomType(), 2 * order + Tr.OrderW());
}

/**
 * @brief Element-level setup: size the work arrays.
 */
template <class VARS>
void AdvectionBoundaryNLFormIntegrator<VARS>::prepare_element(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr) {
  const int nb = el.Size();
  this->shape_.resize(nb);
  this->u_values_.assign(2 * nb, 0.0);
  this->aux_values_.assign(this->vaux_gf_.size(), 0.0);
  for (int b = 0; b < nb; ++b) {
    this->shape_[b].SetSize(el[b]->GetDof());
  }
  this->nor_.SetSize(Tr.GetSpaceDim());
  this->v_.SetSize(this->sdim_);
}

/**
 * @brief Evaluate shapes, unknowns, auxiliary fields (if needed), normal and velocity at ip.
 *
 * @return weight * (v.n). The normal from CalcOrtho is not normalized and carries the surface
 *         measure, hence no Tr.Weight() factor. Axisymmetric factor included.
 */
template <class VARS>
double AdvectionBoundaryNLFormIntegrator<VARS>::evaluate_at_ip(
    const mfem::Array<const mfem::FiniteElement*>& el,
    const mfem::Array<const mfem::Vector*>& elfun, mfem::ElementTransformation& Tr,
    const mfem::IntegrationPoint& ip) {
  this->set_int_point(Tr, ip);

  // Normal first, while Tr is set at ip
  if (Tr.GetSpaceDim() > 1) {
    mfem::CalcOrtho(Tr.Jacobian(), this->nor_);
  } else {
    // 1D: boundary points, assumption: attribute 1 = left, otherwise right
    this->nor_[0] = (Tr.Attribute == 1) ? -1.0 : 1.0;
  }

  const int nb = el.Size();
  for (int b = 0; b < nb; ++b) {
    el[b]->CalcShape(ip, this->shape_[b]);
    this->u_values_[b] = (*elfun[b]) * this->shape_[b];
    this->u_values_[b + nb] = this->u_old_[b].GetValue(Tr, ip);
  }
  if (this->needs_aux_) {
    for (std::size_t k = 0; k < this->vaux_gf_.size(); ++k) {
      this->aux_values_[k] = this->vaux_gf_[k].GetValue(Tr, ip);
    }
  }

  const std::span<const double> u(this->u_values_);
  const std::span<const double> aux(this->aux_values_);
  const auto coef = [this, &u, &aux](const Coefficient& c) {
    Coefficient copy = c;  // eval_coefficient may set the coordinates of the coefficient
    return this->eval_coefficient(copy, u, aux);
  };
  this->velocity_.eval(Tr, this->vaux_gf_, coef, this->v_);

  double weight = ip.weight;
  if (this->isAxisymmetric()) {
    weight *= AxiCylindricalCoefficient().Eval(Tr, ip);
  }
  return weight * (this->v_ * this->nor_);
}

/**
 * @brief Assemble the boundary residual: R_r = ( (v.n)^+ u_c + (v.n)^- u_in, psi_r ).
 */
template <class VARS>
void AdvectionBoundaryNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int nb = el.Size();
  this->prepare_element(el, Tr);

  for (int b = 0; b < nb; ++b) {
    elvect[b]->SetSize(el[b]->GetDof());
    *elvect[b] = 0.0;
  }

  const std::span<const double> u(this->u_values_);
  const std::span<const double> aux(this->aux_values_);
  const mfem::IntegrationRule& ir = get_integration_rule(el, Tr);

  for (int q = 0; q < ir.GetNPoints(); ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    const double w_vn = this->evaluate_at_ip(el, elfun, Tr, ip);  // weight * (v.n)
    if (w_vn == 0.0) {
      continue;  // wall
    }
    // Outflow: interior value; inflow: imposed value
    const double advected = (w_vn > 0.0) ? this->u_values_[this->var_blk_]
                                         : this->eval_coefficient(this->inflow_, u, aux);
    elvect[this->row_blk_]->Add(w_vn * advected, this->shape_[this->row_blk_]);
  }
}

/**
 * @brief Assemble the boundary Jacobian: J_rc = ( (v.n)^+ phi_c, psi_r ).
 *
 * Only block (r, c) is non-zero, but all blocks are sized and reset (no stale contribution
 * from other integrators).
 */
template <class VARS>
void AdvectionBoundaryNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int nb = el.Size();
  this->prepare_element(el, Tr);

  for (int i = 0; i < nb; ++i) {
    for (int j = 0; j < nb; ++j) {
      elmats(i, j)->SetSize(el[i]->GetDof(), el[j]->GetDof());
      *elmats(i, j) = 0.0;
    }
  }

  const mfem::IntegrationRule& ir = get_integration_rule(el, Tr);

  for (int q = 0; q < ir.GetNPoints(); ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    const double w_vn = this->evaluate_at_ip(el, elfun, Tr, ip);  // weight * (v.n)
    if (w_vn <= 0.0) {
      continue;  // wall or inflow: no dependence on the unknowns
    }
    mfem::AddMult_a_VWt(w_vn, this->shape_[this->row_blk_], this->shape_[this->var_blk_],
                        *elmats(this->row_blk_, this->var_blk_));
  }
}
