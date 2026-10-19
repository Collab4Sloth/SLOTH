/**
 * @file AdvectionNLFormIntegrator.tpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the conservative advection term div(v u)
 * @version 0.5
 * @date 2026-10-17
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
#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AdvectionNLFormIntegrator object.
 *
 * @param geometry     Geometry of the problem (Cartesian, axisymmetric, ...).
 * @param time_step    Current time step.
 * @param u_old        Vector of previous-time-step solution fields.
 * @param aux_old      Vector of previous-time-step auxiliary fields.
 * @param params       Parameters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the integrator.
 * @param coefficients List of coefficients (Velocity components).
 */
template <class VARS>
AdvectionNLFormIntegrator<VARS>::AdvectionNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, std::move(auxvars),
                                  coefficients) {
  this->integrator_name_ = "Advection";
}

/**
 * @brief Initialize the integrator: velocity coefficients, auxiliary fields, terms.
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::init() {
  this->get_coefficients();

  this->vaux_gf_.reserve(this->vaux_gf_.size() + this->aux_infos_.size());
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }

  this->build_terms();
  MFEM_VERIFY(!this->terms_.empty(), this->integrator_name_ +
                                         "NLFormIntegrator: no advected variable (Velocity or "
                                         "VelocityPrefactor coefficients are expected). Please "
                                         "check your data.");
}

/**
 * @brief Configure the velocity of each block (see VelocityField).
 *
 * A block is advected if it holds all the components (Velocity_X, Velocity_Y[, Velocity_Z]),
 * or a VelocityPrefactor coefficient (v = alpha grad(a), parameter "velocity_gradient").
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::get_coefficients() {
  this->sdim_ = this->u_old_[0].ParFESpace()->GetParMesh()->SpaceDimension();
  this->grad_symbol_.clear();
  if (this->params_.has_parameter("velocity_gradient")) {
    this->grad_symbol_ = this->params_.template get_param_value<std::string>("velocity_gradient");
  }

  this->velocity_.assign(this->nb_blk_, VelocityField());
  this->needs_aux_ = false;
  for (unsigned int i = 0; i < this->nb_blk_; ++i) {
    const auto get_coef = [this, i](GlossaryType type, int id) {
      return this->get_coefficient(i, type, id);
    };
    this->velocity_[i].configure(this->integrator_name_, i, this->sdim_, get_coef,
                                 this->grad_symbol_, this->aux_infos_);
    this->needs_aux_ = this->needs_aux_ || this->velocity_[i].needs_aux_values();
  }
}

/**
 * @brief Default terms: (b, b) for each block b having Velocity coefficients.
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::build_terms() {
  this->terms_.clear();
  for (unsigned int b = 0; b < this->nb_blk_; ++b) {
    if (this->velocity_[b].active()) {
      const int ib = static_cast<int>(b);
      this->terms_.push_back({ib, ib, ib});
    }
  }
}

/**
 * @brief Velocity of block 'blk' at the current integration point (Tr set at ip).
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::compute_velocity(int blk, mfem::ElementTransformation& Tr,
                                                       const std::span<const double>& u,
                                                       const std::span<const double>& aux,
                                                       mfem::Vector& v) {
  const auto coef = [this, &u, &aux](const Coefficient& c) {
    return this->compute_coefficient(c, u, aux);
  };
  this->velocity_[blk].eval(Tr, this->vaux_gf_, coef, v);
}

/**
 * @brief Integration rule shared by all blocks (exact for the highest-order block).
 */
template <class VARS>
const mfem::IntegrationRule& AdvectionNLFormIntegrator<VARS>::get_integration_rule(
    const mfem::Array<const mfem::FiniteElement*>& el, const mfem::ElementTransformation& Tr) {
  int order = 0;
  for (int b = 0; b < el.Size(); ++b) {
    order = std::max(order, el[b]->GetOrder());
  }
  return mfem::IntRules.Get(el[0]->GetGeomType(), 2 * order + Tr.OrderW());
}

/**
 * @brief Element-level setup: size the work arrays and extract the local dofs of u_old once.
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::prepare_element(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr) {
  const int nb = el.Size();
  this->shape_.resize(nb);
  this->u_old_dofs_.resize(nb);
  this->u_values_.assign(2 * nb, 0.0);
  this->aux_values_.assign(this->vaux_gf_.size(), 0.0);
  this->v_.SetSize(this->sdim_);

  for (int b = 0; b < nb; ++b) {
    this->shape_[b].SetSize(el[b]->GetDof());
    this->u_old_[b].GetElementDofValues(Tr.ElementNo, this->u_old_dofs_[b]);
  }
}

/**
 * @brief Evaluate the shape functions, the unknowns and (if needed) the auxiliary fields at ip.
 *
 * @return The quadrature weight (including the axisymmetric factor if needed).
 */
template <class VARS>
double AdvectionNLFormIntegrator<VARS>::evaluate_at_ip(
    const mfem::Array<const mfem::FiniteElement*>& el,
    const mfem::Array<const mfem::Vector*>& elfun, mfem::ElementTransformation& Tr,
    const mfem::IntegrationPoint& ip) {
  this->set_int_point(Tr, ip);

  const int nb = el.Size();
  for (int b = 0; b < nb; ++b) {
    el[b]->CalcShape(ip, this->shape_[b]);
    this->u_values_[b] = (*elfun[b]) * this->shape_[b];
    this->u_values_[b + nb] = this->u_old_dofs_[b] * this->shape_[b];
  }

  if (this->needs_aux_) {
    for (std::size_t k = 0; k < this->vaux_gf_.size(); ++k) {
      this->aux_values_[k] = this->vaux_gf_[k].GetValue(Tr, ip);
    }
  }

  double weight = ip.weight * Tr.Weight();
  if (this->isAxisymmetric()) {
    weight *= AxiCylindricalCoefficient().Eval(Tr, ip);
  }
  return weight;
}

/**
 * @brief adv_ = v . grad(psi_r) for the row-block finite element.
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::compute_advection_shape(const mfem::FiniteElement& fe,
                                                              mfem::ElementTransformation& Tr) {
  this->dshape_.SetSize(fe.GetDof(), this->sdim_);
  this->adv_.SetSize(fe.GetDof());
  fe.CalcPhysDShape(Tr, this->dshape_);
  this->dshape_.Mult(this->v_, this->adv_);
}

/**
 * @brief Assemble the element-level residual: R_r = - ( u_c v, grad(psi_r) ).
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  // Catch_Time_Section("AdvectionNLFormIntegrator::AssembleElementVector");
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
    const double weight = this->evaluate_at_ip(el, elfun, Tr, ip);

    for (const auto& term : this->terms_) {
      this->compute_velocity(term.vel, Tr, u, aux, this->v_);
      this->compute_advection_shape(*el[term.row], Tr);
      elvect[term.row]->Add(-weight * this->u_values_[term.var], this->adv_);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian: J_rc = - ( v.grad(psi_r), phi_c ).
 *
 * Only blocks (r, c) are non-zero, but all blocks are sized and reset so that no stale
 * contribution left by another integrator in the shared elmats can be re-assembled.
 */
template <class VARS>
void AdvectionNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("AdvectionNLFormIntegrator::AssembleElementGrad");
  const int nb = el.Size();
  this->prepare_element(el, Tr);

  for (int i = 0; i < nb; ++i) {
    for (int j = 0; j < nb; ++j) {
      elmats(i, j)->SetSize(el[i]->GetDof(), el[j]->GetDof());
      *elmats(i, j) = 0.0;
    }
  }

  const std::span<const double> u(this->u_values_);
  const std::span<const double> aux(this->aux_values_);
  const mfem::IntegrationRule& ir = get_integration_rule(el, Tr);

  for (int q = 0; q < ir.GetNPoints(); ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    const double weight = this->evaluate_at_ip(el, elfun, Tr, ip);

    for (const auto& term : this->terms_) {
      this->compute_velocity(term.vel, Tr, u, aux, this->v_);
      this->compute_advection_shape(*el[term.row], Tr);
      mfem::AddMult_a_VWt(-weight, this->adv_, this->shape_[term.var],
                          *elmats(term.row, term.var));
    }
  }
}
