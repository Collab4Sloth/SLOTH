/**
 * @file SourceTermNLFormIntegrator.tpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @brief VF of a user-defined source term
 * @version 0.4
 * @date 2026-10-16
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
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new SourceTermNLFormIntegrator object.
 *
 * @param geometry     Geometry of the problem (Cartesian, axisymmetric, ...).
 * @param time_step    Current time step.
 * @param u_old        Vector of previous-time-step solution fields.
 * @param aux_old      Vector of previous-time-step auxiliary fields.
 * @param params       Parameters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the integrator.
 * @param coefficients List of coefficients defining the source terms.
 */
template <class VARS>
SourceTermNLFormIntegrator<VARS>::SourceTermNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, std::move(auxvars),
                                  coefficients) {
  this->integrator_name_ = "SourceTerm";
}

/**
 * @brief Initialize the integrator: check and store the coefficients, take ownership of the
 * auxiliary grid functions.
 */
template <class VARS>
void SourceTermNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  // reserve: avoids reallocations (and thus copies of ParGridFunction) while filling
  this->vaux_gf_.reserve(this->vaux_gf_.size() + this->aux_infos_.size());
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief Retrieve the SourceTerm coefficient (ID 0) of each block.
 *
 * Positional storage: source_term_[i] is the source term of block i. Also flags:
 *  - the source terms whose derivative can be non-zero (neither constant, explicit nor
 *    file-based), consistently with SlothNLFormIntegrator::eval_gradient_coefficient;
 *  - whether the unknowns and auxiliary fields must be evaluated (at least one non-constant
 *    source term).
 */
template <class VARS>
void SourceTermNLFormIntegrator<VARS>::get_coefficients() {
  this->has_derivative_.clear();
  this->any_derivative_ = false;
  this->needs_fields_ = false;
  for (unsigned int i = 0; i < this->nb_blk_; ++i) {
    auto coef = this->get_coefficient(i, GlossaryType::SourceTerm, 0);
    MFEM_VERIFY(coef.has_value(),
                "SourceTermNLFormIntegrator: a SourceTerm coefficient (ID 0) is required for "
                "every block (use a zero constant for blocks without source term).");
    const bool has_derivative = !coef->is_scalar() && !coef->is_explicit() && !coef->is_file();
    this->has_derivative_.push_back(has_derivative);
    this->any_derivative_ = this->any_derivative_ || has_derivative;
    this->needs_fields_ = this->needs_fields_ || !coef->is_scalar();
    this->source_term_.add(*coef);
  }
}

/**
 * @brief Integration rule shared by all blocks (exact for the highest-order block).
 */
template <class VARS>
const mfem::IntegrationRule& SourceTermNLFormIntegrator<VARS>::get_integration_rule(
    const mfem::Array<const mfem::FiniteElement*>& el, const mfem::ElementTransformation& Tr) {
  int order = 0;
  for (int b = 0; b < el.Size(); ++b) {
    order = std::max(order, el[b]->GetOrder());
  }
  return mfem::IntRules.Get(el[0]->GetGeomType(), 2 * order + Tr.OrderW());
}

/**
 * @brief Element-level setup: size the work arrays and, if needed, extract the local dofs of
 * u_old once (instead of one GridFunction::GetValue per block and integration point).
 */
template <class VARS>
void SourceTermNLFormIntegrator<VARS>::prepare_element(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr) {
  const int nb = el.Size();
  this->shape_.resize(nb);
  this->u_old_dofs_.resize(nb);
  this->u_values_.assign(2 * nb, 0.0);
  this->aux_values_.assign(this->vaux_gf_.size(), 0.0);

  for (int b = 0; b < nb; ++b) {
    this->shape_[b].SetSize(el[b]->GetDof());
    if (this->needs_fields_) {
      this->u_old_[b].GetElementDofValues(Tr.ElementNo, this->u_old_dofs_[b]);
    }
  }
}

/**
 * @brief Evaluate the shape functions at ip and, if needed, the current, previous and
 * auxiliary values (once for all blocks).
 *
 * @return The quadrature weight at ip (including the axisymmetric factor if needed).
 */
template <class VARS>
double SourceTermNLFormIntegrator<VARS>::evaluate_at_ip(
    const mfem::Array<const mfem::FiniteElement*>& el,
    const mfem::Array<const mfem::Vector*>& elfun, mfem::ElementTransformation& Tr,
    const mfem::IntegrationPoint& ip) {
  this->set_int_point(Tr, ip);

  const int nb = el.Size();
  for (int b = 0; b < nb; ++b) {
    el[b]->CalcShape(ip, this->shape_[b]);
  }

  if (this->needs_fields_) {
    for (int b = 0; b < nb; ++b) {
      this->u_values_[b] = (*elfun[b]) * this->shape_[b];
      this->u_values_[b + nb] = this->u_old_dofs_[b] * this->shape_[b];
    }
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
 * @brief Assemble the element-level residual vector: R_i = - (f_i, psi_i).
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elvect  Array of vectors where the element residual contributions are stored.
 */
template <class VARS>
void SourceTermNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  // Catch_Time_Section("SourceTermNLFormIntegrator::AssembleElementVector");
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

    for (int b = 0; b < nb; ++b) {
      const double s = this->eval_coefficient(this->source_term_[b], u, aux, b);
      elvect[b]->Add(-s * weight, this->shape_[b]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian: J_ij = - (df_i/du_j phi_j, psi_i).
 *
 * All blocks (including off-diagonal ones) are sized and reset, so that no stale contribution
 * left by another integrator in the shared elmats can be re-assembled. Rows whose source term
 * cannot depend on the current solution are skipped; if none can, no integration is performed.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elmats  Array of dense matrices where the element Jacobian contributions are stored.
 */
template <class VARS>
void SourceTermNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("SourceTermNLFormIntegrator::AssembleElementGrad");
  const int nb = el.Size();

  for (int i = 0; i < nb; ++i) {
    for (int j = 0; j < nb; ++j) {
      elmats(i, j)->SetSize(el[i]->GetDof(), el[j]->GetDof());
      *elmats(i, j) = 0.0;
    }
  }
  if (!this->any_derivative_) {
    return;  // constant / explicit / file-based source terms: zero Jacobian
  }

  this->prepare_element(el, Tr);

  const std::span<const double> u(this->u_values_);
  const std::span<const double> aux(this->aux_values_);
  const mfem::IntegrationRule& ir = get_integration_rule(el, Tr);

  for (int q = 0; q < ir.GetNPoints(); ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    const double weight = this->evaluate_at_ip(el, elfun, Tr, ip);

    for (int i = 0; i < nb; ++i) {
      if (!this->has_derivative_[i]) {
        continue;
      }
      for (int j = 0; j < nb; ++j) {
        const double ds = this->eval_gradient_coefficient(this->source_term_[i], j, u, aux, i);
        if (ds == 0.0) {
          continue;
        }
        mfem::AddMult_a_VWt(-ds * weight, this->shape_[i], this->shape_[j], *elmats(i, j));
      }
    }
  }
}
