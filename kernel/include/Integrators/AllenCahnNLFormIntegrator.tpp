/**
 * @file AllenCahnNLFormIntegrator.tpp
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
#include <typeinfo>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AllenCahnNLFormIntegrator object.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 *
 */
template <class VARS>
AllenCahnNLFormIntegrator<VARS>::AllenCahnNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {
  this->integrator_name_ = "AllenCahn";

  this->check_variables_consistency();
}

/**
 * @brief Initialize the AllenCahn integrator.
 *
 * 1. Checks that all coefficients contain the expected types.
 * 2. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::check_variables_consistency() {}

//////////////////////////////////////////////////////////////////////////////
// Residual
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * Blocks sharing the same finite element (blocks with different finite elements: see
 * AssembleElementVectorGeneric):
 * - integration points in the outer loop, blocks in the inner loop;
 * - shape functions, physical gradients and weight computed once per integration point;
 * - {u, u_old} and auxiliary values computed once per element for all the blocks
 *   (see SlothNLFormIntegrator::prepare_element_values);
 * - coefficients evaluated by reference (no copy).
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int num_blocks = el.Size();
  if (!this->same_finite_elements(el)) {
    this->AssembleElementVectorGeneric(el, Tr, elfun, elvect);
    return;
  }

  const mfem::FiniteElement& fe = *el[0];
  const int nd = fe.GetDof();
  const int dim = fe.GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  gradU.SetSize(dim);
  for (int blk = 0; blk < num_blocks; ++blk) {
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
  }

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
  this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, *ir);

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    fe.CalcShape(ip, Psi);
    this->set_int_point(Tr, ip);
    const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

    fe.CalcPhysDShape(Tr, gradPsi);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    for (int blk = 0; blk < num_blocks; ++blk) {
      const double coef_mobi = get_mob_at_ip(blk, ipv.u, ipv.aux);
      const double lamb = this->eval_coefficient(lambda[blk], ipv.u, ipv.aux, blk);

      // Laplacian : (grad(u), grad(psi))
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= coef_mobi * lamb * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // (w'(u), psi)
      double ww = coef_mobi * this->eval_gradient_energy_coefficient(double_well_energy[blk], blk,
                                                                     ipv.u, ipv.aux, blk);
      ww *= weight_coef;
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Jacobian
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * Blocks sharing the same finite element (blocks with different finite elements: see
 * AssembleElementGradGeneric). Organised as AssembleElementVector: everything that does not
 * depend on the block is computed once per integration point, and the mobility of each block
 * once per integration point.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int num_blocks = el.Size();
  if (!this->same_finite_elements(el)) {
    this->AssembleElementGradGeneric(el, Tr, elfun, elmats);
    return;
  }

  const mfem::FiniteElement& fe = *el[0];
  const int nd = fe.GetDof();
  const int dim = fe.GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  for (int blk = 0; blk < num_blocks; ++blk) {
    for (int jblk = 0; jblk < num_blocks; ++jblk) {
      elmats(blk, jblk)->SetSize(nd);
      *elmats(blk, jblk) = 0.0;
    }
  }

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
  this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, *ir);
  this->mob_at_ip_.assign(static_cast<std::size_t>(num_blocks), 0.);

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    fe.CalcShape(ip, Psi);
    this->set_int_point(Tr, ip);
    const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

    fe.CalcPhysDShape(Tr, gradPsi);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    for (int b = 0; b < num_blocks; ++b) {
      this->mob_at_ip_[static_cast<std::size_t>(b)] = get_mob_at_ip(b, ipv.u, ipv.aux);
    }

    for (int blk = 0; blk < num_blocks; ++blk) {
      // Diagonal block: Laplacian + w''(u)
      const double coef_mobi = this->mob_at_ip_[static_cast<std::size_t>(blk)];
      const double lamb = this->eval_coefficient(lambda[blk], ipv.u, ipv.aux, blk);
      const double mob_lamb = coef_mobi * lamb * weight_coef;
      AddMult_a_AAt(mob_lamb, gradPsi, *elmats(blk, blk));

      double fun_val = coef_mobi * this->eval_hessian_coefficient(double_well_energy[blk], blk, blk,
                                                                  ipv.u, ipv.aux, blk);
      fun_val *= weight_coef;
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));

      // Off-diagonal blocks: d/du_j (M_b dW_b/du_b) = M_b d2W_b/du_b du_j (mobility and free
      // energy of the row block blk; the derivative of M_b is added below)
      for (int jblk = 0; jblk < num_blocks; ++jblk) {
        if (jblk == blk) continue;
        double fun_val_off =
            coef_mobi *
            this->eval_hessian_coefficient(double_well_energy[blk], blk, jblk, ipv.u, ipv.aux, blk);
        fun_val_off *= weight_coef;
        AddMult_a_VVt(fun_val_off, Psi, *elmats(blk, jblk));
      }

      // Mobility depending on the unknowns:
      //   d/du_j (M (lambda grad u . grad psi_i + W' psi_i))
      //     += dM/du_j psi_j (lambda grad u . grad psi_i + W' psi_i)
      this->add_mobility_derivative(blk, lamb, gradPsi, Psi, *elfun[blk], ipv.u, ipv.aux,
                                    weight_coef, elmats);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// General case: blocks with different finite elements
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Integration rule for blocks with different finite elements: the one of the
 * highest-order block, so that the products of shape functions of any two blocks are
 * integrated as in the case of identical finite elements.
 */
template <class VARS>
const mfem::IntegrationRule& AllenCahnNLFormIntegrator<VARS>::generic_integration_rule(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr) const {
  int order = 0;
  for (int b = 0; b < el.Size(); ++b) {
    order = std::max(order, el[b]->GetOrder());
  }
  return mfem::IntRules.Get(el[0]->GetGeomType(), 2 * order + Tr.OrderW());
}

/**
 * @brief At an integration point, for blocks with different finite elements: shape functions
 * of each block (shapes_), value of each block computed with its own shape functions and
 * previous values (u_generic_ = {u, u_old}), auxiliary values (aux_generic_).
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::generic_values_at_ip(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::IntegrationPoint& ip, const mfem::Array<const mfem::Vector*>& elfun) {
  const int num_blocks = el.Size();
  this->set_int_point(Tr, ip);
  this->shapes_.resize(static_cast<std::size_t>(num_blocks));
  this->u_generic_.resize(2 * static_cast<std::size_t>(num_blocks));
  for (int b = 0; b < num_blocks; ++b) {
    auto& shape = this->shapes_[static_cast<std::size_t>(b)];
    shape.SetSize(el[b]->GetDof());
    el[b]->CalcShape(ip, shape);
    this->u_generic_[static_cast<std::size_t>(b)] = (*elfun[b]) * shape;
    this->u_generic_[static_cast<std::size_t>(b + num_blocks)] =
        this->u_old_[static_cast<std::size_t>(b)].GetValue(Tr, ip);
  }
  this->aux_generic_.resize(this->vaux_gf_.size());
  for (std::size_t k = 0; k < this->vaux_gf_.size(); ++k) {
    this->aux_generic_[k] = this->vaux_gf_[k].GetValue(Tr, ip);
  }
}

/**
 * @brief Element residual for blocks with different finite elements (same formulas as
 * AssembleElementVector).
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementVectorGeneric(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    elvect[blk]->SetSize(el[blk]->GetDof());
    *elvect[blk] = 0.;
  }

  const mfem::IntegrationRule& ir = this->generic_integration_rule(el, Tr);
  for (int i = 0; i < ir.GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(i);
    this->generic_values_at_ip(el, Tr, ip, elfun);
    const std::span<const double> u(this->u_generic_);
    const std::span<const double> aux(this->aux_generic_);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    for (int blk = 0; blk < num_blocks; ++blk) {
      gradPsi.SetSize(el[blk]->GetDof(), el[blk]->GetDim());
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradU.SetSize(el[blk]->GetDim());

      const double coef_mobi = get_mob_at_ip(blk, u, aux);
      const double lamb = this->eval_coefficient(lambda[blk], u, aux, blk);

      // Laplacian : (grad(u), grad(psi))
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= coef_mobi * lamb * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // (w'(u), psi)
      const double ww =
          coef_mobi * this->eval_gradient_energy_coefficient(double_well_energy[blk], blk, u, aux, blk) *
          weight_coef;
      elvect[blk]->Add(ww, this->shapes_[static_cast<std::size_t>(blk)]);
    }
  }
}

/**
 * @brief Element Jacobian for blocks with different finite elements (same formulas as
 * AssembleElementGrad). The block (blk, jblk) is a nd_blk x nd_jblk matrix: test functions
 * of block blk, trial functions of block jblk.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementGradGeneric(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    for (int jblk = 0; jblk < num_blocks; ++jblk) {
      elmats(blk, jblk)->SetSize(el[blk]->GetDof(), el[jblk]->GetDof());
      *elmats(blk, jblk) = 0.0;
    }
  }
  this->mob_at_ip_.assign(static_cast<std::size_t>(num_blocks), 0.);

  const mfem::IntegrationRule& ir = this->generic_integration_rule(el, Tr);
  for (int i = 0; i < ir.GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(i);
    this->generic_values_at_ip(el, Tr, ip, elfun);
    const std::span<const double> u(this->u_generic_);
    const std::span<const double> aux(this->aux_generic_);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    for (int b = 0; b < num_blocks; ++b) {
      this->mob_at_ip_[static_cast<std::size_t>(b)] = get_mob_at_ip(b, u, aux);
    }

    for (int blk = 0; blk < num_blocks; ++blk) {
      const mfem::Vector& Psi_blk = this->shapes_[static_cast<std::size_t>(blk)];
      gradPsi.SetSize(el[blk]->GetDof(), el[blk]->GetDim());
      el[blk]->CalcPhysDShape(Tr, gradPsi);

      // Diagonal block: Laplacian + w''(u)
      const double coef_mobi = this->mob_at_ip_[static_cast<std::size_t>(blk)];
      const double lamb = this->eval_coefficient(lambda[blk], u, aux, blk);
      AddMult_a_AAt(coef_mobi * lamb * weight_coef, gradPsi, *elmats(blk, blk));
      const double fun_val =
          coef_mobi *
          this->eval_hessian_coefficient(double_well_energy[blk], blk, blk, u, aux, blk) * weight_coef;
      AddMult_a_VVt(fun_val, Psi_blk, *elmats(blk, blk));

      // Off-diagonal blocks (as AssembleElementGrad): M_b d2W_b/du_b du_j
      for (int jblk = 0; jblk < num_blocks; ++jblk) {
        if (jblk == blk) continue;
        const double fun_val_off =
            coef_mobi * this->eval_hessian_coefficient(double_well_energy[blk], blk, jblk, u, aux, blk) *
            weight_coef;
        AddMult_a_VWt(fun_val_off, Psi_blk, this->shapes_[static_cast<std::size_t>(jblk)],
                      *elmats(blk, jblk));
      }

      // Mobility depending on the unknowns (trial functions of each block)
      this->add_mobility_derivative(blk, lamb, gradPsi, Psi_blk, *elfun[blk], u, aux, weight_coef,
                                    elmats, &this->shapes_);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Coefficients
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Retrieve and store the coefficients required by the Allen-Cahn integrator.
 *
 * This method collects the coefficients of type `Capillary`, `Mobility`, and `FreeEnergy`
 * from each coefficients and adds them to the corresponding internal storage:
 * - 'lambda' stores the capillary coefficients,
 * - 'mobility' stores the mobility coefficients,
 * - 'double_well_energy' stores the free energy coefficients.
 *
 * Only coefficients with ID 0 are considered for each block.
 *
 * @tparam VARS Template parameter defining the variabls used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Capillary, 0).has_value()) {
      lambda.add(*(this->get_coefficient(i, GlossaryType::Capillary, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::FreeEnergy, 0).has_value()) {
      double_well_energy.add(*(this->get_coefficient(i, GlossaryType::FreeEnergy, 0)));
    }
  }
}

/**
 * @brief Add the derivative of the mobility to the rows of block blk:
 *   elmats(blk, j) += dM/du_j psi_j (lambda grad u . grad psi_i + W' psi_i)
 * Nothing is added for a mobility independent of the unknowns (dM/du_j = 0).
 * psi_j: shape functions of block j, given by col_shapes (blocks with different finite
 * elements), Psi otherwise.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::add_mobility_derivative(
    const int blk, const double lamb, const mfem::DenseMatrix& gradPsi, const mfem::Vector& Psi,
    const mfem::Vector& elfun_blk, const std::span<const double>& u,
    const std::span<const double>& aux, const double weight_coef,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats, const std::vector<mfem::Vector>* col_shapes) {
  const int num_blocks = static_cast<int>(this->nb_blk_);
  bool has_derivative = false;
  this->dmob_at_ip_.assign(static_cast<std::size_t>(num_blocks), 0.);
  for (int j = 0; j < num_blocks; ++j) {
    const double d = this->eval_gradient_coefficient(mobility[blk], j, u, aux, blk);
    this->dmob_at_ip_[static_cast<std::size_t>(j)] = d;
    has_derivative = has_derivative || (d != 0.0);
  }
  if (!has_derivative) return;

  // g_i = lambda grad u . grad psi_i + W' psi_i
  mfem::Vector grad_u(gradPsi.Width()), g(gradPsi.Height());
  gradPsi.MultTranspose(elfun_blk, grad_u);
  gradPsi.Mult(grad_u, g);
  g *= lamb;
  g.Add(this->eval_gradient_energy_coefficient(double_well_energy[blk], blk, u, aux, blk), Psi);

  for (int j = 0; j < num_blocks; ++j) {
    const double d = this->dmob_at_ip_[static_cast<std::size_t>(j)];
    if (d == 0.0) continue;
    const mfem::Vector& Psi_j =
        col_shapes == nullptr ? Psi : (*col_shapes)[static_cast<std::size_t>(j)];
    AddMult_a_VWt(d * weight_coef, g, Psi_j, *elmats(blk, j));
  }
}

/**
 * @brief Compute the value of the mobility at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param blk   Index of the block.
 * @param values {u, u_old} at the integration point
 * @param aux_values auxiliary values at the integration point
 *
 * @return The computed mobility at integration point
 */
template <class VARS>
double AllenCahnNLFormIntegrator<VARS>::get_mob_at_ip(
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  return this->eval_coefficient(mobility[blk], values, aux_values, static_cast<int>(blk));
}
