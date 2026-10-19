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
 * Same computation as the historical implementation (AssembleElementVectorGeneric), reorganised:
 * - integration points in the outer loop, blocks in the inner loop;
 * - shape functions, physical gradients and weight computed once per integration point;
 * - {u, u_old} and auxiliary values computed once per element for all the blocks
 *   (see SlothNLFormIntegrator::prepare_element_values);
 * - coefficients evaluated by reference (no copy).
 * The contributions to each elvect[blk] are accumulated in the same order as before.
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
      const double lamb = this->eval_coefficient(lambda[blk], ipv.u, ipv.aux);

      // Laplacian : (grad(u), grad(psi))
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= coef_mobi * lamb * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // (w'(u), psi)
      double ww = coef_mobi * this->eval_gradient_energy_coefficient(double_well_energy[blk], blk,
                                                                     ipv.u, ipv.aux);
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
 * Same computation as the historical implementation (AssembleElementGradGeneric), reorganised
 * as AssembleElementVector: everything that does not depend on the block is computed once per
 * integration point, and the mobility of each block once per integration point. The off-diagonal
 * blocks no longer re-run the whole quadrature loop for each (blk, jblk) pair.
 * The contributions to each elmats(blk, jblk) are accumulated in the same order as before.
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
      const double lamb = this->eval_coefficient(lambda[blk], ipv.u, ipv.aux);
      const double mob_lamb = coef_mobi * lamb * weight_coef;
      AddMult_a_AAt(mob_lamb, gradPsi, *elmats(blk, blk));

      double fun_val = coef_mobi * this->eval_hessian_coefficient(double_well_energy[blk], blk, blk,
                                                                  ipv.u, ipv.aux);
      fun_val *= weight_coef;
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));

      // Off-diagonal blocks (same formula as the historical implementation:
      // mobility and free energy of block jblk, Hessian component (blk, jblk))
      for (int jblk = 0; jblk < num_blocks; ++jblk) {
        if (jblk == blk) continue;
        double fun_val_off =
            this->mob_at_ip_[static_cast<std::size_t>(jblk)] *
            this->eval_hessian_coefficient(double_well_energy[jblk], blk, jblk, ipv.u, ipv.aux);
        fun_val_off *= weight_coef;
        AddMult_a_VVt(fun_val_off, Psi, *elmats(blk, jblk));
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Historical implementations (blocks with different finite elements)
//////////////////////////////////////////////////////////////////////////////

template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementVectorGeneric(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);
  std::vector<double> vaux_gf_at_ip(this->vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      for (size_t k = 0; k < this->vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = this->vaux_gf_[k].GetValue(Tr, ip);
      }
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);

      const double coef_mobi = get_mob_at_ip(blk, std::span<const double>(u_values),
                                             std::span<const double>(vaux_gf_at_ip));
      const double lamb = this->compute_coefficient(lambda[blk], std::span<const double>(u_values),
                                                    std::span<const double>(vaux_gf_at_ip));

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      gradU *= coef_mobi * lamb * weight_coef;

      gradPsi.AddMult(gradU, *elvect[blk]);

      double ww = coef_mobi * this->compute_gradient_energy_coefficient(
                                  double_well_energy[blk], blk, std::span<const double>(u_values),
                                  std::span<const double>(vaux_gf_at_ip));
      ww *= weight_coef;
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementGradGeneric(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(this->vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    // Diag
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      for (size_t k = 0; k < this->vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = this->vaux_gf_[k].GetValue(Tr, ip);
      }
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      const double coef_mobi = get_mob_at_ip(blk, std::span<const double>(u_values),
                                             std::span<const double>(vaux_gf_at_ip));

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      const double lamb = this->compute_coefficient(lambda[blk], std::span<const double>(u_values),
                                                    std::span<const double>(vaux_gf_at_ip));

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      const double mob_lamb = coef_mobi * lamb * weight_coef;
      AddMult_a_AAt(mob_lamb, gradPsi, *elmats(blk, blk));

      double fun_val =
          coef_mobi * this->compute_hessian_coefficient(double_well_energy[blk], blk, blk,
                                                        std::span<const double>(u_values),
                                                        std::span<const double>(vaux_gf_at_ip));
      fun_val *= weight_coef;

      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));
    }

    // Off-diag
    for (int jblk = 0; jblk < num_blocks; ++jblk) {
      if (jblk == blk) continue;
      elmats(blk, jblk)->SetSize(nd);
      *elmats(blk, jblk) = 0.0;
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const mfem::IntegrationPoint& ip = ir->IntPoint(i);
        el[blk]->CalcShape(ip, Psi);
        this->set_int_point(Tr, ip);

        for (size_t k = 0; k < this->vaux_gf_.size(); ++k) {
          vaux_gf_at_ip[k] = this->vaux_gf_[k].GetValue(Tr, ip);
        }
        for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
          u_values[off_blk] = (*elfun[off_blk]) * Psi;
          u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
        }

        const double coef_mobi = get_mob_at_ip(jblk, std::span<const double>(u_values),
                                               std::span<const double>(vaux_gf_at_ip));

        double fun_val =
            coef_mobi * this->compute_hessian_coefficient(double_well_energy[jblk], blk, jblk,
                                                          std::span<const double>(u_values),
                                                          std::span<const double>(vaux_gf_at_ip));
        double weight_coef = ip.weight * Tr.Weight();
        if (this->isAxisymmetric()) {
          weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
        }

        fun_val *= weight_coef;
        AddMult_a_VVt(fun_val, Psi, *elmats(blk, jblk));
      }
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
  return this->eval_coefficient(mobility[blk], values, aux_values);
}
