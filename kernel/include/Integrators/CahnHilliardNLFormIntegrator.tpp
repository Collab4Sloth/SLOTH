/**
 * @file CahnHilliardNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the CahnHilliard equations
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
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new CahnHilliardNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for Cahn-Hilliard-type
 * problems. It forwards the provided previous solution fields, simulation
 * parameters, auxiliary variables, and coefficients to the base
 * SLOTH nonlinear form integrator.
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
CahnHilliardNLFormIntegrator<VARS>::CahnHilliardNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {
  this->integrator_name_ = "CahnHilliard";
  this->check_variables_consistency();
}

/**
 * @brief Initialize the CahnHilliard integrator.
 *
 * This method performs all necessary setup steps for a CahnHilliard-type
 * nonlinear form integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types: Mobility, FreeEnergy and
 *    Capillary, or Mobility and FreeEnergy if the FreeEnergy of block 0 is variational (no
 *    Capillary then: the gradient terms are part of the free energy).
 * 3. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' is populated with the
 *      required coefficient types for this CahnHilliard integrator.
 *
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::init() {
  const auto energy = this->coefficients_.empty()
                          ? std::nullopt
                          : this->get_coefficient(0, GlossaryType::FreeEnergy, 0);
  this->variational_ = energy.has_value() && energy->is_variational();
  if (this->variational_) {
    this->check_coefficient_types(this->expected_variational_list_);
    for (unsigned int i = 0; i < this->nb_blk_; i++) {
      MFEM_VERIFY(!this->get_coefficient(i, GlossaryType::Capillary, 0).has_value(),
                  "Error in " + this->integrator_name_ +
                      ": the variational FreeEnergy coefficient contains the gradient terms: a "
                      "Capillary coefficient must not be given.");
    }
  } else {
    this->check_coefficient_types(this->expected_list_);
  }
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
void CahnHilliardNLFormIntegrator<VARS>::check_variables_consistency() {
  // TODO(cci): order of variables?
}

/**
 * @brief Retrieve and store the coefficients required by the Cahn-Hilliard integrator.
 *
 * This method collects the coefficients of type `Capillary`, `Mobility`, and `FreeEnergy`
 * from each coefficients and adds them to the corresponding internal storage:
 * - 'lambda' stores the capillary coefficients,
 * - 'mobility' stores the mobility coefficients,
 * - 'double_well_energy' stores the free energy coefficients.
 *
 * Only coefficients with ID 0 are considered for each block.
 *
 * @remark This method can be overridden in derived classes to provide
 *         custom behavior for retrieving coefficients.
 *
 * @tparam VARS Template parameter defining the variabls used in the integrator.
 *
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::get_coefficients() {
  if (this->variational_) {
    // Mobility of each block and variational free energy of phi (block 0)
    for (unsigned int i = 0; i < this->nb_blk_; i++) {
      if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
        mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
      }
    }
    this->free_energy = this->get_coefficient(0, GlossaryType::FreeEnergy, 0);
    this->isotropic_ = this->free_energy->is_isotropic() && !this->free_energy->is_explicit();
    return;
  }
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
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector of the CahnHilliard-type nonlinear problem by element
 *
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elvect  Array of vectors where the computed element residual contributions
 *                will be stored.
 *
 * @note Users typically do not call this function directly; it is invoked
 *       internally during the assembly of the global nonlinear form.
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  if (this->variational_) {
    this->variational_element_vector(el, Tr, elfun, elvect);
    return;
  }
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
  this->prepare_theta_dofs(Tr, elfun, 1);

  //
  // Block 0 R(phi)=mu - w' + div lambda grad phi
  //
  {
    int blk = 0;
    int off_blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // Coefficients at all the integration points of the element
    this->eval_coefficient_at_points(lambda[blk], 0, SlothNLBase::CoefValue, *el[blk], Tr, *ir,
                                     elfun, this->vaux_gf_, this->lambda_at_points_);
    this->eval_coefficient_at_points(double_well_energy[blk], 0, SlothNLBase::CoefFirstDerivative,
                                     *el[blk], Tr, *ir, elfun, this->vaux_gf_,
                                     this->energy_at_points_);
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      const auto& mu = *elfun[off_blk] * Psi;

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      // phi = theta phi + (1 - theta) phi_old in the gradient term (theta-scheme)
      gradPsi.MultTranspose(this->gradient_term_dofs(elfun, blk), gradU);

      gradU *= -weight_coef * this->lambda_at_points_.value(i);
      gradPsi.AddMult(gradU, *elvect[blk]);

      // Given u, compute (w'(u), psi), psi is shape function
      const double ww = weight_coef * (mu - this->energy_at_points_.energy_derivative(i, blk));

      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
  //
  // Block 1 R(mu) = div M grad mu
  //
  {
    int blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    this->eval_coefficient_at_points(mobility[blk], 0, SlothNLBase::CoefValue, *el[blk], Tr, *ir,
                                     elfun, this->vaux_gf_, this->mobility_at_points_);

    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= this->mobility_at_points_.value(i) * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix of a CahHilliard-type nonlinear problem by element.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elmat   Array of dense matrices where the computed element Jacobian contributions
 *                will be stored.
 *
 * @note Users typically do not call this function directly; it is invoked
 *       internally during the assembly of the global nonlinear form.
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  if (this->variational_) {
    this->variational_element_grad(el, Tr, elfun, elmats);
    return;
  }
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
  this->prepare_theta_dofs(Tr, elfun, 1);

  // Block 0  0 dR(phi)dphi = d(mu - w' + div lambda grad phi)/dphi
  {
    int blk = 0;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // Coefficients at all the integration points of the element
    this->eval_coefficient_at_points(lambda[blk], 0, SlothNLBase::CoefValue, *el[blk], Tr, *ir,
                                     elfun, this->vaux_gf_, this->lambda_at_points_);
    this->eval_coefficient_at_points(double_well_energy[blk], 0, SlothNLBase::CoefSecondDerivative,
                                     *el[blk], Tr, *ir, elfun, this->vaux_gf_,
                                     this->energy_at_points_);
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);

      // Gradient term: weight theta of the theta-scheme
      const double coef_lambda = this->lambda_at_points_.value(i) * this->theta_;

      AddMult_a_AAt(-weight_coef * coef_lambda, gradPsi, *elmats(blk, blk));

      double fun_val = -weight_coef * this->energy_at_points_.second_derivative(i, blk, blk);
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));
    }
  }

  // Block 0 1 dR(phi)dmu=d(mu - w' + div lambda grad phi)/dmu
  {
    int blk = 0;
    int off_blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, off_blk)->SetSize(nd);
    *elmats(blk, off_blk) = 0.0;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      AddMult_a_VVt(weight_coef, Psi, *elmats(blk, off_blk));
    }
  }
  // Block 1 0  dR(mu)dphi=d(div M grad mu)/dphi
  {
    int blk = 1;
    int off_blk = 0;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, off_blk)->SetSize(nd, el[off_blk]->GetDof());
    *elmats(blk, off_blk) = 0.0;

    // Mobility depending on phi:
    //   dR(mu)_i/dphi_j = (dM/dphi) psi_j grad(mu) . grad(psi_i)
    // (zero, and skipped, for a mobility independent of phi)
    mfem::Vector Psi_phi(el[off_blk]->GetDof());
    mfem::Vector grad_mu(dim), adv(nd);
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // Mobility and its derivatives at all the integration points (also used by block 1 1)
    this->eval_coefficient_at_points(
        mobility[blk], 0, SlothNLBase::CoefValue | SlothNLBase::CoefFirstDerivative, *el[blk], Tr,
        *ir, elfun, this->vaux_gf_, this->mobility_at_points_);
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      const double dmob_dphi = this->mobility_at_points_.derivative(i, off_blk) * weight_coef;
      if (dmob_dphi == 0.0) continue;

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], grad_mu);
      gradPsi.Mult(grad_mu, adv);  // adv_i = grad(mu) . grad(psi_i)
      el[off_blk]->CalcShape(ip, Psi_phi);
      AddMult_a_VWt(dmob_dphi, adv, Psi_phi, *elmats(blk, off_blk));
    }
  }

  // Block 1 1  dR(mu)dmu=dR(mu)dphi=d(div M grad mu)/dmu
  {
    int blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // Mobility at the integration points: computed for block 1 0 (same element and rule)
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      const double coef_mob = this->mobility_at_points_.value(i) * weight_coef;
      el[blk]->CalcPhysDShape(Tr, gradPsi);

      AddMult_a_AAt(coef_mob, gradPsi, *elmats(blk, blk));
    }
  }
}

/**
 * @brief Theta-scheme: dofs of theta c + (1 - theta) c_old on the current element for the n
 * concentration blocks (nothing for theta = 1)
 *
 * The theta-scheme applies to the gradient term of the equation of mu only, for a classic or an
 * isotropic variational free energy. The chemical potentials are not interpolated.
 *
 * @param Tr Element transformation
 * @param elfun Values of the dofs of the blocks
 * @param n Number of concentration blocks (1 for the classic form: phi)
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::prepare_theta_dofs(
    mfem::ElementTransformation& Tr, const mfem::Array<const mfem::Vector*>& elfun, const int n) {
  if (this->theta_ == 1.) return;
  MFEM_VERIFY(!this->variational_ || this->isotropic_,
              "Error in " + this->integrator_name_ +
                  ": the theta-scheme requires an isotropic gradient energy.");
  MFEM_VERIFY(this->u_old_.size() >= static_cast<std::size_t>(n),
              "Error in " + this->integrator_name_ +
                  ": the theta-scheme requires the values at the beginning of the time step.");
  this->theta_dofs_.resize(static_cast<std::size_t>(n));
  for (int blk = 0; blk < n; ++blk) {
    const auto b = static_cast<std::size_t>(blk);
    const mfem::Vector& c = *elfun[blk];
    this->u_old_[b].GetElementDofValues(Tr.ElementNo, this->uold_dofs_);
    MFEM_VERIFY(this->uold_dofs_.Size() == c.Size(),
                "Error in " + this->integrator_name_ +
                    ": c and c_old must be on the same finite element space.");
    this->theta_dofs_[b].SetSize(c.Size());
    add(this->theta_, c, 1. - this->theta_, this->uold_dofs_, this->theta_dofs_[b]);
  }
}

/**
 * @brief Dofs used in the gradient term of the concentration block blk: c (theta = 1),
 * theta c + (1 - theta) c_old otherwise (see prepare_theta_dofs)
 */
template <class VARS>
const mfem::Vector& CahnHilliardNLFormIntegrator<VARS>::gradient_term_dofs(
    const mfem::Array<const mfem::Vector*>& elfun, const int blk) const {
  return this->theta_ == 1. ? *elfun[blk] : this->theta_dofs_[static_cast<std::size_t>(blk)];
}

//////////////////////////////////////////////////////////////////////////////
// Variational free energy: n concentrations c_i (blocks 0..n-1) and n potentials mu_i (blocks
// n..2n-1)
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Number of concentrations n of the variational form: blocks c_0..c_{n-1}, then
 * mu_0..mu_{n-1}
 */
template <class VARS>
int CahnHilliardNLFormIntegrator<VARS>::nb_concentrations(const int nb_blocks) const {
  MFEM_VERIFY(nb_blocks % 2 == 0 && nb_blocks > 0,
              "Error in " + this->integrator_name_ +
                  ": the blocks are n concentrations followed by n chemical potentials (" +
                  std::to_string(nb_blocks) + " blocks).");
  return nb_blocks / 2;
}

/**
 * @brief Evaluates the mobility of each chemical potential mu_i (block n + i) at all the
 * integration points of the element, a coefficient of one variable being a function of c_i.
 * Jacobian: derivatives with respect to the blocks.
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::potential_mobility_at_points(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const int n, const bool jacobian) {
  this->potential_mobility_.resize(static_cast<std::size_t>(n));
  const unsigned int outputs = jacobian
                                   ? (SlothNLBase::CoefValue | SlothNLBase::CoefFirstDerivative)
                                   : SlothNLBase::CoefValue;
  for (int i = 0; i < n; ++i) {
    this->eval_coefficient_at_points(mobility[n + i], i, outputs, fe, Tr, ir, elfun, this->vaux_gf_,
                                     this->potential_mobility_[static_cast<std::size_t>(i)]);
  }
}

/**
 * @brief Element residual with a variational free energy f(c, grad c) of the n concentrations:
 *   R_{c_i}(psi)  = int mu_i psi - df/dc_i psi - df/d(grad c_i) . grad psi,
 *   R_{mu_i}(psi) = int M_i grad mu_i . grad psi.
 * Isotropic f (kappa_i/2 |grad c_i|^2 + g(c)): operations of the classic form.
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::variational_element_vector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int n = this->nb_concentrations(el.Size());
  MFEM_VERIFY(this->same_finite_elements(el), "Error in " + this->integrator_name_ +
                                                  ": all the blocks must be on the same finite "
                                                  "element space.");
  const mfem::FiniteElement& fe = *el[0];
  const bool generic = !this->isotropic_;
  const int nd = fe.GetDof();
  const int dim = generic ? Tr.GetSpaceDim() : fe.GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  gradU.SetSize(dim);
  for (int blk = 0; blk < el.Size(); ++blk) {
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
  }
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
  this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, *ir);
  this->prepare_theta_dofs(Tr, elfun, n);
  this->variational_at_points(*this->free_energy, n, this->isotropic_, fe, Tr, *ir, elfun,
                              this->vaux_gf_, false, this->energy_);
  this->potential_mobility_at_points(fe, Tr, *ir, elfun, n, false);
  const double* kappa = generic ? nullptr : this->energy_.kappa.HostRead();
  const double* df_dvar = this->energy_.df_dvar.HostRead();

  // Equation of mu_i (block c_i): mu_i - df/dc_i + div df/d(grad c_i)
  for (int blk = 0; blk < n; ++blk) {
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      fe.CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      const auto& mu = *elfun[n + blk] * Psi;

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      fe.CalcPhysDShape(Tr, gradPsi);
      if (generic) {
        // (mu_i, psi) - (df/dc_i, psi) - (df/d(grad c_i), grad psi)
        elvect[blk]->Add(weight_coef * mu, Psi);
        this->add_variational_vector(this->energy_, i, blk, -weight_coef, Psi, gradPsi,
                                     *elvect[blk]);
        continue;
      }
      gradPsi.MultTranspose(this->gradient_term_dofs(elfun, blk), gradU);
      gradU *= -weight_coef * kappa[i * n + blk];
      gradPsi.AddMult(gradU, *elvect[blk]);

      const double ww = weight_coef * (mu - df_dvar[i * n + blk]);
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }

  // Equation of c_i (block mu_i): div M_i grad mu_i
  for (int c = 0; c < n; ++c) {
    const int blk = n + c;
    const CoefficientAtPoints& mob = this->potential_mobility_[static_cast<std::size_t>(c)];
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      fe.CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      fe.CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= mob.value(i) * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);
    }
  }
}

/**
 * @brief Element Jacobian with a variational free energy f(c, grad c) of the n concentrations:
 *   (c_i, c_j):   - second derivatives of int f (see add_variational_grad), none for an explicit f
 *   (c_i, mu_i):  int psi_j psi
 *   (mu_i, u_j):  int dM_i/du_j psi_j grad mu_i . grad psi   (all the blocks u_j)
 *   (mu_i, mu_i): int M_i grad psi_j . grad psi
 * Isotropic f: operations of the classic form.
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::variational_element_grad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int nb = el.Size();
  const int n = this->nb_concentrations(nb);
  MFEM_VERIFY(this->same_finite_elements(el), "Error in " + this->integrator_name_ +
                                                  ": all the blocks must be on the same finite "
                                                  "element space.");
  const mfem::FiniteElement& fe = *el[0];
  const bool generic = !this->isotropic_;
  const int nd = fe.GetDof();
  const int dim = generic ? Tr.GetSpaceDim() : fe.GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  for (int blk = 0; blk < nb; ++blk) {
    for (int jblk = 0; jblk < nb; ++jblk) {
      elmats(blk, jblk)->SetSize(nd);
      *elmats(blk, jblk) = 0.0;
    }
  }
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
  this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, *ir);
  this->prepare_theta_dofs(Tr, elfun, n);
  this->variational_at_points(*this->free_energy, n, this->isotropic_, fe, Tr, *ir, elfun,
                              this->vaux_gf_, true, this->energy_);
  this->potential_mobility_at_points(fe, Tr, *ir, elfun, n, true);
  const bool second = this->energy_.with_second_derivatives;
  const double* kappa = generic ? nullptr : this->energy_.kappa.HostRead();
  const double* d2f_dvar2 = second && !generic ? this->energy_.d2f_dvar2.HostRead() : nullptr;

  // Blocks (c_i, c_j): free energy
  for (int blk = 0; blk < n; ++blk) {
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      fe.CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      fe.CalcPhysDShape(Tr, gradPsi);
      if (generic) {
        for (int jblk = 0; jblk < n && second; ++jblk) {
          this->add_variational_grad(this->energy_, i, blk, jblk, -weight_coef, Psi, gradPsi,
                                     *elmats(blk, jblk));
        }
        continue;
      }

      const double coef_lambda = kappa[i * n + blk] * this->theta_;
      AddMult_a_AAt(-weight_coef * coef_lambda, gradPsi, *elmats(blk, blk));

      for (int jblk = 0; jblk < n && d2f_dvar2 != nullptr; ++jblk) {
        double fun_val = -weight_coef * d2f_dvar2[(i * n + blk) * n + jblk];
        AddMult_a_VVt(fun_val, Psi, *elmats(blk, jblk));
      }
    }
  }

  // Blocks (c_i, mu_i): mass
  for (int blk = 0; blk < n; ++blk) {
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      fe.CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      AddMult_a_VVt(weight_coef, Psi, *elmats(blk, n + blk));
    }
  }

  // Blocks (mu_i, u_j): mobility depending on the unknowns (in practice the concentrations)
  //   dR(mu_i)/du_j = (dM_i/du_j) psi_j grad(mu_i) . grad(psi)
  this->grad_mu_.SetSize(dim);
  this->adv_.SetSize(nd);
  for (int c = 0; c < n; ++c) {
    const int blk = n + c;
    const CoefficientAtPoints& mob = this->potential_mobility_[static_cast<std::size_t>(c)];
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      fe.CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      bool computed = false;
      for (int jblk = 0; jblk < nb; ++jblk) {
        const double dmob = mob.derivative(i, jblk) * weight_coef;
        if (dmob == 0.0) continue;
        if (!computed) {
          fe.CalcPhysDShape(Tr, gradPsi);
          gradPsi.MultTranspose(*elfun[blk], this->grad_mu_);
          gradPsi.Mult(this->grad_mu_, this->adv_);  // adv_k = grad(mu_i) . grad(psi_k)
          computed = true;
        }
        AddMult_a_VWt(dmob, this->adv_, Psi, *elmats(blk, jblk));
      }
    }
  }

  // Blocks (mu_i, mu_i): mobility
  for (int c = 0; c < n; ++c) {
    const int blk = n + c;
    const CoefficientAtPoints& mob = this->potential_mobility_[static_cast<std::size_t>(c)];
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      fe.CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      const double coef_mob = mob.value(i) * weight_coef;
      fe.CalcPhysDShape(Tr, gradPsi);

      AddMult_a_AAt(coef_mob, gradPsi, *elmats(blk, blk));
    }
  }
}
