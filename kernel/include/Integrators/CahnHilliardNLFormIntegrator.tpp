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
#include <string>
#include <utility>
#include <vector>

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
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "CahnHilliard";
  this->check_variables_consistency();
}

/**
 * @brief Initialize the CahnHilliard integrator.
 *
 * This method performs all necessary setup steps for a CahnHilliard-type
 * nonlinear form integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
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
void CahnHilliardNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "CahnHilliardNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      temperature_found = true;
      break;
    }
  }
  if (this->params_.has_parameter("ScaleMobilityByTemperature")) {
    this->scale_mobility_by_temperature_ =
        this->params_.template get_param_value<bool>("ScaleMobilityByTemperature");
    if (this->scale_mobility_by_temperature_) {
      MFEM_VERIFY(
          temperature_found,
          "CahnHilliardNLFormIntegrator: "
          "Temperature variable required to scale mobility, but not found in auxiliary variables");
    }
  }
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
  std::vector<double> vaux_gf_at_ip;

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
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      const auto& mu = *elfun[off_blk] * Psi;
      const auto& phi = *elfun[blk] * Psi;
      const auto& phin = this->u_old_[blk].GetValue(Tr, ip);
      // Get aux values at ip TODO(cci) (move in method)
      vaux_gf_at_ip.clear();
      for (const auto& aux_gf : vaux_gf_) {
        vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
      }

      const double xx = ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);

      gradU *= -xx * this->compute_coefficient(lambda[blk], std::span<const double>({phi, phin}),
                                               std::span<const double>(vaux_gf_at_ip));
      gradPsi.AddMult(gradU, *elvect[blk]);

      // Given u, compute (w'(u), psi), psi is shape function
      const double ww =
          xx * (mu - this->compute_gradient_coefficient(double_well_energy[blk], blk,
                                                        std::span<const double>({phi, phin}),
                                                        std::span<const double>(vaux_gf_at_ip)));

      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
  //
  // Block 1 R(mu) = div M grad mu
  //
  {
    int blk = 1;
    int off_blk = 0;
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

    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      const auto& phi = *elfun[off_blk] * Psi;
      const auto& phin = this->u_old_[off_blk].GetValue(Tr, ip);
      vaux_gf_at_ip.clear();
      for (const auto& aux_gf : vaux_gf_) {
        vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= this->compute_coefficient(mobility[off_blk], std::span<const double>({phi, phin}),
                                         std::span<const double>(vaux_gf_at_ip)) *
               ip.weight * Tr.Weight();
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
  std::vector<double> vaux_gf_at_ip;

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
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      Tr.SetIntPoint(&ip);

      const auto& phi = *elfun[blk] * Psi;
      const auto& phin = this->u_old_[blk].GetValue(Tr, ip);

      // Get aux values at ip TODO(cci) (move in method)
      vaux_gf_at_ip.clear();
      for (const auto& aux_gf : vaux_gf_) {
        vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
      }

      const double xx = -ip.weight * Tr.Weight();

      el[blk]->CalcPhysDShape(Tr, gradPsi);

      const double coef_lambda =
          this->compute_coefficient(lambda[blk], std::span<const double>({phi, phin}),
                                    std::span<const double>(vaux_gf_at_ip));

      AddMult_a_AAt(xx * coef_lambda, gradPsi, *elmats(blk, blk));

      double fun_val =
          xx * this->compute_hessian_coefficient(double_well_energy[blk], blk, blk,
                                                 std::span<const double>({phi, phin}),
                                                 std::span<const double>(vaux_gf_at_ip));
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
      Tr.SetIntPoint(&ip);

      double ww = ip.weight * Tr.Weight();
      AddMult_a_VVt(ww, Psi, *elmats(blk, off_blk));
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

    elmats(blk, off_blk)->SetSize(nd);
    *elmats(blk, off_blk) = 0.0;
  }

  // Block 1 1  dR(mu)dmu=dR(mu)dphi=d(div M grad mu)/dmu
  {
    int blk = 1;
    int off_blk = 0;
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
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      vaux_gf_at_ip.clear();
      for (const auto& aux_gf : vaux_gf_) {
        vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
      }
      const auto& phi = *elfun[off_blk] * Psi;
      const auto& phin = this->u_old_[off_blk].GetValue(Tr, ip);

      const double coef_mob =
          this->compute_coefficient(mobility[blk], std::span<const double>({phi, phin}),
                                    std::span<const double>(vaux_gf_at_ip)) *
          ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);

      AddMult_a_AAt(coef_mob, gradPsi, *elmats(blk, blk));
    }
  }
}

/**
 * @brief Return the value of the coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_

 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).

 * @return The computed scalar value of the coefficient.
 */
template <class VARS>
double CahnHilliardNLFormIntegrator<VARS>::compute_coefficient(
    Coefficient coef, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  if (coef.is_scalar()) {
    return coef.compute();
  } else {
    std::span<const double> u(values.begin(), values.begin() + this->nb_blk_ - 1);
    std::span<const double> un(values.begin() + this->nb_blk_ - 1, values.end());

    const std::span<const double>& input = coef.is_implicit() ? u : un;
    return coef.compute(input, aux_values);
  }
}

/**
 * @brief Compute the value of a specific component of the gradient of a coefficient.
 *
 * This method evaluates the gradient of the given coefficient with respect to the
 * variable corresponding to the specified index 'blk'. By default, the 'values' vector
 * contains {u, u_old}, and auxiliary variables remain accessible via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose gradient is to be computed.
 * @param iblk    Index of the gradient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the gradient component of the coefficient.
 *
 * @note This function can be overridden in derived classes.
 */
template <class VARS>
double CahnHilliardNLFormIntegrator<VARS>::compute_gradient_coefficient(
    Coefficient coef, const int blk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  if (coef.is_scalar()) {
    return 0.0;
  } else {
    std::span<const double> u(values.begin(), values.begin() + this->nb_blk_ - 1);
    std::span<const double> un(values.begin() + this->nb_blk_ - 1, values.end());

    const std::span<const double>& input = coef.is_implicit() ? u : un;
    return coef.compute_gradient(blk, input, aux_values);
  }
}

/**
 * @brief Compute the value of a specific component of the Hessian of a  coefficient.
 *
 * This method evaluates the Hessian of the given coefficient with respect to the
 * variables corresponding to the specified indexes 'iblk' and 'jblk'. By default,
 * the 'values' vector contains {u, u_old}, and auxiliary variables remain accessible
 * via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose Hessian is to be computed.
 * @param iblk   Index of the row  of the Hessian component.
 * @param jblk   Index of the column of the Hessian component.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the Hessian component for the given coefficient.
 *
 * @note This function can be overridden in derived classes to implement
 *       more complex or nonlinear Hessian evaluations.
 */
template <class VARS>
double CahnHilliardNLFormIntegrator<VARS>::compute_hessian_coefficient(
    Coefficient coef, const int iblk, const int jblk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  if (coef.is_scalar()) {
    return 0.0;
  } else {
    std::span<const double> u(values.begin(), values.begin() + this->nb_blk_ - 1);
    std::span<const double> un(values.begin() + this->nb_blk_ - 1, values.end());

    const std::span<const double>& input = coef.is_implicit() ? u : un;
    return coef.compute_hessian(iblk, jblk, input, aux_values);
  }
}
