/**
 * @file TimeCHNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative for CH equations
 * @version 0.1
 * @date 2025-09-05
 *
 * Copyright CEA (C) 2025
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
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new TimeCHNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for a splitted time derivative operator.
 * It forwards the provided previous solution fields, simulation
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
TimeCHNLFormIntegrator<VARS>::TimeCHNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "SplitTimeDerivative";

  this->check_variables_consistency();
}

/**
 * @brief Initialize the SplitTimeDerivative integrator.
 *
 * This method performs all necessary setup steps for the SplitTimeDerivative
 * nonlinear form integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
 * 3. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' is populated with the
 *      required coefficient types for this integrator.
 *
 */
template <class VARS>
void TimeCHNLFormIntegrator<VARS>::init() {
  if (this->expected_list_.empty()) {
    this->get_coefficients();
  } else {
    this->check_coefficient_types(this->expected_list_);
    this->get_coefficients();
  }
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void TimeCHNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "TimeCHNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      break;
    }
  }
}

/**
 * @brief Retrieve and store the coefficients required by the TimeDerivative integrator.
 *
 * This method collects one default coefficient:
 * - 'coefficient_A' stores the first coefficient,
 *
 *
 * @remark By default, A is equal to one. This method can be overridden in derived classes to
 * provide custom behavior for retrieving coefficients.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 */
template <class VARS>
void TimeCHNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    this->coefficient_A.add(Coefficient(Glossary::Default, 1.0));
  }
}

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector  by element
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
void TimeCHNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  //////////////////////
  // Block 0 R(phi) on mu term
  {
    int blk = 0;

    int nd = el[blk]->GetDof();
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
  }

  //////////////////////
  // Block 1 R(mu) on dphi/dt term
  {
    int blk = 1;

    int off_blk = 0;

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
      const auto& phin = this->u_old_[blk].GetValue(Tr, ip);

      double coef_a =
          this->compute_coefficient(this->coefficient_A[blk], std::span<const double>({phi, phin}));
      const double ww = coef_a * phi * ip.weight * Tr.Weight();
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix by element.
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
void TimeCHNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("TimeCHNLFormIntegrator::AssembleElementGrad");
  // loop over diagonal entries
  // block 0 0  dR(phi)dphi
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
  }
  // block 0 1  dR(phi)dmu
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
  }

  // block 1 0   dR(mu)dPhi
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
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      el[blk]->CalcPhysShape(Tr, Psi);

      const auto& phi = *elfun[blk] * Psi;
      const auto& phin = this->u_old_[blk].GetValue(Tr, ip);

      double coef_a =
          this->compute_coefficient(this->coefficient_A[blk], std::span<const double>({phi, phin}));
      double w = coef_a * Tr.Weight() * ip.weight;
      AddMult_a_VVt(w, Psi, *elmats(blk, off_blk));
    }
  }
  // block 1 1 dR(mu)dmu
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
  }
}

/**
 * @brief Return the value of the coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 * @tparam VARS
 * @param coef
 * @param values
 * @return double
 */
template <class VARS>
double TimeCHNLFormIntegrator<VARS>::compute_coefficient(Coefficient coef,
                                                         const std::span<const double>& values) {
  std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
  std::span<const double> un(values.begin() + this->nb_blk_, values.end());
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    // coef_value = coef.compute({u});
    MFEM_VERIFY(false,
                "Implicit coefficient for TimeDerivative integrator not implemented yet. Please "
                "check your data.");
  } else if (coef.is_explicit()) {
    coef_value = coef.compute(un);
  } else if (coef.is_scalar()) {
    coef_value = coef.compute();
  }
  return coef_value;
}
