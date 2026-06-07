/**
 * @file RobinNLFormIntegrator.tpp
 * @author Marine Harel (marine.harel@cea.fr)
 * @brief Robin integrator
 * @version 0.1
 * @date 2026-02-23
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

#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Integrators/RobinNLFormIntegrator.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Retrieve and store the coefficients required by the integrator
 *        or raises an error if a coefficient is missing.
 *
 * This method collects coefficients of type RobinA and RobinB and stores them:
 * - 'robin_a' stores the a-coefficient for the Robin boundary condition,
 * - 'robin_b' stores the b-coefficient for the Robin boundary condition.
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
void RobinNLFormIntegrator<VARS>::get_coefficients() {
  auto robin_a_opt = this->get_coefficient(blk_, GlossaryType::RobinA, 0, this->bdr_id_);
  auto robin_b_opt = this->get_coefficient(blk_, GlossaryType::RobinB, 0, this->bdr_id_);

  if (robin_a_opt.has_value() && robin_b_opt.has_value()) {
    this->robin_a = *robin_a_opt;
    this->robin_b = *robin_b_opt;
  } else {
    std::string msg =
        "RobinNLFormIntegrator<VARS>::get_coefficients(): "
        "At least one Robin coefficient not found for boundary " +
        std::to_string(this->bdr_id_);
    mfem::mfem_error(msg.c_str());
  }
}

/**
 * @brief Construct a new RobinNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base SLOTH
 * nonlinear form integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 * @param block        Block to which apply Robin condition
 * @param bdr_id       Boundary id
 *
 */
template <class VARS>
RobinNLFormIntegrator<VARS>::RobinNLFormIntegrator(
    Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients,
    const unsigned int block, const unsigned int bdr_id)
    : SlothNLFormIntegrator<VARS>(geometry, u_old, aux_old, params, auxvars, coefficients),
      blk_(block),
      bdr_id_(bdr_id),
      robin_a(Coefficient(Glossary::Default, 0.0)),
      robin_b(Coefficient(Glossary::Default, 0.0)) {
  this->integrator_name_ = "Robin";
}

/**
 * @brief Initialize the integrator.
 *
 * This method retrieves and stores the coefficients, or raises an error if they were not provided.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void RobinNLFormIntegrator<VARS>::init() {
  this->get_coefficients();
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
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
void RobinNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);
  std::vector<double> vaux_gf_at_ip;
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Set all blocks to zero
    int nd = el[blk]->GetDof();
    Psi.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    // Add Robin on the right block
    if (blk == static_cast<int>(blk_)) {
      const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

      // Loop over integration points
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const mfem::IntegrationPoint& ip = ir->IntPoint(i);

        el[blk]->CalcShape(ip, Psi);
        Tr.SetIntPoint(&ip);

        // Get aux values at ip TODO(cci) (move in method)
        vaux_gf_at_ip.clear();
        for (const auto& aux_gf : vaux_gf_) {
          vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
        }
        // Get values
        for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
          u_values[off_blk] = (*elfun[off_blk]) * Psi;
          u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
        }

        const double val =
            Tr.Weight() * (this->compute_coefficient(robin_a, std::span<const double>(u_values),
                                                     std::span<const double>(vaux_gf_at_ip)) *
                               u_values[blk] -
                           this->compute_coefficient(robin_b, std::span<const double>(u_values),
                                                     std::span<const double>(vaux_gf_at_ip)));

        add(*elvect[blk], ip.weight * val, Psi, *elvect[blk]);
      }
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
void RobinNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array2D<mfem::DenseMatrix*>& elmat) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);
  std::vector<double> vaux_gf_at_ip;
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Set all blocks to zero
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    Psi.SetSize(nd);
    gradPsi.SetSize(nd, dim);
    elmat(blk, blk)->SetSize(nd);
    *elmat(blk, blk) = 0.0;

    // Add Robin on the right block
    if (blk == static_cast<int>(blk_)) {
      const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

      for (int i = 0; i < ir->GetNPoints(); i++) {
        const mfem::IntegrationPoint& ip = ir->IntPoint(i);
        el[blk]->CalcShape(ip, Psi);

        // Get aux values at ip TODO(cci) (move in method)
        vaux_gf_at_ip.clear();
        for (const auto& aux_gf : vaux_gf_) {
          vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
        }
        // Get values
        for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
          u_values[off_blk] = (*elfun[off_blk]) * Psi;
          u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
        }

        Tr.SetIntPoint(&ip);
        double weight_coef = ip.weight * Tr.Weight();
        if (this->isAxisymmetric()) {
          weight_coef *= mfem::CylindricalRadialCoefficient().Eval(Tr, ip);
        }

        const double coeff_robin =
            (this->compute_coefficient(robin_a, std::span<const double>(u_values),
                                       std::span<const double>(vaux_gf_at_ip)) +
             this->compute_gradient_coefficient(robin_a, blk, std::span<const double>(u_values),
                                                std::span<const double>(vaux_gf_at_ip)) *
                 u_values[blk]) *
            weight_coef;
        AddMult_a_VVt(coeff_robin, Psi, *elmat(blk, blk));
      }
    }
  }
}
