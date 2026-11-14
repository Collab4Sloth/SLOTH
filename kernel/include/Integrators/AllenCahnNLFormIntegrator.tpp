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
 * All the blocks share the same finite element:
 * - integration points in the outer loop, blocks in the inner loop;
 * - shape functions, physical gradients and weight computed once per integration point;
 * - {u, u_old} and auxiliary values computed once per element for all the blocks
 *   (see SlothNLFormIntegrator::prepare_element_values);
 * - coefficients evaluated once per element at all the integration points (see
 *   SlothNLFormIntegrator::eval_coefficient_at_points): by batch for the coefficients generated
 *   from a JSON file.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int num_blocks = el.Size();
  MFEM_VERIFY(this->same_finite_elements(el),
              "Error in AllenCahn: all the phases must be on the same finite element space.");

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
  this->coefficients_at_points(fe, Tr, *ir, elfun, false);

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    fe.CalcShape(ip, Psi);
    this->set_int_point(Tr, ip);

    fe.CalcPhysDShape(Tr, gradPsi);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    for (int blk = 0; blk < num_blocks; ++blk) {
      const auto b = static_cast<std::size_t>(blk);
      const double coef_mobi = this->mobility_at_points_[b].value(i);
      const double lamb = this->lambda_at_points_[b].value(i);

      // Laplacian : (grad(u), grad(psi))
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= coef_mobi * lamb * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // (w'(u), psi)
      double ww = coef_mobi * this->energy_at_points_[b].energy_derivative(i, blk);
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
 * All the blocks share the same finite element. Organised as AssembleElementVector: everything
 * that does not
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
  MFEM_VERIFY(this->same_finite_elements(el),
              "Error in AllenCahn: all the phases must be on the same finite element space.");

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
  this->coefficients_at_points(fe, Tr, *ir, elfun, true);
  this->mob_at_ip_.assign(static_cast<std::size_t>(num_blocks), 0.);
  this->dmob_at_ip_.assign(static_cast<std::size_t>(num_blocks), 0.);

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    fe.CalcShape(ip, Psi);
    this->set_int_point(Tr, ip);

    fe.CalcPhysDShape(Tr, gradPsi);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    for (int b = 0; b < num_blocks; ++b) {
      const auto ub = static_cast<std::size_t>(b);
      this->mob_at_ip_[ub] = this->mobility_at_points_[ub].value(i);
    }

    for (int blk = 0; blk < num_blocks; ++blk) {
      const auto b = static_cast<std::size_t>(blk);
      const CoefficientAtPoints& energy = this->energy_at_points_[b];
      // Diagonal block: Laplacian + w''(u)
      const double coef_mobi = this->mob_at_ip_[b];
      const double lamb = this->lambda_at_points_[b].value(i);
      const double mob_lamb = coef_mobi * lamb * weight_coef;
      AddMult_a_AAt(mob_lamb, gradPsi, *elmats(blk, blk));

      double fun_val = coef_mobi * energy.second_derivative(i, blk, blk);
      fun_val *= weight_coef;
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));

      // Off-diagonal blocks: d/du_j (M_b dW_b/du_b) = M_b d2W_b/du_b du_j (mobility and free
      // energy of the row block blk; the derivative of M_b is added below)
      for (int jblk = 0; jblk < num_blocks; ++jblk) {
        if (jblk == blk) continue;
        double fun_val_off = coef_mobi * energy.second_derivative(i, blk, jblk);
        fun_val_off *= weight_coef;
        AddMult_a_VVt(fun_val_off, Psi, *elmats(blk, jblk));
      }

      // Mobility depending on the unknowns:
      //   d/du_j (M (lambda grad u . grad psi_i + W' psi_i))
      //     += dM/du_j psi_j (lambda grad u . grad psi_i + W' psi_i)
      bool has_derivative = false;
      for (int j = 0; j < num_blocks; ++j) {
        const double d = this->mobility_at_points_[b].derivative(i, j);
        this->dmob_at_ip_[static_cast<std::size_t>(j)] = d;
        has_derivative = has_derivative || (d != 0.0);
      }
      if (has_derivative) {
        this->add_mobility_derivative(blk, lamb, energy.energy_derivative(i, blk), gradPsi, Psi,
                                      *elfun[blk], weight_coef, elmats);
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
 * @brief Adds the derivative of the mobility to the Jacobian, from known values at the current
 * point:
 *   d/du_j (M (lambda grad u . grad psi_i + W' psi_i))
 *     += dM/du_j psi_j (lambda grad u . grad psi_i + W' psi_i)
 *
 * @param blk Row block
 * @param lamb Capillary coefficient of the block
 * @param w_prime Derivative W' of the free energy of the block
 * @param gradPsi Physical gradients of the shape functions
 * @param Psi Shape functions
 * @param elfun_blk Values of the dofs of the block
 * @param weight_coef Integration weight
 * @param elmats Element matrices
 * @remark dM/du_j are read in dmob_at_ip_
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::add_mobility_derivative(
    const int blk, const double lamb, const double w_prime, const mfem::DenseMatrix& gradPsi,
    const mfem::Vector& Psi, const mfem::Vector& elfun_blk, const double weight_coef,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int num_blocks = static_cast<int>(this->nb_blk_);

  // g_i = lambda grad u . grad psi_i + W' psi_i
  this->grad_u_mob_.SetSize(gradPsi.Width());
  this->g_mob_.SetSize(gradPsi.Height());
  gradPsi.MultTranspose(elfun_blk, this->grad_u_mob_);
  gradPsi.Mult(this->grad_u_mob_, this->g_mob_);
  this->g_mob_ *= lamb;
  this->g_mob_.Add(w_prime, Psi);

  for (int j = 0; j < num_blocks; ++j) {
    const double d = this->dmob_at_ip_[static_cast<std::size_t>(j)];
    if (d == 0.0) continue;
    AddMult_a_VWt(d * weight_coef, this->g_mob_, Psi, *elmats(blk, j));
  }
}

/**
 * @brief Evaluates the coefficients of each block at all the integration points of the element
 * (blocks with the same finite element).
 *
 * Residual: mobility M, capillary coefficient lambda and W'. Jacobian: in addition dM/du_j and
 * W''. See SlothNLFormIntegrator::eval_coefficient_at_points.
 *
 * @param fe Finite element of the blocks
 * @param Tr Transformation of the element
 * @param ir Integration rule
 * @param elfun Values of the dofs of each block
 * @param jacobian Derivatives required by the Jacobian
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::coefficients_at_points(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian) {
  const int num_blocks = elfun.Size();
  const auto nb = static_cast<std::size_t>(num_blocks);
  this->mobility_at_points_.resize(nb);
  this->lambda_at_points_.resize(nb);
  this->energy_at_points_.resize(nb);

  const unsigned int mobility_outputs =
      jacobian ? (SlothNLBase::CoefValue | SlothNLBase::CoefFirstDerivative)
               : SlothNLBase::CoefValue;
  const unsigned int energy_outputs =
      jacobian ? (SlothNLBase::CoefFirstDerivative | SlothNLBase::CoefSecondDerivative)
               : SlothNLBase::CoefFirstDerivative;
  for (int blk = 0; blk < num_blocks; ++blk) {
    const auto b = static_cast<std::size_t>(blk);
    this->eval_coefficient_at_points(mobility[blk], blk, mobility_outputs, fe, Tr, ir, elfun,
                                     this->vaux_gf_, this->mobility_at_points_[b]);
    this->eval_coefficient_at_points(lambda[blk], blk, SlothNLBase::CoefValue, fe, Tr, ir, elfun,
                                     this->vaux_gf_, this->lambda_at_points_[b]);
    this->eval_coefficient_at_points(double_well_energy[blk], blk, energy_outputs, fe, Tr, ir,
                                     elfun, this->vaux_gf_, this->energy_at_points_[b]);
  }
}
