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
#include "Coefficients/FieldsAtPoints.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Coefficients/VariationalEval.hpp"
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
 * 1. Checks that all coefficients contain the expected types: Mobility, FreeEnergy and Capillary,
 *    or Mobility and FreeEnergy if the FreeEnergy of the first block is variational (no Capillary
 *    then: the gradient terms are part of the free energy).
 * 2. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::init() {
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
  // Variational free energy: kappa_b and df/du_b (isotropic), all the derivatives otherwise
  const double* kappa = nullptr;
  const double* df_dvar = nullptr;
  if (this->variational_) {
    this->mobility_at_points(fe, Tr, *ir, elfun, false);
    this->variational_energy_at_points(fe, Tr, *ir, elfun, false);
    if (!this->isotropic_) {
      this->variational_element_vector(fe, Tr, *ir, elvect);
      return;
    }
    kappa = this->kappa_.HostRead();
    df_dvar = this->df_dvar_.HostRead();
  } else {
    this->coefficients_at_points(fe, Tr, *ir, elfun, false);
  }

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
      const double lamb =
          kappa != nullptr ? kappa[i * num_blocks + blk] : this->lambda_at_points_[b].value(i);
      const double w_prime = df_dvar != nullptr
                                 ? df_dvar[i * num_blocks + blk]
                                 : this->energy_at_points_[b].energy_derivative(i, blk);

      // Laplacian : (grad(u), grad(psi))
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= coef_mobi * lamb * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // (w'(u), psi)
      double ww = coef_mobi * w_prime;
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
  // Variational free energy: kappa_b, df/du_b and d2f/du_b du_j (isotropic), all the derivatives
  // otherwise
  const double* kappa = nullptr;
  const double* df_dvar = nullptr;
  const double* d2f_dvar2 = nullptr;
  if (this->variational_) {
    this->mobility_at_points(fe, Tr, *ir, elfun, true);
    this->variational_energy_at_points(fe, Tr, *ir, elfun, true);
    if (!this->isotropic_) {
      this->variational_element_grad(fe, Tr, *ir, elmats);
      return;
    }
    kappa = this->kappa_.HostRead();
    df_dvar = this->df_dvar_.HostRead();
    if (this->with_second_derivatives_) d2f_dvar2 = this->d2f_dvar2_.HostRead();
  } else {
    this->coefficients_at_points(fe, Tr, *ir, elfun, true);
  }
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
      // W'(u) and W''(u) of the block: classic free energy of the block, or derivatives of the
      // variational free energy (no second derivative for an explicit energy)
      const bool with_w_second = !this->variational_ || d2f_dvar2 != nullptr;
      auto w_second = [&](const int jblk) {
        return d2f_dvar2 != nullptr ? d2f_dvar2[(i * num_blocks + blk) * num_blocks + jblk]
                                    : this->energy_at_points_[b].second_derivative(i, blk, jblk);
      };
      const double w_prime = df_dvar != nullptr
                                 ? df_dvar[i * num_blocks + blk]
                                 : this->energy_at_points_[b].energy_derivative(i, blk);
      // Diagonal block: Laplacian + w''(u)
      const double coef_mobi = this->mob_at_ip_[b];
      const double lamb =
          kappa != nullptr ? kappa[i * num_blocks + blk] : this->lambda_at_points_[b].value(i);
      const double mob_lamb = coef_mobi * lamb * weight_coef;
      AddMult_a_AAt(mob_lamb, gradPsi, *elmats(blk, blk));

      if (with_w_second) {
        double fun_val = coef_mobi * w_second(blk);
        fun_val *= weight_coef;
        AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));

        // Off-diagonal blocks: d/du_j (M_b dW_b/du_b) = M_b d2W_b/du_b du_j (mobility and free
        // energy of the row block blk; the derivative of M_b is added below)
        for (int jblk = 0; jblk < num_blocks; ++jblk) {
          if (jblk == blk) continue;
          double fun_val_off = coef_mobi * w_second(jblk);
          fun_val_off *= weight_coef;
          AddMult_a_VVt(fun_val_off, Psi, *elmats(blk, jblk));
        }
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
        this->add_mobility_derivative(blk, lamb, w_prime, gradPsi, Psi, *elfun[blk], weight_coef,
                                      elmats);
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
 * Variational free energy: mobility of each block and FreeEnergy of the first block (free_energy).
 * Otherwise, this method collects the coefficients of type `Capillary`, `Mobility`, and
 * `FreeEnergy` from each coefficients and adds them to the corresponding internal storage:
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
  if (this->variational_) {
    // Mobility of each block and variational free energy of the first block
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

//////////////////////////////////////////////////////////////////////////////
// Variational free energy
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Evaluates the mobility of each block at all the integration points of the element (see
 * SlothNLFormIntegrator::eval_coefficient_at_points): value, and its derivatives with respect to
 * the blocks for the Jacobian.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::mobility_at_points(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian) {
  const int num_blocks = elfun.Size();
  this->mobility_at_points_.resize(static_cast<std::size_t>(num_blocks));
  const unsigned int outputs = jacobian
                                   ? (SlothNLBase::CoefValue | SlothNLBase::CoefFirstDerivative)
                                   : SlothNLBase::CoefValue;
  for (int blk = 0; blk < num_blocks; ++blk) {
    this->eval_coefficient_at_points(mobility[blk], blk, outputs, fe, Tr, ir, elfun, this->vaux_gf_,
                                     this->mobility_at_points_[static_cast<std::size_t>(blk)]);
  }
}

/**
 * @brief Evaluates the derivatives of the variational free energy at all the integration points of
 * the element, by batch (eval_variational_batch): first derivatives with respect to the variables
 * and to the components of their gradients and, for the Jacobian of an energy which is not
 * explicit, second derivatives. Isotropic energy: kappa_b in place of the derivatives with respect
 * to the gradients (kappa_b grad u_b and kappa_b I).
 *
 * Values given to the energy: values and gradients of u (implicit energy) or of u_old (explicit
 * energy); for a semi-implicit energy, variables u and explicit variables u_old. Auxiliary
 * variables: values and gradients of the auxiliary fields.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 * @param fe Finite element of the blocks
 * @param Tr Transformation of the element
 * @param ir Integration rule
 * @param elfun Values of the dofs of each block
 * @param jacobian Second derivatives required by the Jacobian
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::variational_energy_at_points(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian) {
  const Coefficient& f = *this->free_energy;
  const int nb = elfun.Size();
  const auto unb = static_cast<std::size_t>(nb);
  const int nip = ir.GetNPoints();
  const int sdim = Tr.GetSpaceDim();
  const int elem = Tr.ElementNo;
  MFEM_VERIFY(f.variational_dimension() == sdim,
              "Error in " + this->integrator_name_ + ": dimension of the free energy (" +
                  std::to_string(f.variational_dimension()) +
                  ") different from the dimension of the mesh (" + std::to_string(sdim) + ").");
  MFEM_VERIFY(f.variational_nb_variables() == nb,
              "Error in " + this->integrator_name_ + ": free energy of " +
                  std::to_string(f.variational_nb_variables()) + " variables for " +
                  std::to_string(nb) + " blocks: the variables of the free energy are the blocks.");
  const bool explicit_energy = f.is_explicit();
  const bool semi_implicit = f.is_semi_implicit();
  const bool explicit_gradients = semi_implicit && f.uses_explicit_gradients();
  const bool auxiliary_gradients = f.uses_auxiliary_gradients();

  // Element dofs of u and, if needed, of u_old
  this->dofs_u_.resize(unb);
  for (std::size_t b = 0; b < unb; ++b) {
    this->dofs_u_[b] = *elfun[static_cast<int>(b)];
  }
  if (explicit_energy || semi_implicit) {
    this->dofs_u_old_.resize(unb);
    for (std::size_t b = 0; b < unb; ++b) {
      const mfem::FiniteElementSpace* fes = this->u_old_[b].FESpace();
      MFEM_VERIFY(elem < fes->GetNE() && this->same_finite_element(*fes->GetFE(elem), fe),
                  "Error in " + this->integrator_name_ +
                      ": the previous values of the blocks must be on the finite element space of "
                      "the blocks.");
      this->u_old_[b].GetElementDofValues(elem, this->dofs_u_old_[b]);
    }
  }

  // Variables (u or u_old) and explicit variables (u_old) with their gradients
  fields_at_points(fe, Tr, ir, explicit_energy ? this->dofs_u_old_ : this->dofs_u_, this->values_,
                   &this->gradients_);
  if (semi_implicit) {
    fields_at_points(fe, Tr, ir, this->dofs_u_old_, this->explicit_values_,
                     explicit_gradients ? &this->explicit_gradients_ : nullptr);
  } else {
    this->explicit_values_.SetSize(0);
  }
  if (!explicit_gradients) this->explicit_gradients_.SetSize(0);

  // Auxiliary variables with their gradients
  this->auxiliary_fields_at_points(fe, Tr, ir, this->vaux_gf_, this->auxiliary_values_,
                                   auxiliary_gradients ? &this->auxiliary_gradients_ : nullptr);
  if (!auxiliary_gradients) this->auxiliary_gradients_.SetSize(0);

  // Physical coordinates of the points
  if (f.uses_coordinates()) {
    Tr.Transform(ir, this->points_);  // sdim x nip, column q = point q
    this->coordinates_.SetSize(this->points_.Height() * nip);
    std::copy(this->points_.GetData(), this->points_.GetData() + this->coordinates_.Size(),
              this->coordinates_.HostWrite());
  } else {
    this->coordinates_.SetSize(0);
  }

  // Derivatives of the free energy (second ones: no contribution of an explicit energy)
  const int nd = nb * sdim;
  const bool isotropic = this->isotropic_;
  this->with_second_derivatives_ = jacobian && !explicit_energy;
  this->df_dvar_.SetSize(nip * nb);
  VariationalBatch outputs;
  outputs.df_dvar = &this->df_dvar_;
  if (isotropic) {
    this->kappa_.SetSize(nip * nb);
    outputs.d2f_dgrad2_isotropic = &this->kappa_;
  } else {
    this->df_dgrad_.SetSize(nip * nd);
    outputs.df_dgrad = &this->df_dgrad_;
  }
  if (this->with_second_derivatives_) {
    this->d2f_dvar2_.SetSize(nip * nb * nb);
    outputs.d2f_dvar2 = &this->d2f_dvar2_;
    if (!isotropic) {
      this->d2f_dvar_dgrad_.SetSize(nip * nb * nd);
      this->d2f_dgrad2_.SetSize(nip * nd * nd);
      outputs.d2f_dvar_dgrad = &this->d2f_dvar_dgrad_;
      outputs.d2f_dgrad2 = &this->d2f_dgrad2_;
    }
  }
  MFEM_VERIFY(
      f.eval_variational_batch(nip, this->values_, this->explicit_values_, this->auxiliary_values_,
                               this->gradients_, this->explicit_gradients_,
                               this->auxiliary_gradients_, this->coordinates_, outputs),
      "Error in " + this->integrator_name_ + ": evaluation of the free energy.");
}

/**
 * @brief Integration weight at the point ip (the integration point of Tr is set to ip)
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
double AllenCahnNLFormIntegrator<VARS>::integration_weight(mfem::ElementTransformation& Tr,
                                                           const mfem::IntegrationPoint& ip) {
  this->set_int_point(Tr, ip);
  double weight = ip.weight * Tr.Weight();
  if (this->isAxisymmetric()) {
    weight *= AxiCylindricalCoefficient().Eval(Tr, ip);
  }
  return weight;
}

/**
 * @brief Element residual with a variational free energy which is not isotropic, from the
 * mobility and the derivatives of f at the points (mobility_at_points,
 * variational_energy_at_points):
 *   R_b(psi) = int M_b (df/du_b psi + df/d(grad u_b) . grad psi).
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::variational_element_vector(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<mfem::Vector*>& elvect) {
  const int num_blocks = elvect.Size();
  const int nd = fe.GetDof();
  const int sdim = Tr.GetSpaceDim();
  const int ng = num_blocks * sdim;  // components of the gradients of the blocks
  this->Psi.SetSize(nd);
  this->gradPsi.SetSize(nd, sdim);
  this->flux_.SetSize(sdim);
  const double* df_dvar = this->df_dvar_.HostRead();
  const double* df_dgrad = this->df_dgrad_.HostRead();

  for (int q = 0; q < ir.GetNPoints(); q++) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    fe.CalcShape(ip, this->Psi);
    const double weight = this->integration_weight(Tr, ip);
    fe.CalcPhysDShape(Tr, this->gradPsi);

    for (int blk = 0; blk < num_blocks; ++blk) {
      const double a = this->mobility_at_points_[static_cast<std::size_t>(blk)].value(q) * weight;
      // (df/d(grad u_b), grad psi)
      for (int k = 0; k < sdim; ++k) this->flux_(k) = df_dgrad[q * ng + blk * sdim + k];
      this->gradPsi.AddMult_a(a, this->flux_, *elvect[blk]);
      // (df/du_b, psi)
      elvect[blk]->Add(a * df_dvar[q * num_blocks + blk], this->Psi);
    }
  }
}

/**
 * @brief Element Jacobian with a variational free energy which is not isotropic: block (b, j),
 * derivative of R_b with respect to the dofs of the block j (shape functions psi, psi_j):
 *   int M_b (d2f/du_b du_j psi psi_j + (d2f/du_b d(grad u_j) . grad psi_j) psi
 *            + (d2f/d(grad u_b) du_j . grad psi) psi_j
 *            + grad psi . d2f/d(grad u_b) d(grad u_j) grad psi_j)
 *   + int dM_b/du_j psi_j (df/du_b psi + df/d(grad u_b) . grad psi).
 * The terms of the free energy are not added for an explicit energy.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::variational_element_grad(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int nb = elmats.NumRows();
  const int nd = fe.GetDof();
  const int sdim = Tr.GetSpaceDim();
  const int ng = nb * sdim;  // components of the gradients of the blocks
  this->Psi.SetSize(nd);
  this->gradPsi.SetSize(nd, sdim);
  this->flux_.SetSize(sdim);
  this->row_.SetSize(nd);
  this->column_.SetSize(nd);
  this->block_hessian_.SetSize(sdim, sdim);
  this->gradPsi_hessian_.SetSize(nd, sdim);
  this->residual_integrand_.resize(static_cast<std::size_t>(nb));
  const double* df_dvar = this->df_dvar_.HostRead();
  const double* df_dgrad = this->df_dgrad_.HostRead();
  const double* d2f_dvar2 = this->with_second_derivatives_ ? this->d2f_dvar2_.HostRead() : nullptr;
  const double* d2f_dvar_dgrad =
      this->with_second_derivatives_ ? this->d2f_dvar_dgrad_.HostRead() : nullptr;
  const double* d2f_dgrad2 =
      this->with_second_derivatives_ ? this->d2f_dgrad2_.HostRead() : nullptr;

  for (int q = 0; q < ir.GetNPoints(); q++) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    fe.CalcShape(ip, this->Psi);
    const double weight = this->integration_weight(Tr, ip);
    fe.CalcPhysDShape(Tr, this->gradPsi);

    for (int b = 0; b < nb; ++b) {
      const CoefficientAtPoints& mob = this->mobility_at_points_[static_cast<std::size_t>(b)];
      const double a = mob.value(q) * weight;

      // Second derivatives of the free energy
      for (int j = 0; j < nb && d2f_dvar2 != nullptr; ++j) {
        mfem::DenseMatrix& mat = *elmats(b, j);
        // d2f/du_b du_j psi psi_j
        mfem::AddMult_a_VVt(a * d2f_dvar2[(q * nb + b) * nb + j], this->Psi, mat);
        // (d2f/du_b d(grad u_j) . grad psi_j) psi
        for (int k = 0; k < sdim; ++k) {
          this->flux_(k) = d2f_dvar_dgrad[(q * nb + b) * ng + j * sdim + k];
        }
        this->gradPsi.Mult(this->flux_, this->column_);
        mfem::AddMult_a_VWt(a, this->Psi, this->column_, mat);
        // (d2f/d(grad u_b) du_j . grad psi) psi_j
        for (int k = 0; k < sdim; ++k) {
          this->flux_(k) = d2f_dvar_dgrad[(q * nb + j) * ng + b * sdim + k];
        }
        this->gradPsi.Mult(this->flux_, this->row_);
        mfem::AddMult_a_VWt(a, this->row_, this->Psi, mat);
        // grad psi . d2f/d(grad u_b) d(grad u_j) grad psi_j
        for (int k = 0; k < sdim; ++k) {
          for (int l = 0; l < sdim; ++l) {
            this->block_hessian_(k, l) = d2f_dgrad2[(q * ng + b * sdim + k) * ng + j * sdim + l];
          }
        }
        mfem::Mult(this->gradPsi, this->block_hessian_, this->gradPsi_hessian_);
        mfem::AddMult_a_ABt(a, this->gradPsi_hessian_, this->gradPsi, mat);
      }

      // Mobility depending on the blocks: dM_b/du_j psi_j (df/du_b psi + df/d(grad u_b) . grad psi)
      bool has_derivative = false;
      for (int j = 0; j < nb && !has_derivative; ++j) has_derivative = mob.derivative(q, j) != 0.0;
      if (!has_derivative) continue;
      mfem::Vector& g = this->residual_integrand_[static_cast<std::size_t>(b)];
      g.SetSize(nd);
      for (int k = 0; k < sdim; ++k) this->flux_(k) = df_dgrad[q * ng + b * sdim + k];
      this->gradPsi.Mult(this->flux_, g);
      g.Add(df_dvar[q * nb + b], this->Psi);
      for (int j = 0; j < nb; ++j) {
        const double d = mob.derivative(q, j);
        if (d == 0.0) continue;
        mfem::AddMult_a_VWt(d * weight, g, this->Psi, *elmats(b, j));
      }
    }
  }
}
