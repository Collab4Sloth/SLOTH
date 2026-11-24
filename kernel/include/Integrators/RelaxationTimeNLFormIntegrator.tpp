/**
 * @file RelaxationTimeNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative with a relaxation time tau(u, grad u)
 * @version 0.1
 * @date 2026-11-23
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
#include <string>
#include <vector>

#include "Coefficients/FieldsAtPoints.hpp"
#include "Coefficients/VariationalEval.hpp"
#include "Integrators/RelaxationTimeNLFormIntegrator.hpp"

/**
 * @brief Construct a new RelaxationTimeNLFormIntegrator object.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 */
template <class VARS>
RelaxationTimeNLFormIntegrator<VARS>::RelaxationTimeNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction> u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : TimeNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                 coefficients) {
  this->integrator_name_ = "RelaxationTimeDerivative";
}

/**
 * @brief Coefficients of the integrator: A = tau (Glossary::RelaxationTime, one by default) and
 * B = 1 for each block.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void RelaxationTimeNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    auto tau = this->get_coefficient(i, GlossaryType::RelaxationTime, 0);
    this->coefficient_A.add(tau.has_value() ? *tau : Coefficient(Glossary::Default, 1.0));
    this->coefficient_B.add(Coefficient(Glossary::Default, 1.0));
  }
}

/**
 * @brief Relaxation time of block blk at all the integration points of the element: variational
 * coefficient (see variational_tau_at_points) or coefficient of the base class.
 */
template <class VARS>
void RelaxationTimeNLFormIntegrator<VARS>::eval_a_at_points(
    const int blk, const mfem::Array<const mfem::FiniteElement*>& el,
    mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian) {
  this->tau_.active = this->coefficient_A[blk].is_variational();
  if (this->tau_.active) {
    this->variational_tau_at_points(this->coefficient_A[blk], el, Tr, ir, elfun, jacobian);
  } else {
    TimeNLFormIntegrator<VARS>::eval_a_at_points(blk, el, Tr, ir, elfun, jacobian);
  }
}

/**
 * @brief Value of tau at the integration point i
 */
template <class VARS>
double RelaxationTimeNLFormIntegrator<VARS>::a_value(const int i) const {
  return this->tau_.active ? this->tau_.F(i) : TimeNLFormIntegrator<VARS>::a_value(i);
}

/**
 * @brief Derivative of tau with respect to block jblk at the integration point i (variational
 * tau: see add_a_jacobian)
 */
template <class VARS>
double RelaxationTimeNLFormIntegrator<VARS>::a_derivative(const int i, const int jblk) const {
  return this->tau_.active ? 0.0 : TimeNLFormIntegrator<VARS>::a_derivative(i, jblk);
}

/**
 * @brief Jacobian terms of a variational tau(u, grad u) at the integration point i (d/dk_j =
 * dt d/du_j, factor = dt B k w):
 *   factor (dtau/du_j psi_j + dtau/d(grad u_j) . grad psi_j) psi
 */
template <class VARS>
void RelaxationTimeNLFormIntegrator<VARS>::add_a_jacobian(
    const int blk, const int i, const mfem::Array<const mfem::FiniteElement*>& el,
    mfem::ElementTransformation& Tr, const double factor, const mfem::Vector& psi,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const VariationalTau& v = this->tau_;
  if (!v.active || !v.with_derivatives || factor == 0.0) return;
  const double* dtau_du = v.dtau_du.HostRead();
  const double* dtau_dgrad = v.dtau_dgrad.HostRead();
  const mfem::IntegrationPoint& ip = Tr.GetIntPoint();
  for (int jblk = 0; jblk < v.nb; ++jblk) {
    const int nd_j = el[jblk]->GetDof();
    shape_j_.SetSize(nd_j);
    el[jblk]->CalcShape(ip, shape_j_);
    gradPsi_j_.SetSize(nd_j, v.dim);
    el[jblk]->CalcPhysDShape(Tr, gradPsi_j_);
    row_j_.SetSize(nd_j);
    row_j_.Set(dtau_du[i * v.nb + jblk], shape_j_);
    const mfem::Vector dgrad(const_cast<double*>(dtau_dgrad) + (i * v.nb + jblk) * v.dim, v.dim);
    gradPsi_j_.AddMult(dgrad, row_j_);
    AddMult_a_VWt(factor, psi, row_j_, *elmats(blk, jblk));
  }
}

/**
 * @brief Variational relaxation time tau(u, grad u) at all the integration points of the element,
 * with the values of the variables at the end of the step u = u_old + dt k (the unknown of the
 * integrator is k = du/dt): value and, for the Jacobian of an implicit or semi-implicit tau, first
 * derivatives with respect to the variables and to the components of their gradients.
 * Explicit tau: values and gradients of u_old. Semi-implicit tau: explicit variables from u_old.
 * The variables of tau are the first blocks, on the same finite element.
 */
template <class VARS>
void RelaxationTimeNLFormIntegrator<VARS>::variational_tau_at_points(
    const Coefficient& tau, const mfem::Array<const mfem::FiniteElement*>& el,
    mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian) {
  VariationalTau& v = this->tau_;
  const int nb = tau.variational_nb_variables();
  const int sdim = Tr.GetSpaceDim();
  const int nip = ir.GetNPoints();
  const int elem = Tr.ElementNo;
  MFEM_VERIFY(tau.variational_dimension() == sdim,
              "Error in " + this->integrator_name_ +
                  ": dimension of the variational relaxation time (" +
                  std::to_string(tau.variational_dimension()) +
                  ") different from the dimension of the mesh (" + std::to_string(sdim) + ").");
  MFEM_VERIFY(nb >= 1 && nb <= el.Size(),
              "Error in " + this->integrator_name_ + ": variational relaxation time of " +
                  std::to_string(nb) + " variables for " + std::to_string(el.Size()) +
                  " blocks: the variables of the relaxation time are the first blocks.");
  const mfem::FiniteElement& fe = *el[0];
  for (int b = 1; b < nb; ++b) {
    MFEM_VERIFY(this->same_finite_element(*el[b], fe),
                "Error in " + this->integrator_name_ +
                    ": the variables of the variational relaxation time must share the same "
                    "finite element.");
  }
  const bool explicit_tau = tau.is_explicit();
  const bool semi_implicit = tau.is_semi_implicit();
  const bool explicit_gradients = semi_implicit && tau.uses_explicit_gradients();
  const bool auxiliary_gradients = tau.uses_auxiliary_gradients();
  v.nb = nb;
  v.dim = sdim;

  // Element dofs of u_old and of u = u_old + dt k
  const auto unb = static_cast<std::size_t>(nb);
  v.dofs_old.resize(unb);
  v.dofs.resize(unb);
  for (std::size_t b = 0; b < unb; ++b) {
    const mfem::FiniteElementSpace* fes = this->u_old_[b].FESpace();
    MFEM_VERIFY(elem < fes->GetNE() && this->same_finite_element(*fes->GetFE(elem), fe),
                "Error in " + this->integrator_name_ +
                    ": the previous values of the blocks must be on the finite element space of "
                    "the blocks.");
    this->u_old_[b].GetElementDofValues(elem, v.dofs_old[b]);
    v.dofs[b] = v.dofs_old[b];
    v.dofs[b].Add(this->time_step_, *elfun[static_cast<int>(b)]);
  }

  // Variables (u or u_old) and explicit variables (u_old) with their gradients
  fields_at_points(fe, Tr, ir, explicit_tau ? v.dofs_old : v.dofs, v.values, &v.gradients);
  if (semi_implicit) {
    fields_at_points(fe, Tr, ir, v.dofs_old, v.explicit_values,
                     explicit_gradients ? &v.explicit_gradients : nullptr);
  } else {
    v.explicit_values.SetSize(0);
  }
  if (!explicit_gradients) v.explicit_gradients.SetSize(0);

  // Auxiliary variables with their gradients
  this->auxiliary_fields_at_points(fe, Tr, ir, this->vaux_gf_, v.auxiliary_values,
                                   auxiliary_gradients ? &v.auxiliary_gradients : nullptr);
  if (!auxiliary_gradients) v.auxiliary_gradients.SetSize(0);

  // Physical coordinates of the points
  if (tau.uses_coordinates()) {
    Tr.Transform(ir, v.points);  // sdim x nip, column q = point q
    v.coordinates.SetSize(v.points.Height() * nip);
    std::copy(v.points.GetData(), v.points.GetData() + v.coordinates.Size(),
              v.coordinates.HostWrite());
  } else {
    v.coordinates.SetSize(0);
  }

  // Value and first derivatives (no derivative for an explicit tau)
  v.with_derivatives = jacobian && !explicit_tau;
  VariationalBatch outputs;
  v.F.SetSize(nip);
  outputs.F = &v.F;
  if (v.with_derivatives) {
    v.dtau_du.SetSize(nip * nb);
    v.dtau_dgrad.SetSize(nip * nb * sdim);
    outputs.df_dvar = &v.dtau_du;
    outputs.df_dgrad = &v.dtau_dgrad;
  }
  MFEM_VERIFY(
      tau.eval_variational_batch(nip, v.values, v.explicit_values, v.auxiliary_values, v.gradients,
                                 v.explicit_gradients, v.auxiliary_gradients, v.coordinates,
                                 outputs),
      "Error in " + this->integrator_name_ + ": evaluation of the variational relaxation time.");
}
