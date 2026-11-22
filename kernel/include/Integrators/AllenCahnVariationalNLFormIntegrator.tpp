/**
 * @file AllenCahnVariationalNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Allen-Cahn integrator with a variational free energy f(u, grad u)
 * @version 0.1
 * @date 2026-11-22
 *
 * Copyright CEA (C) 2026
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
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Coefficients/FieldsAtPoints.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Coefficients/VariationalEval.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AllenCahnVariationalNLFormIntegrator object.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param geometry     Geometry of the problem (Cartesian, axisymmetric).
 * @param time_step    Time step.
 * @param u_old        Fields of the blocks at the previous time step.
 * @param aux_old      Auxiliary fields at the previous time step.
 * @param params       Parameters of the integrator.
 * @param auxvars      Auxiliary variables required by the integrator.
 * @param coefficients Coefficients of each block (Mobility, variational FreeEnergy).
 */
template <class VARS>
AllenCahnVariationalNLFormIntegrator<VARS>::AllenCahnVariationalNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {
  this->integrator_name_ = "AllenCahnVariational";
}

/**
 * @brief Initialize the integrator: checks the types of the coefficients and stores them.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnVariationalNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief Retrieve the coefficients: the mobility of each block and the variational free energy
 * of the first block, function of all the blocks.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnVariationalNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
  }
  const auto energy = this->get_coefficient(0, GlossaryType::FreeEnergy, 0);
  MFEM_VERIFY(
      energy.has_value() && energy->is_variational(),
      "Error in " + this->integrator_name_ +
          ": the FreeEnergy coefficient must be a variational coefficient (\"variational\": "
          "true in the JSON file), function of the blocks and of their gradients.");
  this->free_energy = *energy;
}

/**
 * @brief Integration weight at the point ip (the integration point of Tr is set to ip)
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
double AllenCahnVariationalNLFormIntegrator<VARS>::integration_weight(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip) {
  this->set_int_point(Tr, ip);
  double weight = ip.weight * Tr.Weight();
  if (this->isAxisymmetric()) {
    weight *= AxiCylindricalCoefficient().Eval(Tr, ip);
  }
  return weight;
}

//////////////////////////////////////////////////////////////////////////////
// Coefficients at the integration points
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Evaluates the mobility of each block at all the integration points of the element (see
 * SlothNLFormIntegrator::eval_coefficient_at_points): value, and its derivatives with respect to
 * the blocks for the Jacobian.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnVariationalNLFormIntegrator<VARS>::mobility_at_points(
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
 * @brief Evaluates the derivatives of the free energy at all the integration points of the
 * element, by batch (eval_variational_batch): first derivatives with respect to the variables and
 * to the components of their gradients and, for the Jacobian of an energy which is not explicit,
 * second derivatives.
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
void AllenCahnVariationalNLFormIntegrator<VARS>::energy_at_points(
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

  // Auxiliary variables with their gradients: from their element dofs if they share the finite
  // element of the blocks, point by point otherwise
  const std::size_t naux = this->vaux_gf_.size();
  this->dofs_aux_.resize(naux);
  bool aux_dofs = true;
  for (std::size_t k = 0; k < naux && aux_dofs; ++k) {
    const mfem::FiniteElementSpace* fes = this->vaux_gf_[k].FESpace();
    aux_dofs = elem < fes->GetNE() && this->same_finite_element(*fes->GetFE(elem), fe);
    if (aux_dofs) this->vaux_gf_[k].GetElementDofValues(elem, this->dofs_aux_[k]);
  }
  if (aux_dofs) {
    fields_at_points(fe, Tr, ir, this->dofs_aux_, this->auxiliary_values_,
                     auxiliary_gradients ? &this->auxiliary_gradients_ : nullptr);
  } else {
    const auto unip = static_cast<std::size_t>(nip);
    const auto usdim = static_cast<std::size_t>(sdim);
    this->auxiliary_values_.SetSize(static_cast<int>(unip * naux));
    this->auxiliary_gradients_.SetSize(auxiliary_gradients ? static_cast<int>(unip * naux * usdim)
                                                           : 0);
    double* v = this->auxiliary_values_.HostWrite();
    double* g = auxiliary_gradients ? this->auxiliary_gradients_.HostWrite() : nullptr;
    mfem::Vector grad(sdim);
    for (int q = 0; q < nip; ++q) {
      const mfem::IntegrationPoint& ip = ir.IntPoint(q);
      this->set_int_point(Tr, ip);
      const auto uq = static_cast<std::size_t>(q);
      for (std::size_t k = 0; k < naux; ++k) {
        v[uq * naux + k] = this->vaux_gf_[k].GetValue(Tr, ip);
        if (g == nullptr) continue;
        this->vaux_gf_[k].GetGradient(Tr, grad);
        for (std::size_t c = 0; c < usdim; ++c) {
          g[(uq * naux + k) * usdim + c] = grad(static_cast<int>(c));
        }
      }
    }
  }
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
  this->with_second_derivatives_ = jacobian && !explicit_energy;
  this->df_dvar_.SetSize(nip * nb);
  this->df_dgrad_.SetSize(nip * nd);
  VariationalBatch outputs;
  outputs.df_dvar = &this->df_dvar_;
  outputs.df_dgrad = &this->df_dgrad_;
  if (this->with_second_derivatives_) {
    this->d2f_dvar2_.SetSize(nip * nb * nb);
    this->d2f_dvar_dgrad_.SetSize(nip * nb * nd);
    this->d2f_dgrad2_.SetSize(nip * nd * nd);
    outputs.d2f_dvar2 = &this->d2f_dvar2_;
    outputs.d2f_dvar_dgrad = &this->d2f_dvar_dgrad_;
    outputs.d2f_dgrad2 = &this->d2f_dgrad2_;
  }
  MFEM_VERIFY(
      f.eval_variational_batch(nip, this->values_, this->explicit_values_, this->auxiliary_values_,
                               this->gradients_, this->explicit_gradients_,
                               this->auxiliary_gradients_, this->coordinates_, outputs),
      "Error in " + this->integrator_name_ + ": evaluation of the free energy.");
}

//////////////////////////////////////////////////////////////////////////////
// Residual
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Assemble the element residual of each block b:
 *   R_b(psi) = int M_b (df/du_b psi + df/d(grad u_b) . grad psi).
 *
 * All the blocks share the same finite element. The mobility and the derivatives of the free
 * energy are evaluated once per element at all the integration points.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnVariationalNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int num_blocks = el.Size();
  MFEM_VERIFY(this->same_finite_elements(el), "Error in " + this->integrator_name_ +
                                                  ": all the blocks must be on the same finite "
                                                  "element space.");

  const mfem::FiniteElement& fe = *el[0];
  const int nd = fe.GetDof();
  const int sdim = Tr.GetSpaceDim();
  const int ng = num_blocks * sdim;  // components of the gradients of the blocks
  this->Psi_.SetSize(nd);
  this->gradPsi_.SetSize(nd, sdim);
  this->flux_.SetSize(sdim);
  for (int blk = 0; blk < num_blocks; ++blk) {
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
  }

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
  this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, *ir);
  this->mobility_at_points(fe, Tr, *ir, elfun, false);
  this->energy_at_points(fe, Tr, *ir, elfun, false);
  const double* df_dvar = this->df_dvar_.HostRead();
  const double* df_dgrad = this->df_dgrad_.HostRead();

  for (int q = 0; q < ir->GetNPoints(); q++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(q);
    fe.CalcShape(ip, this->Psi_);
    const double weight = this->integration_weight(Tr, ip);
    fe.CalcPhysDShape(Tr, this->gradPsi_);

    for (int blk = 0; blk < num_blocks; ++blk) {
      const double a = this->mobility_at_points_[static_cast<std::size_t>(blk)].value(q) * weight;
      // (df/d(grad u_b), grad psi)
      for (int k = 0; k < sdim; ++k) this->flux_(k) = df_dgrad[q * ng + blk * sdim + k];
      this->gradPsi_.AddMult_a(a, this->flux_, *elvect[blk]);
      // (df/du_b, psi)
      elvect[blk]->Add(a * df_dvar[q * num_blocks + blk], this->Psi_);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Jacobian
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Assemble the element Jacobian: block (b, j), derivative of R_b with respect to the dofs
 * of the block j (shape functions psi, psi_j):
 *   int M_b (d2f/du_b du_j psi psi_j + (d2f/du_b d(grad u_j) . grad psi_j) psi
 *            + (d2f/d(grad u_b) du_j . grad psi) psi_j
 *            + grad psi . d2f/d(grad u_b) d(grad u_j) grad psi_j)
 *   + int dM_b/du_j psi_j (df/du_b psi + df/d(grad u_b) . grad psi).
 * The terms of the free energy are not added for an explicit energy.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void AllenCahnVariationalNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int num_blocks = el.Size();
  MFEM_VERIFY(this->same_finite_elements(el), "Error in " + this->integrator_name_ +
                                                  ": all the blocks must be on the same finite "
                                                  "element space.");

  const mfem::FiniteElement& fe = *el[0];
  const int nd = fe.GetDof();
  const int sdim = Tr.GetSpaceDim();
  const int nb = num_blocks;
  const int ng = nb * sdim;  // components of the gradients of the blocks
  this->Psi_.SetSize(nd);
  this->gradPsi_.SetSize(nd, sdim);
  this->flux_.SetSize(sdim);
  this->row_.SetSize(nd);
  this->column_.SetSize(nd);
  this->block_hessian_.SetSize(sdim, sdim);
  this->gradPsi_hessian_.SetSize(nd, sdim);
  this->residual_integrand_.resize(static_cast<std::size_t>(nb));
  for (int blk = 0; blk < nb; ++blk) {
    for (int jblk = 0; jblk < nb; ++jblk) {
      elmats(blk, jblk)->SetSize(nd);
      *elmats(blk, jblk) = 0.0;
    }
  }

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
  this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, *ir);
  this->mobility_at_points(fe, Tr, *ir, elfun, true);
  this->energy_at_points(fe, Tr, *ir, elfun, true);
  const double* df_dvar = this->df_dvar_.HostRead();
  const double* df_dgrad = this->df_dgrad_.HostRead();
  const double* d2f_dvar2 = this->with_second_derivatives_ ? this->d2f_dvar2_.HostRead() : nullptr;
  const double* d2f_dvar_dgrad =
      this->with_second_derivatives_ ? this->d2f_dvar_dgrad_.HostRead() : nullptr;
  const double* d2f_dgrad2 =
      this->with_second_derivatives_ ? this->d2f_dgrad2_.HostRead() : nullptr;

  for (int q = 0; q < ir->GetNPoints(); q++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(q);
    fe.CalcShape(ip, this->Psi_);
    const double weight = this->integration_weight(Tr, ip);
    fe.CalcPhysDShape(Tr, this->gradPsi_);

    for (int b = 0; b < nb; ++b) {
      const CoefficientAtPoints& mob = this->mobility_at_points_[static_cast<std::size_t>(b)];
      const double a = mob.value(q) * weight;

      // Second derivatives of the free energy
      for (int j = 0; j < nb && d2f_dvar2 != nullptr; ++j) {
        mfem::DenseMatrix& mat = *elmats(b, j);
        // d2f/du_b du_j psi psi_j
        mfem::AddMult_a_VVt(a * d2f_dvar2[(q * nb + b) * nb + j], this->Psi_, mat);
        // (d2f/du_b d(grad u_j) . grad psi_j) psi
        for (int k = 0; k < sdim; ++k) {
          this->flux_(k) = d2f_dvar_dgrad[(q * nb + b) * ng + j * sdim + k];
        }
        this->gradPsi_.Mult(this->flux_, this->column_);
        mfem::AddMult_a_VWt(a, this->Psi_, this->column_, mat);
        // (d2f/d(grad u_b) du_j . grad psi) psi_j
        for (int k = 0; k < sdim; ++k) {
          this->flux_(k) = d2f_dvar_dgrad[(q * nb + j) * ng + b * sdim + k];
        }
        this->gradPsi_.Mult(this->flux_, this->row_);
        mfem::AddMult_a_VWt(a, this->row_, this->Psi_, mat);
        // grad psi . d2f/d(grad u_b) d(grad u_j) grad psi_j
        for (int k = 0; k < sdim; ++k) {
          for (int l = 0; l < sdim; ++l) {
            this->block_hessian_(k, l) = d2f_dgrad2[(q * ng + b * sdim + k) * ng + j * sdim + l];
          }
        }
        mfem::Mult(this->gradPsi_, this->block_hessian_, this->gradPsi_hessian_);
        mfem::AddMult_a_ABt(a, this->gradPsi_hessian_, this->gradPsi_, mat);
      }

      // Mobility depending on the blocks: dM_b/du_j psi_j (df/du_b psi + df/d(grad u_b) . grad psi)
      bool has_derivative = false;
      for (int j = 0; j < nb && !has_derivative; ++j) has_derivative = mob.derivative(q, j) != 0.0;
      if (!has_derivative) continue;
      mfem::Vector& g = this->residual_integrand_[static_cast<std::size_t>(b)];
      g.SetSize(nd);
      for (int k = 0; k < sdim; ++k) this->flux_(k) = df_dgrad[q * ng + b * sdim + k];
      this->gradPsi_.Mult(this->flux_, g);
      g.Add(df_dvar[q * nb + b], this->Psi_);
      for (int j = 0; j < nb; ++j) {
        const double d = mob.derivative(q, j);
        if (d == 0.0) continue;
        mfem::AddMult_a_VWt(d * weight, g, this->Psi_, *elmats(b, j));
      }
    }
  }
}
