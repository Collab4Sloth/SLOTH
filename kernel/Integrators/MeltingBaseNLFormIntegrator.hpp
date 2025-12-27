/**
 * @file MeltingBaseNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for melting integrators
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
#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Base class for melting integrators
 *
 * @tparam VARS
 */
template <class VARS>
class MeltingBaseNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  std::list<GlossaryType> expected_list_{GlossaryType::Mobility, GlossaryType::PhaseFieldPotential};
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

 protected:
  Coefficients mobility;
  Coefficients interpolation_potential;

  virtual double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                        const mfem::IntegrationPoint& ir, unsigned int blk,
                                        const double u, const double un) = 0;

  virtual double compute_coefficient(Coefficient coef, const std::vector<double>& values);

  virtual double compute_gradient_coefficient(Coefficient coef, const int blk,
                                              const std::vector<double>& values);

  virtual double compute_hessian_coefficient(Coefficient coef, const int iblk, const int jblk,
                                             const std::vector<double>& values);

  void get_coefficients() override;

 public:
  MeltingBaseNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                              const Parameters& params, std::vector<VARS*> auxvars,
                              const std::vector<Coefficients>& coefficients);

   void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvec) override;

   void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
  void init() override;
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Melting Base N L Form Integrator< V A R S>:: Melting Base N L Form
 * Integrator object
 *
 * @tparam VARS
 * @param u_old
 * @param params
 * @param auxvars
 * @param coefficients
 */
template <class VARS>
MeltingBaseNLFormIntegrator<VARS>::MeltingBaseNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {}

/**
 * @brief Initialization
 *
 * @tparam VARS
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();
}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam VARS
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Catch_Time_Section("MeltingBaseNLFormIntegrator:AssembleElementVector");
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    // elvect.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      const auto& u = *elfun[blk] * Psi;
      const auto& un = this->u_old_[blk].GetValue(Tr, ip);

      const double coef_mobi = mobility[blk].compute() * ip.weight * Tr.Weight();
      const double alpha = this->get_phase_change_at_ip(Tr, ip, blk, u, un);
      const double ww =
          coef_mobi * alpha *
          this->compute_gradient_coefficient(interpolation_potential[blk], blk, {u, un});
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam VARS
 * @param el
 * @param Tr
 * @param elfun
 * @param elmats
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("MeltingBaseNLFormIntegrator::AssembleElementGrad");
  int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    // elmat.SetSize(nd);
    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    // elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      el[blk]->CalcShape(ip, Psi);
      const auto& u = *elfun[blk] * Psi;
      const auto& un = this->u_old_[blk].GetValue(Tr, ip);

      const double coef_mobi = mobility[blk].compute() * ip.weight * Tr.Weight();
      const double alpha = this->get_phase_change_at_ip(Tr, ip, blk, u, un);
      double fun_val =
          coef_mobi * alpha *
          this->compute_hessian_coefficient(interpolation_potential[blk], blk, blk,
                                            {u, un});  // this->energy_derivatives(2, Tr, ip)(u);
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));  // w'(u)*(du, psi)
    }
  }
}

/**
 * @brief Get coefficients
 * @remark Could be overridden by child classes
 *
 * @tparam VARS
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::PhaseFieldPotential, 0).has_value()) {
      interpolation_potential.add(
          *(this->get_coefficient(i, GlossaryType::PhaseFieldPotential, 0)));
    }
    // if (this->get_coefficient(i, GlossaryType::Temperature, 0).has_value()) {
    //   melting_temperature.add(*(this->get_coefficient(i, GlossaryType::Temperature, 0)));
    // }
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
double MeltingBaseNLFormIntegrator<VARS>::compute_coefficient(Coefficient coef,
                                                              const std::vector<double>& values) {
  const double u = values[0];
  const double un = values[1];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute({u});
  } else if (coef.is_explicit()) {
    coef_value = coef.compute({un});
  }
  return coef_value;
}

/**
 * @brief Return the value of the component blk of the gradient of the coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 *
 * @tparam VARS
 * @param coef
 * @param blk
 * @param values
 * @return double
 */
template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::compute_gradient_coefficient(
    Coefficient coef, const int blk, const std::vector<double>& values) {
  const double u = values[0];
  const double un = values[1];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_gradient(blk, {u});
  } else if (coef.is_explicit()) {
    coef_value = coef.compute_gradient(blk, {un});
  }
  return coef_value;
}

template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::compute_hessian_coefficient(
    Coefficient coef, const int iblk, const int jblk, const std::vector<double>& values) {
  const double u = values[0];
  const double un = values[1];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_hessian(iblk, jblk, {u});
  }
  return coef_value;
}
