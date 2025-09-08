/**
 * @file CahnHilliardNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the CahnHilliard equations
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
#include <utility>
#include <vector>

#include "Coefficients/LambdaCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/OmegaCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the Allen Cahn equation
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
class CahnHilliardNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator,
                                     public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;

  template <typename... Args>
  double mobility(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                  const Parameters& parameters);

  template <typename... Args>
  double omega(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
               const Parameters& parameters);

  template <typename... Args>
  double lambda(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                const Parameters& parameters);

  FType double_well_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                               const mfem::IntegrationPoint& ir);

  void check_variables_consistency();

 protected:
  std::vector<mfem::ParGridFunction> u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  std::vector<mfem::Vector> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_gf_infos_;
  std::vector<mfem::ParGridFunction> temp_gf_;
  bool scale_mobility_by_temperature_{false};

  virtual FType energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                                   const mfem::IntegrationPoint& ir);

 public:
  CahnHilliardNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                               const Parameters& params, std::vector<VARS*> auxvars);
  ~CahnHilliardNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvec);

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmats);

  std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>> get_energy(
      std::vector<mfem::ParGridFunction*> gfu, const double omega);
  std::unique_ptr<GradientEnergyCoefficient> get_grad_energy(
      std::vector<mfem::ParGridFunction*> gfu, const double lambda);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Energy derivatives contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return FuncType
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::energy_derivatives(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return [this, order_derivative, &Tr, &ir](const double& u) {
    return this->double_well_derivative(order_derivative, Tr, ir)(u);
  };
}

/**
 * @brief Return the value of the mobility coefficient at integration point
 *
 * @remark Actually, only explicit mobility is available. (See diffusion integrator  to extend
 * towards implicit mobility. )
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::mobility(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  const int nElement = Tr.ElementNo;
  // TODO(cci) return vector of double
  MobilityCoefficient<0, MOBI> mobi_coeff(&this->u_old_[0], parameters);
  double mob_coeff = mobi_coeff.Eval(Tr, ip);

  return mob_coeff;
}

/**
 * @brief Return the value of the omega coefficient at integration point
 *
 * @remark Actually, only constant omega is available.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::omega(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  OmegaCoefficient<0, Omega::Constant> omega_coeff(&this->u_old_[0], parameters);
  return omega_coeff.Eval(Tr, ip);
}
/**
 * @brief Return the value of the lambda coefficient at integration point
 *
 * @remark Actually, only constant lambda is available.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::lambda(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  LambdaCoefficient<0, Lambda::Constant> lambda_coeff(&this->u_old_[0], parameters);
  return lambda_coeff.Eval(Tr, ip);
}

/**
 * @brief double_well_derivative contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return std::function<double(const double&, const double&)>
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::double_well_derivative(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return FType([this, order_derivative, &Tr, &ir](const double& u) {
    const auto& un = this->u_old_[0].GetValue(Tr, ir);
    FType W_derivative;
    if (order_derivative == 1) {
      W_derivative = this->energy_first_derivative_potential_.getPotentialFunction(un);
    } else if (order_derivative == 2) {
      W_derivative = this->energy_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const double omega = this->omega(Tr, ir, u, this->params_);

    const auto& w_prime = omega * W_derivative(u);
    return w_prime;
  });
}

/**
 * @brief Construct a new CahnHilliardNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::CahnHilliardNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->check_variables_consistency();
}

/**
 * @brief Check variables consistency
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::check_variables_consistency() {
  this->aux_gf_ = this->get_aux_gf();
  this->aux_old_gf_ = this->get_aux_old_gf();
  this->aux_gf_infos_ = this->get_aux_infos();

  // Temperature scaling for mobility
  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>: at least "
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
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
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

      const double xx = ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      const double lambda = this->lambda(Tr, ip, phi, this->params_);

      gradU *= -xx * lambda;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // Given u, compute (w'(u), psi), psi is shape function

      const double ww = xx * (mu - this->energy_derivatives(1, Tr, ip)(phi));

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

      const auto& mu = *elfun[blk] * Psi;
      const auto& phi = *elfun[off_blk] * Psi;

      const double coef_mob = this->mobility(Tr, ip, phi, this->params_) * ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= coef_mob;
      gradPsi.AddMult(gradU, *elvect[blk]);
    }
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // loop over diagonal entries
  int num_blocks = el.Size();

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

      const double xx = -ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);

      const double lambda = this->lambda(Tr, ip, phi, this->params_);

      AddMult_a_AAt(xx * lambda, gradPsi, *elmats(blk, blk));

      double fun_val = xx * this->energy_derivatives(2, Tr, ip)(phi);
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

      const auto& phi = *elfun[off_blk] * Psi;
      const double coef_mob = this->mobility(Tr, ip, phi, this->params_) * ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);

      AddMult_a_AAt(coef_mob, gradPsi, *elmats(blk, blk));
    }
  }
}

/**
 * @brief Return the energy coefficient associated with the integrator
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param gfu
 * @param lambda
 * @param omega
 * @return EnergyCoefficient<SCHEME, ENERGY>
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>>
CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::get_energy(
    std::vector<mfem::ParGridFunction*> gfu, const double omega) {
  return std::make_unique<HomogeneousEnergyCoefficient<ENERGY>>(gfu[0], omega);
}

/**
 * @brief  Return the gradient energy coefficient associated with the integrator
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param gfu
 * @param lambda
 * @return GradientEnergyCoefficient
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<GradientEnergyCoefficient>
CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::get_grad_energy(
    std::vector<mfem::ParGridFunction*> gfu, const double lambda) {
  return std::make_unique<GradientEnergyCoefficient>(gfu[0], lambda);
}

/**
 * @brief Destroy the CahnHilliardNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
CahnHilliardNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::~CahnHilliardNLFormIntegrator() {}
