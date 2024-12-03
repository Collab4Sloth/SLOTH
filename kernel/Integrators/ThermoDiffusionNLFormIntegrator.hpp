/**
 * @file ThermoDiffusionNLFormIntegrator.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief FV for the mass diffusion equation
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <memory>
#include <tuple>
#include <vector>

#include "Coefficients/DiffusionCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the FV of the mass diffusion equation
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
class ThermoDiffusionNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  const Parameters diffusion_params_;
  mfem::ParGridFunction u_old_, mu_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  mfem::DenseMatrix gradPsi, gradPsiOld;
  mfem::Vector Psi, gradU, gradUOld, gradmu, elfunOld;

  template <typename... Args>
  double diffusion(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                   const double u, const Parameters& parameters);
  template <typename... Args>
  double diffusion_prime(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                         const double u, const Parameters& parameters);

 protected:
  virtual void get_parameters(const Parameters& vectr_param);

 public:
  ThermoDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                  const std::vector<mfem::ParGridFunction>& aux_gf);
  ~ThermoDiffusionNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun, mfem::Vector& elvect);

  virtual void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun, mfem::DenseMatrix& elmat);
  
  void computeMu();

  std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>> get_energy(
      mfem::ParGridFunction* gfu, const double diffu);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::get_parameters(const Parameters& params) {
}
/**
 * @brief Return the value of the diffusion coefficient at integration point
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 * @tparam Args
 * @param Tr
 * @param ip
 * @param gfu
 * @param parameters
 * @return double
 */
template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
template <typename... Args>
double ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::diffusion(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  if (SCHEME == CoefficientDiscretization::Implicit) {
    DiffusionCoefficient<0, DIFFU_NAME> diff_coeff(u, parameters);
    return diff_coeff.Eval(Tr, ip);
  } else {
    DiffusionCoefficient<0, DIFFU_NAME> diff_coeff(&this->u_old_, parameters);
    return diff_coeff.Eval(Tr, ip);
  }
}

/**
 * @brief Return the first derivative of the diffusion coefficient vs the variable at integration
 * point
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 * @tparam Args
 * @param Tr
 * @param ip
 * @param gfu
 * @param parameters
 * @return double
 */
template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
template <typename... Args>
double ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::diffusion_prime(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  double coef = 0.;
  if (SCHEME == CoefficientDiscretization::Implicit) {
    DiffusionCoefficient<1, DIFFU_NAME> diff_coeff(u, parameters);
    coef = diff_coeff.Eval(Tr, ip);
  }
  return coef;
}

/**
 * @brief Construct a new ThermoDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::ThermoDiffusionNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::ThermoDiffusionNLFormIntegrator(
    const mfem::ParGridFunction& u_old, const Parameters& params,
    const std::vector<mfem::ParGridFunction>& aux_gf)
    : u_old_(u_old), diffusion_params_(params), aux_gf_(aux_gf), mu_(u_old) {
  // this->get_parameters(params);
}

/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */

template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::computeMu() {
    for (int i = 0; i < this->u_old_.Size(); i++)
    {
      this->mu_[i] = std::log(std::max(this->u_old_[i],1e-10)) + 1;
      this->mu_[i] = this->u_old_[i];
    }
    }


/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */

template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  int nd = el.GetDof();
  int dim = el.GetDim();
  int nElement = Tr.ElementNo;

  this->gradPsi.SetSize(nd, dim);
  this->Psi.SetSize(nd);
  this->gradPsiOld.SetSize(nd, dim);
  this->gradU.SetSize(dim);
  this->gradUOld.SetSize(dim);
  elvect.SetSize(nd);
  this->gradmu.SetSize(dim);
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());
  elvect = 0.0;
  this->computeMu();
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    this->mu_.GetGradient(Tr,this->gradmu);

    el.CalcShape(ip, Psi);
    this->u_old_.GetElementDofValues(nElement, this->elfunOld);
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * Psi;

  
    el.CalcPhysDShape(Tr, this->gradPsi);  // Compute derivative of shape function

    this->gradPsi.MultTranspose(elfun, this->gradU);
    this->gradPsi.MultTranspose(this->elfunOld, this->gradUOld);   // possible de passer par GetGradient pour pouvoir eviter de definir elfunUOld

    // std::cout  << std::fixed << std::setprecision(15) << "Element n° = " << nElement << " var_u_old = " << this->u_old_.GetValue(Tr,ip) << "    vector_u_old * Psi  = "  << elfunOld * Psi  << "   u = " << u  << std::endl;
    // std::cout  << std::fixed << std::setprecision(15) << "Element n° = " << nElement << "     " <<  this->gradUOld[0] << "     " << this->gradU[0] << "   " << this->gradUOld[0] - this->gradU[0] <<  std::endl;
    // std::cout << "  uold   " << this->u_old_ << std::endl;
    // std::cout << "  mu   "<< this->mu_ << std::endl;
    // std::cout << std::fixed << std::setprecision(15) << "   grad_mu" << gradmu[0] << "   gradUOld  " << this->gradUOld[0] << "   gradU  " << this->gradU[0] << std::endl;
    const double coeff_diffu =
        this->diffusion(Tr, ip, u, this->diffusion_params_) * ip.weight * Tr.Weight();

    //this->gradU.Add(-1.0,this->gradUOld);
    this->gradU *= 1e-8* ip.weight * Tr.Weight();
    this->gradUOld *= 1e-8* ip.weight * Tr.Weight();
    this->gradmu *= coeff_diffu;

    std::cout << this->gradmu[0] << "   " << this->gradUOld[0] << "   " << this->gradU[0] << std::endl;
    this->gradPsi.AddMult(this->gradU, elvect,1.0);
    this->gradPsi.AddMult(this->gradUOld, elvect,-1.0);
    this->gradPsi.AddMult(this->gradmu, elvect,1.0);
    // this->gradPsi.AddMult(this->gradU, elvect,1.0);
  }
}

/**
 * @brief  Jacobian part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::DenseMatrix& elmat) {
  int nd = el.GetDof();
  int dim = el.GetDim();

  Psi.SetSize(nd);
  gradPsi.SetSize(nd, dim);
  elmat.SetSize(nd);
  mfem::Vector vec;
  vec.SetSize(nd);

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  elmat = 0.0;
  vec = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcShape(ip, Psi);
    const auto& u = elfun * Psi;
    // Laplacian : compute (grad(phi), grad(psi)), phi is shape function.

    Tr.SetIntPoint(&ip);
    const double coeff_diffu =
        this->diffusion(Tr, ip, u, this->diffusion_params_) * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, elmat);

    gradPsi.MultTranspose(elfun, gradU);
    gradPsi.AddMult(gradU, vec);
    const auto coef_diffu_derivative = this->diffusion_prime(Tr, ip, u, this->diffusion_params_);
    AddMult_a_VWt(coef_diffu_derivative * ip.weight * Tr.Weight(), Psi, vec, elmat);
  }
}

template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::get_energy(mfem::ParGridFunction* gfu,
                                                                const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu, diffu);
}

/**
 * @brief Destroy the ThermoDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::ThermoDiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */

template <CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
ThermoDiffusionNLFormIntegrator<SCHEME, DIFFU_NAME>::~ThermoDiffusionNLFormIntegrator() {}
