/**
 * @file InterDiffusionNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief inter-diffusion integrator
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/DiffusionCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the VF of an inter-diffusion equation
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
class InterDiffusionNLFormIntegrator : public mfem::NonlinearFormIntegrator,
                                       public SlothNLFormIntegrator<VARS> {
 private:
  std::vector<std::tuple<std::string, double>> inter_diffusion_coeff_;
  double coeff_stab_;
  std::set<std::string> list_components_;
  void get_parameters();

 protected:
  const double x_tol_ = 1.e-6;
  std::string current_component_;
  std::string last_component_;
  int number_of_components_{NBCOMPONENT};
  SlothGridFunction u_old_;
  std::map<std::string, mfem::ParGridFunction> mu_gf_;
  std::map<std::string, mfem::ParGridFunction> x_gf_;
  std::map<std::string, mfem::ParGridFunction> mob_gf_;
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradMu_;
  virtual void add_interdiffusion_flux(mfem::ElementTransformation& Tr, const int nElement,
                                       const mfem::IntegrationPoint& ip, const int dim) = 0;

 public:
  InterDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                 std::vector<VARS*> auxvars);
  ~InterDiffusionNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun, mfem::Vector& elvect);

  virtual void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun, mfem::DenseMatrix& elmat);

  std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>> get_energy(
      mfem::ParGridFunction* gfu, const double diffu);
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Get parameters
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
void InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, NBCOMPONENT>::get_parameters() {
  std::set<std::string> mob_components, x_components;
  // Get the last component : x_last = 1 - Sum[i!=last] x_i
  this->last_component_ = this->params_.template get_param_value<std::string>("last_component");
  x_components.insert(toUpperCase(this->last_component_));

  // Get stabilization coefficient
  this->coeff_stab_ = this->params_.template get_param_value<double>("D");

  // Number of auxiliary variables for a n-ary system:
  // - chemical potentials : n (expected symbol "mu")
  // - mobilities : n if not defined as parameters named "Ma", "Mb",.... Priority to auxialiary
  //   variables (expected symbol "mob")
  // - molar fractions : n-2 because the primary variable is implicitly
  //   taken into account (expected symbol "x") knowing last_component_
  // Identification of the components must be uniform between these variables. Consequently, list of
  // components coming from chemical potentials must equal to the list of components coming from
  // mobilities

  // Get Chemical potentials and try to get Mobilities from Variables
  std::vector<mfem::ParGridFunction> aux_gf = this->get_aux_gf();
  std::vector<std::vector<std::string>> aux_infos = this->get_aux_infos();

  for (size_t i = 0; i < aux_gf.size(); ++i) {
    if (!aux_infos[i].empty()) {
      const auto& info = aux_infos[i];
      const int size = info.size();
      if (size >= 3) {
        const std::string_view type_info_view(info[size - 2]);
        const std::string elem_info(info[size - 3]);
        if (!elem_info.empty() && !type_info_view.empty() && type_info_view.starts_with("mu")) {
          this->mu_gf_.emplace(toUpperCase(elem_info), std::move(aux_gf[i]));
          this->list_components_.insert(toUpperCase(elem_info));
        }
        if (!elem_info.empty() && !type_info_view.empty() && type_info_view.starts_with("mob")) {
          this->mob_gf_.emplace(toUpperCase(elem_info), std::move(aux_gf[i]));
          mob_components.insert(toUpperCase(elem_info));
        }
        if (!elem_info.empty() && !type_info_view.empty() && type_info_view.starts_with("x")) {
          this->x_gf_.emplace(toUpperCase(elem_info), std::move(aux_gf[i]));
          x_components.insert(toUpperCase(elem_info));
        }
      }
    }
  }

  // Get Mobilities from Parameters if not
  if (this->mob_gf_.empty()) {
    const auto mugf = this->mu_gf_.begin()->second;
    mfem::ParGridFunction mob_gf(mugf.ParFESpace());
    for (const std::string& component : this->list_components_) {
      std::string upper_name = "M" + component;
      std::string lower_name = "M" + toLowerCase(component);
      std::string param_name;
      if (this->params_.has_parameter(upper_name)) {
        param_name = upper_name;
      } else if (this->params_.has_parameter(lower_name)) {
        param_name = lower_name;
      } else {
        std::string error_mess = "Mobility parameter for component " + component +
                                 " must be defined. Please check the data.";
        mfem::mfem_error(error_mess.c_str());
      }
      const double mob = this->params_.template get_param_value<double>(param_name);
      mfem::ConstantCoefficient cc(mob);
      mob_gf.ProjectCoefficient(cc);
      this->mob_gf_.emplace(toUpperCase(component), std::move(mob_gf));
      mob_components.insert(toUpperCase(component));
    }
  }
  std::string msg_pot = "InterDiffusionNLFormIntegrator requires " +
                        std::to_string(this->number_of_components_) +
                        " chemical potentials "
                        "among auxiliary variables";
  MFEM_VERIFY(this->mu_gf_.size() == this->number_of_components_, msg_pot.c_str());
  std::string msg_mob = "InterDiffusionNLFormIntegrator requires " +
                        std::to_string(this->number_of_components_) +
                        " mobilities either fully defined "
                        "from auxiliary variables or with Parameters";
  MFEM_VERIFY(this->mob_gf_.size() == this->number_of_components_, msg_mob.c_str());

  MFEM_VERIFY(mob_components == this->list_components_,
              "List of components for chemical potentials and mobilities must be the same. Please "
              "check your data.");

  std::set<std::string> result;
  std::set_difference(this->list_components_.begin(), this->list_components_.end(),
                      x_components.begin(), x_components.end(),
                      std::inserter(result, result.end()));

  MFEM_VERIFY(x_components.size() == this->list_components_.size() - 1 && result.size() == 1,
              "List of components for molar fraction must  those defined for chemical "
              "potentials except one. Please check your data.");
  this->current_component_ = *result.begin();
  this->x_gf_.emplace(*result.begin(), std::move(this->u_old_));
  x_components.insert(*result.begin());
  std::string msg_gf = "InterDiffusionNLFormIntegrator requires " +
                       std::to_string(this->number_of_components_ - 1) +
                       " grid function for molar fractions ";
  MFEM_VERIFY(this->x_gf_.size() == this->number_of_components_ - 1, msg_gf.c_str());
}

/**
 * @brief Construct a new InterDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::InterDiffusionNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, NBCOMPONENT>::
    InterDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                   std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->get_parameters();
}

/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
void InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, NBCOMPONENT>::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  int nd = el.GetDof();
  int dim = el.GetDim();
  int nElement = Tr.ElementNo;

  this->Psi.SetSize(nd);
  this->gradPsi.SetSize(nd, dim);

  this->gradMu_.SetSize(dim);

  elvect.SetSize(nd);
  mfem::Vector grad_uold;
  grad_uold.SetSize(dim);

  // Initialization
  elvect = 0.0;

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  // Loop over integration points
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);

    el.CalcShape(ip, Psi);
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * Psi;

    // Stabilization contribution : D_stab * (Grad u - Grad un)
    el.CalcPhysDShape(Tr, this->gradPsi);
    this->gradPsi.MultTranspose(elfun, this->gradMu_);
    this->u_old_.GetGradient(Tr, this->gradPsi, grad_uold);

    this->gradMu_.Add(-1, grad_uold);
    this->gradMu_ *= this->coeff_stab_;

    // Interdiffusion flux (see child classes)
    this->add_interdiffusion_flux(Tr, nElement, ip, dim);

    this->gradMu_ *= ip.weight * Tr.Weight();

    this->gradPsi.AddMult(this->gradMu_, elvect, 1.0);
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
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
void InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, NBCOMPONENT>::AssembleElementGrad(
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

    Tr.SetIntPoint(&ip);
    const double coeff_diffu = this->coeff_stab_ * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, elmat);
  }
}

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, NBCOMPONENT>::get_energy(
    mfem::ParGridFunction* gfu, const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu, diffu);
}

/**
 * @brief Destroy the InterDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::InterDiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME, int NBCOMPONENT>
InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME,
                               NBCOMPONENT>::~InterDiffusionNLFormIntegrator() {}
