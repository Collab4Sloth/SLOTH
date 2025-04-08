/**
 * @file TernaryInterDiffusionNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief ternary inter-diffusion
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
#include "Integrators/InterDiffusionNLFormIntegrator.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the VF of an inter-diffusion equation for a ternary system
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
class TernaryInterDiffusionNLFormIntegrator
    : public InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, 3> {
 protected:
  void add_interdiffusion_flux(mfem::ElementTransformation& Tr, const int nElement,
                               const mfem::IntegrationPoint& ip, const int dim) final;

 public:
  TernaryInterDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                        const Parameters& params, std::vector<VARS*> auxvars);
  ~TernaryInterDiffusionNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new TernaryInterDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::TernaryInterDiffusionNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
TernaryInterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::
    TernaryInterDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                          const Parameters& params, std::vector<VARS*> auxvars)
    : InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, 3>(u_old, params, auxvars) {}

/**
 * @brief Build the interdiffusion flux
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 * @param nElement
 * @param ip
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void TernaryInterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::add_interdiffusion_flux(
    mfem::ElementTransformation& Tr, const int nElement, const mfem::IntegrationPoint& ip,
    const int dim) {
  std::array<mfem::Vector, 3> grad_mu;
  grad_mu[0].SetSize(dim);
  grad_mu[1].SetSize(dim);
  grad_mu[2].SetSize(dim);
  std::array<double, 3> M;
  std::array<double, 3> X;

  // First (current) component
  X[0] = this->x_gf_.at(this->current_component_).GetValue(nElement, ip);
  M[0] = this->mob_gf_.at(this->current_component_).GetValue(nElement, ip);
  auto mu = SlothGridFunction(this->mu_gf_.at(this->current_component_));
  mu.GetGradient(Tr, this->gradPsi, grad_mu[0]);
  auto sum_xi = X[0];
  int i = 1;
  // All components except the first and the last ones
  for (const auto& [component, compo_gf] : this->x_gf_) {
    if (component != this->current_component_) {
      X[i] = this->x_gf_.at(component).GetValue(nElement, ip);
      M[i] = this->mob_gf_.at(component).GetValue(nElement, ip);
      auto mu = SlothGridFunction(this->mu_gf_.at(component));
      mu.GetGradient(Tr, this->gradPsi, grad_mu[i]);

      sum_xi += X[i];
      ++i;
    }
  }
  // Last component
  MFEM_VERIFY(i == this->number_of_components_ - 1,
              "Error: the last component should not be defined yet.");
  MFEM_VERIFY(sum_xi < 1 + this->x_tol_, "Error: non-physical sum of molar fraction");
  MFEM_VERIFY(-this->x_tol_ < sum_xi, "Error: non-physical sum of molar fraction");
  // To avoid unphysical fraction
  X[i] = std::max(0., std::min(1., 1. - sum_xi));

  M[i] = this->mob_gf_.at(this->last_component_).GetValue(nElement, ip);
  auto mun = SlothGridFunction(this->mu_gf_.at(this->last_component_));
  mun.GetGradient(Tr, this->gradPsi, grad_mu[i]);

  for (int j = 1; j < this->number_of_components_; j++) {
    double inter_coefficient = X[0] * X[j] * (M[0] * X[j] + M[j] * X[0]);
    if (this->interdiffusion_scaling_) {
      inter_coefficient /= (Physical::R * this->temp_gf_.at("TEMPERATURE").GetValue(nElement, ip));
    }
    this->gradMu_.Add(inter_coefficient, grad_mu[0]);
    this->gradMu_.Add(-inter_coefficient, grad_mu[j]);
  }
}

/**
 * @brief Destroy the TernaryInterDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::TernaryInterDiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
TernaryInterDiffusionNLFormIntegrator<VARS, SCHEME,
                                      DIFFU_NAME>::~TernaryInterDiffusionNLFormIntegrator() {}
