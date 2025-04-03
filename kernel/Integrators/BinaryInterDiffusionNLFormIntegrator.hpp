/**
 * @file BinaryInterDiffusionNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief binary inter-diffusion
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
 * @brief  Class dedicated to the VF of an inter-diffusion equation for a binary system
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
class BinaryInterDiffusionNLFormIntegrator
    : public InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, 2> {
 protected:
  void add_interdiffusion_flux(mfem::ElementTransformation& Tr, const int nElement,
                               const mfem::IntegrationPoint& ip, const int dim) final;

 public:
  BinaryInterDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                       std::vector<VARS*> auxvars);
  ~BinaryInterDiffusionNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new BinaryInterDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::BinaryInterDiffusionNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
BinaryInterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::
    BinaryInterDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                         const Parameters& params, std::vector<VARS*> auxvars)
    : InterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME, 2>(u_old, params, auxvars) {}

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
void BinaryInterDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::add_interdiffusion_flux(
    mfem::ElementTransformation& Tr, const int nElement, const mfem::IntegrationPoint& ip,
    const int dim) {
  std::array<mfem::Vector, 2> grad_mu;
  grad_mu[0].SetSize(dim);
  grad_mu[1].SetSize(dim);
  std::array<double, 2> M;

  auto calculate_gradient = [this, &Tr, &nElement, &ip, &dim, &grad_mu, &M](const auto& it,
                                                                            const int i) {
    const std::string& elem = it->first;
    M[i] = this->mob_gf_.at(elem).GetValue(nElement, ip);
    auto mu = SlothGridFunction(this->mu_gf_.at(elem));

    mu.GetGradient(Tr, this->gradPsi, grad_mu[i]);
  };

  auto itx = this->x_gf_.begin();
  const double x1 = itx->second.GetValue(nElement, ip);
  auto it = this->mu_gf_.begin();

  calculate_gradient(it, 0);
  ++it;
  calculate_gradient(it, 1);
  double M12 = x1 * (1.0 - x1) * (M[0] * (1.0 - x1) + M[1] * x1);
  if (this->interdiffusion_scaling_) {
    M12 /= (Physical::R * this->temp_gf_.at("TEMPERATURE").GetValue(nElement, ip));
  }
  this->gradMu_.Add(M12, grad_mu[0]);
  this->gradMu_.Add(-M12, grad_mu[1]);
}

/**
 * @brief Destroy the BinaryInterDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::BinaryInterDiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
BinaryInterDiffusionNLFormIntegrator<VARS, SCHEME,
                                     DIFFU_NAME>::~BinaryInterDiffusionNLFormIntegrator() {}
