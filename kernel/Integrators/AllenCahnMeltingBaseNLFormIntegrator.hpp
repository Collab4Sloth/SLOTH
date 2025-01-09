/**
 * @file AllenCahnMeltingBaseNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <tuple>
#include <vector>

#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/AllenCahnNLFormIntegrator.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
class AllenCahnMeltingBaseNLFormIntegrator
    : public AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI> {
 private:
  PotentialFunctions<1, SCHEME, INTERPOLATION> interpolation_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, INTERPOLATION> interpolation_second_derivative_potential_;

  FType enthalpy_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                            const mfem::IntegrationPoint& ir);

 protected:
  FType energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                           const mfem::IntegrationPoint& ir) override;
  void get_parameters(const Parameters& vectr_param) override;
  virtual double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                        const mfem::IntegrationPoint& ir) = 0;

 public:
  AllenCahnMeltingBaseNLFormIntegrator(const mfem::ParGridFunction& _u_old,
                                       const Parameters& params,
                                       const std::vector<mfem::ParGridFunction>& aux_gf);
  ~AllenCahnMeltingBaseNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Overload energy derivatives to account for enthalpy of melting derivatives
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param order_derivative
 * @return FuncType
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
FType AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::energy_derivatives(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return [this, order_derivative, &Tr, &ir](const double& u) {
    return AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::energy_derivatives(order_derivative, Tr,
                                                                               ir)(u) +
           this->enthalpy_derivative(order_derivative, Tr, ir)(u);
  };
}

/**
 * @brief Enthalpy of melting derivatives
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION

 * @param order_derivative
 * @return FuncType
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
FType AllenCahnMeltingBaseNLFormIntegrator<
    SCHEME, ENERGY, MOBI, INTERPOLATION>::enthalpy_derivative(const int order_derivative,
                                                              mfem::ElementTransformation& Tr,
                                                              const mfem::IntegrationPoint& ir) {
  return FType([this, order_derivative, &Tr, &ir](const double& u) {
    const auto& un = this->u_old_.GetValue(Tr, ir);

    FType H_derivative;
    if (order_derivative == 1) {
      H_derivative = this->interpolation_first_derivative_potential_.getPotentialFunction(un);
    } else if (order_derivative == 2) {
      H_derivative = this->interpolation_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const auto& alpha = this->get_phase_change_at_ip(Tr, ir);
    const auto& h_prime = alpha * H_derivative(u);
    return h_prime;
  });
}
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
void AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(
    const Parameters& params) {
  AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::get_parameters(params);
}
/**
 * @brief Construct AllenCahnMeltingBaseNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION

 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 * @param alpha
 * @return AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnMeltingBaseNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                         const Parameters& params,
                                         const std::vector<mfem::ParGridFunction>& aux_gf)
    : AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>(u_old, params, aux_gf) {
  this->get_parameters(params);
}

/**
 * @brief Destroy the AllenCahnMeltingBaseNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION

 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI,
                                     INTERPOLATION>::~AllenCahnMeltingBaseNLFormIntegrator() {}
