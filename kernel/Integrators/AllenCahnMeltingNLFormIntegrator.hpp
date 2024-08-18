/**
 * @file AllenCahnMeltingNLFormIntegrator.hpp
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

#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Coefficients/PhaseFieldMobilities.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/AllenCahnNLFormIntegrator.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION, PhaseChange PHASECHANGE>
class AllenCahnMeltingNLFormIntegrator final
    : public AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI> {
 private:
  double alpha_;

  MobilityFunctions<MOBI> mobility_function_;
  PhaseChangeFunction<PHASECHANGE> phase_change_function_;

  PotentialFunctions<1, SCHEME, INTERPOLATION> interpolation_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, INTERPOLATION> interpolation_second_derivative_potential_;

  FuncType enthalpy_derivative(const int order_derivative);

  FuncType energy_derivatives(const int order_derivative) override;

 protected:
  void get_parameters(const Parameters& vectr_param) override;

 public:
  AllenCahnMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old, const Parameters& params);
  ~AllenCahnMeltingNLFormIntegrator();
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
 * @tparam PHASECHANGE
 * @param order_derivative
 * @return FuncType
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION, PhaseChange PHASECHANGE>
FuncType
AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION,
                                 PHASECHANGE>::energy_derivatives(const int order_derivative) {
  return [this, order_derivative](const double& u, const double& un) {
    return AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::energy_derivatives(order_derivative)(
               u, un) +
           this->enthalpy_derivative(order_derivative)(u, un);
  };
}

/**
 * @brief Enthalpy of melting derivatives
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @tparam PHASECHANGE
 * @param order_derivative
 * @return FuncType
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION, PhaseChange PHASECHANGE>
FuncType
AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION,
                                 PHASECHANGE>::enthalpy_derivative(const int order_derivative) {
  return FuncType([this, order_derivative](const double& u, const double& un) {
    std::function<double(const double&)> H_derivative;
    if (order_derivative == 1) {
      H_derivative = this->interpolation_first_derivative_potential_.getPotentialFunction(un);
    } else if (order_derivative == 2) {
      H_derivative = this->interpolation_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const auto& Mphi = this->mobility_function_.getMobilityFunction(un);
    const auto& alpha = this->phase_change_function_.getPhaseChangeFunction(un);
    const auto& h_prime = Mphi(this->mob_) * alpha(this->alpha_) * H_derivative(u);
    return h_prime;
  });
}
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION, PhaseChange PHASECHANGE>
void AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION,
                                      PHASECHANGE>::get_parameters(const Parameters& params) {
  this->alpha_ = params.get_param_value_or_default<double>("melting_factor", 0.);
}
/**
 * @brief Construct AllenCahnMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @tparam PHASECHANGE
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 * @param alpha
 * @return AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION, PHASECHANGE>::
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION, PhaseChange PHASECHANGE>
AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION, PHASECHANGE>::
    AllenCahnMeltingNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params)
    : AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>(u_old, params) {
  AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::get_parameters(params);
  this->get_parameters(params);
}

/**
 * @brief Destroy the AllenCahnMeltingNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @tparam PHASECHANGE
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION, PhaseChange PHASECHANGE>
AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION,
                                 PHASECHANGE>::~AllenCahnMeltingNLFormIntegrator() {}
