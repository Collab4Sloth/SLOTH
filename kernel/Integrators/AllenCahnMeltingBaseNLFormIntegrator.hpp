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

template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
class AllenCahnMeltingBaseNLFormIntegrator
    : public AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI> {
 private:
  PotentialFunctions<1, SCHEME, INTERPOLATION> interpolation_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, INTERPOLATION> interpolation_second_derivative_potential_;

  FType enthalpy_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                            const mfem::IntegrationPoint& ir);

 protected:
  FType energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                           const mfem::IntegrationPoint& ir) override;
  virtual double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                        const mfem::IntegrationPoint& ir) = 0;
  virtual double get_seed_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir);

 public:
  AllenCahnMeltingBaseNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                       const Parameters& params, std::vector<VARS*> auxvars);
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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
FType AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                       const mfem::IntegrationPoint& ir) {
  return [this, order_derivative, &Tr, &ir](const double& u) {
    return AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::energy_derivatives(
               order_derivative, Tr, ir)(u) +
           this->enthalpy_derivative(order_derivative, Tr, ir)(u);
  };
}

/**
 * @brief Get the value of the seed used to initiate the phase transformation
 * @remark by default, 0. Child classes may re-implemented it.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param Tr
 * @param ir
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
double
AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::get_seed_at_ip(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return 0.0;
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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
FType AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    enthalpy_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                        const mfem::IntegrationPoint& ir) {
  return FType([this, order_derivative, &Tr, &ir](const double& u) {
    const auto& un = this->u_old_[0].GetValue(Tr, ir);

    const auto& alpha = this->get_phase_change_at_ip(Tr, ir);
    double seed = 0.0;

    FType H_derivative;
    if (order_derivative == 1) {
      H_derivative = this->interpolation_first_derivative_potential_.getPotentialFunction(un);
      // Explicit term, not present in Jacobian (order two for H_derivative)
      seed = this->get_seed_at_ip(Tr, ir);

    } else if (order_derivative == 2) {
      H_derivative = this->interpolation_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const auto& h_prime = alpha * H_derivative(u) + seed;
    return h_prime;
  });
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
 * @return AllenCahnMeltingBaseNLFormIntegrator<VARS,SCHEME, ENERGY, MOBI, INTERPOLATION>::
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnMeltingBaseNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                         const Parameters& params, std::vector<VARS*> auxvars)
    : AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>(u_old, params, auxvars) {}

/**
 * @brief Destroy the AllenCahnMeltingBaseNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION

 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                     INTERPOLATION>::~AllenCahnMeltingBaseNLFormIntegrator() {}
