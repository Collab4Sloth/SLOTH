/**
 * @file AllenCahnConstantMeltingNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-9-3
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <tuple>
#include <vector>

#include "Integrators/AllenCahnMeltingBaseNLFormIntegrator.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
class AllenCahnConstantMeltingNLFormIntegrator final
    : public AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION> {
 private:
  double alpha_;

 protected:
  void get_parameters(const Parameters& vectr_param) override;
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;

 public:
  AllenCahnConstantMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old,
                                           const Parameters& params,
                                           const std::vector<mfem::ParGridFunction>& aux_gf);
  ~AllenCahnConstantMeltingNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct AllenCahnConstantMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param u_old
 * @param params
 * @param aux_gf
 * @return AllenCahnConstantMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnConstantMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnConstantMeltingNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                             const Parameters& params,
                                             const std::vector<mfem::ParGridFunction>& aux_gf)
    : AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                aux_gf) {
  this->get_parameters(params);
}

/**
 * @brief Get parameters
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param params
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
void AllenCahnConstantMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(
    const Parameters& params) {
  AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(params);
  this->alpha_ = params.get_param_value<double>("melting_factor");
}

/**
 * @brief Get the value of the phae change at integration point
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param Tr
 * @param ir
 * @return double
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
double AllenCahnConstantMeltingNLFormIntegrator<
    SCHEME, ENERGY, MOBI, INTERPOLATION>::get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                                                 const mfem::IntegrationPoint& ir) {
  const double phase_change_at_ip = this->alpha_;
  return phase_change_at_ip;
}
/**
 * @brief Destroy the AAllenCahnConstantMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnConstantMeltingNLFormIntegrator<
    SCHEME, ENERGY, MOBI, INTERPOLATION>::~AllenCahnConstantMeltingNLFormIntegrator() {}
