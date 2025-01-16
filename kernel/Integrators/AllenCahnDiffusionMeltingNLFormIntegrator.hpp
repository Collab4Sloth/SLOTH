/**
 * @file AllenCahnDiffusionMeltingNLFormIntegrator.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief FV for the Allen Cahn equation with the driving force induced by mass transfer
 * @version 0.1
 * @date 2025-01-16
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <tuple>
#include <vector>

#include "Integrators/AllenCahnMeltingBaseNLFormIntegrator.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
class AllenCahnDiffusionMeltingNLFormIntegrator final
    : public AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION> {
 private:
  double melting_temperature_;
  double melting_enthalpy_;
  int n;

 protected:
  void get_parameters(const Parameters& vectr_param) override;
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;

 public:
  AllenCahnDiffusionMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old,
                                            const Parameters& params,
                                            const std::vector<mfem::ParGridFunction>& aux_gf);
  ~AllenCahnDiffusionMeltingNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct AllenCahnDiffusionMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param u_old
 * @param params
 * @param aux_gf
 * @return AllenCahnDiffusionMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnDiffusionMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnDiffusionMeltingNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                              const Parameters& params,
                                              const std::vector<mfem::ParGridFunction>& aux_gf)
    : AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                aux_gf) {
  auto n_ = (aux_gf.size() - 2) / 3;
  MFEM_VERIFY((std::is_same<decltype(n), std::size_t>::value),
              "Wrong number of auxiliary variables");
  this->n = n_;
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
void AllenCahnDiffusionMeltingNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(
    const Parameters& params) {
  AllenCahnMeltingBaseNLFormIntegrator<SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(params);
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
double AllenCahnDiffusionMeltingNLFormIntegrator<
    SCHEME, ENERGY, MOBI, INTERPOLATION>::get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                                                 const mfem::IntegrationPoint& ir) {
  const double g_s_at_ip = this->aux_gf_[0].GetValue(Tr, ir);
  const double g_l_at_ip = this->aux_gf_[1].GetValue(Tr, ir);
  mfem::real_t driving_force = g_s_at_ip - g_l_at_ip;
  for (std::size_t i = 2; i < this->aux_gf_.size(); i += 3) {
    driving_force += this->aux_gf_[i].GetValue(Tr, ir) * (this->aux_gf_[i + 1].GetValue(Tr, ir) -
                                                          this->aux_gf_[i + 2].GetValue(Tr, ir));
  }
  return driving_force;
}
/**
 * @brief Destroy the AAllenCahnDiffusionMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnDiffusionMeltingNLFormIntegrator<
    SCHEME, ENERGY, MOBI, INTERPOLATION>::~AllenCahnDiffusionMeltingNLFormIntegrator() {}
