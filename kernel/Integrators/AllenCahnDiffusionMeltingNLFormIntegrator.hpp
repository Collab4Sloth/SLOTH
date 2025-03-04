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

template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
class AllenCahnDiffusionMeltingNLFormIntegrator final
    : public AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION> {
 private:
  double melting_temperature_;
  double melting_enthalpy_;
  int n;
  std::vector<mfem::ParGridFunction> aux_gf_;

 protected:
  void get_parameters();
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;

 public:
  // AllenCahnDiffusionMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old,
  //                                           const Parameters& params,
  //                                           const std::vector<mfem::ParGridFunction>& aux_gf);
  AllenCahnDiffusionMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old,
                                            const Parameters& params, std::vector<VARS*> auxvars);
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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnDiffusionMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnDiffusionMeltingNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                              const Parameters& params, std::vector<VARS*> auxvars)
    : AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                      auxvars) {
  // auto n_ = (aux_gf.size() - 2) / 3;
  // // TO DO : check if the aux_gf size is compatible
  // // MFEM_VERIFY((std::is_same<decltype(n), std::size_t>::value),
  // //             "Wrong number of auxiliary variables");

  // this->n = n_;



  // std::vector<mfem::ParGridFunction> aux_gf = this->get_aux_gf();
  // std::vector<std::vector<std::string>> inf = this->get_aux_infos();
  // std::map<std::tuple<std::string,std::string>,mfem::ParGridFunction> x;
  // std::map<std::tuple<std::string>,mfem::ParGridFunction> mu;
  // std::map<std::tuple<std::string>,mfem::ParGridFunction> g;
  // // for (size_t i = 0; i < inf.size(); i++) {
  // //   for (size_t j = 0; j < inf[i].size(); j++) {
  // //     std::cout << inf[i][j] << std::endl;
  // //   }
  // //   std::cout << "------------" << std::endl;
  // // }

  // std::vector<mfem::ParGridFunction> mu_gf;
  // for (size_t i = 0; i < inf.size(); i++) {
  //   if (inf[i][1] == "mu")
  //   {
  //     mu[std::make_tuple(inf[i][0])] = aux_gf[0];
  //   }
    
  //   for (size_t j = 0; j < inf[i].size(); j++) {
  //     std::cout << inf[i][j] << std::endl;
  //   }
  // }
  this->aux_gf_ = this->get_aux_gf();
  this->n = 2;
  this->get_parameters();
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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
void AllenCahnDiffusionMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                               INTERPOLATION>::get_parameters() {}

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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
double AllenCahnDiffusionMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  const double g_s_at_ip = this->aux_gf_[0].GetValue(Tr, ir);
  const double g_l_at_ip = this->aux_gf_[1].GetValue(Tr, ir);
  mfem::real_t driving_force = g_s_at_ip - g_l_at_ip;
  for (std::size_t i = 2; i < this->aux_gf_.size(); i += 3) {
    driving_force -= this->aux_gf_[i].GetValue(Tr, ir) * (this->aux_gf_[i + 1].GetValue(Tr, ir) -
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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnDiffusionMeltingNLFormIntegrator<
    VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::~AllenCahnDiffusionMeltingNLFormIntegrator() {}
