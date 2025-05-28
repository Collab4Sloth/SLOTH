
/**
 * @file Nucleation.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to compute inter-diffusion coefficients
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Property/PropertyBase.hpp"

#pragma once

class Nucleation : public PropertyBase {
 private:
  std::string secondary_phase_;
  bool nucleation_already_detected_{false};

 protected:
  std::vector<mfem::Vector> dgm_;
  std::vector<mfem::Vector> xph_;
  void check_variables_consistency(
      const std::vector<std::vector<std::string>>& unks_info,
      const std::vector<std::vector<std::string>>& vect_aux_infos) override;
  void get_property(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;

 public:
  explicit Nucleation(const Parameters& params);

  ~Nucleation();
};

/**
 * @brief Construct a new Nucleation::Nucleation object
 *
 */
Nucleation::Nucleation(const Parameters& params) : PropertyBase(params) {
  this->secondary_phase_ = this->params_.template get_param_value<std::string>("secondary_phase");
}

/**
 * @briefCheck consistency of variables (primary and auxiliary)
 *
 */
void Nucleation::check_variables_consistency(
    const std::vector<std::vector<std::string>>& unks_info,
    const std::vector<std::vector<std::string>>& vect_aux_infos) {
  // Check auxiliary variables
  bool secondary_phase_found = false;
  for (size_t i = 0; i < vect_aux_infos.size(); ++i) {
    const auto& variable_info = vect_aux_infos[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");

    size_t vsize = variable_info.size();
    const std::string& symbol = toLowerCase(variable_info.back());
    if (calphad_outputs::from(symbol) != calphad_outputs::dgm) continue;

    MFEM_VERIFY(variable_info.size() == 2,
                "Nucleation::check_variables_consistency(). Error while getting driving forces. "
                "Expected [name of the phase , 'dgm']");
    if (variable_info[0] == this->secondary_phase_) secondary_phase_found = true;
    break;
  }
  MFEM_VERIFY(secondary_phase_found,
              "Nucleation::check_variables_consistency(). Error while getting driving forces. "
              "Expected [secondary_phase , 'dgm']");
}

/**
 * @brief Compute the variables (properties) as functions of auxialiary variables
 *
 * @param vect_unk
 * @param unks_info
 * @param vect_aux_gf
 * @param vect_aux_infos
 */
void Nucleation::get_property(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
        output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  if (!this->nucleation_already_detected_) {
    this->dgm_.clear();
    this->xph_.clear();
    for (const auto& [aux_infos, aux_gf] : input_system) {
      if (!aux_infos.empty()) {
        const int size = aux_infos.size();
        const std::string_view type_info_view(aux_infos.back());
        if (!type_info_view.empty() && type_info_view.starts_with("dgm") &&
            toUpperCase(aux_infos[size - 2]) == toUpperCase(this->secondary_phase_)) {
          this->dgm_.emplace_back(aux_gf);
        }
        if (!type_info_view.empty() && type_info_view == "xph" &&
            toUpperCase(aux_infos[size - 2]) == toUpperCase(this->secondary_phase_)) {
          std::cout << " get prop liquid" << std::endl;

          this->xph_.emplace_back(aux_gf);
        }
      }
    }

    // const int size_gf = this->dgm_[0].Size();
    const int size_gf = this->xph_[0].Size();
    std::vector<double> dgm(size_gf, 0.);
    std::vector<double> xph(size_gf, 0.);
    for (int k = 0; k < size_gf; k++) {
      dgm[k] = this->dgm_[0](k);
      xph[k] = this->xph_[0](k);
    }

    for (auto& [output_infos, output_value] : output_system) {
      mfem::Vector vv(dgm.size());
      for (int k = 0; k < dgm.size(); k++) {
        if (dgm[k] > 0.) {
          // vv(k) = 0.;
          vv(k) = 1. - xph[k];
          // SlothInfo::print(" Nucleation at node  ", k);
          SlothInfo::print(" Nucleation at node  ", k, " with phase fraction ", xph[k]);
          this->nucleation_already_detected_ = true;
        } else {
          vv(k) = output_value.get()[k];
        }
      }
      output_value.get() = vv;
    }
  } else {
    std::cout << "NUCLEATION ALREADY FOUND" << std::endl;
  }
}

/**
 * @brief Destroy the Nucleation::Nucleation object
 *
 */
Nucleation::~Nucleation() {}
