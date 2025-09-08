
/**
 * @file InterDiffusionCoefficient.hpp
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
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Property/PropertyBase.hpp"

#pragma once

/**
 * @brief Compute inter-diffusion coefficients
 *
 * J_i = Sum_j [ M_ij * grad mu_i ] - Sum_j [ M_ij * grad mu_j ]
 *
 * with M_ij = x_i * x_j * (M_i x_j + M_j x_i)
 *
 */
class InterDiffusionCoefficient : public PropertyBase {
 private:
  // List of the chemical elements
  std::set<std::string> list_components_;

  // Flag used to know if the formula must be extended in two-phase
  bool single_phase_ = true;

 protected:
  // Chemical element corresponding to the first element
  std::string first_component_;
  // Chemical element corresponding to the last element (the reference for diffusion chemical
  // potential)
  std::string last_component_;
  // The name of the primary phase
  std::string primary_phase_;
  // The name of the secondary phase (the phase formed during phase change)
  std::string secondary_phase_;
  // The phase-field variable
  std::vector<mfem::Vector> phi_gf_;
  // The molar fraction of the elements
  std::map<std::string, mfem::Vector> x_gf_;
  // The mobilities of the elements for the primary phase
  std::map<std::string, mfem::Vector> primary_mob_gf_;
  // The mobilities of the elements for the secondary phase
  std::map<std::string, mfem::Vector> secondary_mob_gf_;

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;
  void get_property(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
          output_system,
      std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) override;

 public:
  explicit InterDiffusionCoefficient(const Parameters& params);

  ~InterDiffusionCoefficient();
};

/**
 * @brief Construct a new InterDiffusionCoefficient::InterDiffusionCoefficient object
 *
 */
InterDiffusionCoefficient::InterDiffusionCoefficient(const Parameters& params)
    : PropertyBase(params) {
  this->first_component_ =
      toUpperCase(this->params_.template get_param_value<std::string>("first_component"));
  this->last_component_ =
      toUpperCase(this->params_.template get_param_value<std::string>("last_component"));
  this->primary_phase_ =
      toUpperCase(this->params_.template get_param_value<std::string>("primary_phase"));
  this->secondary_phase_ =
      toUpperCase(this->params_.template get_param_value_or_default<std::string>("secondary_phase",
                                                                                 "UNDEFINED"));
}

/**
 * @brief Check the consistency of the inputs and outputs of the property problem.
 *
 * @param output_system The outputs property problem (primary variables).
 * @param input_system  The inputs of the property problem (auxiliary variables).
 */
void InterDiffusionCoefficient::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
        output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  // Check variables
  bool has_inter_mob = false;
  for (const auto& [aux_infos, gf] : output_system) {
    MFEM_VERIFY(!aux_infos.empty(),
                "InterDiffusionCoefficient::check_variables_consistency: error while "
                "getting auxiliary variables. Additional informations must be given for each "
                "variables.");
    const int vsize =
        aux_infos.size() - 1;  // -1 because the first info is the name of the variable
    const std::string_view symbol(toLowerCase(aux_infos.back()));

    if (symbol == "inter_mob") {
      MFEM_VERIFY(vsize == 2,
                  "InterDiffusionCoefficient::check_variables_consistency: error while "
                  "getting molar fractions.Expected [element, 'inter_mob']");
      has_inter_mob = true;
    }
  }
  MFEM_VERIFY(has_inter_mob,
              "InterDiffusionCoefficient::check_variables_consistency: error while "
              "getting inter_mob variables");

  bool has_mob = false;
  bool has_x = false;

  // Check auxiliary variables
  for (const auto& [aux_infos, aux_gf] : input_system) {
    MFEM_VERIFY(!aux_infos.empty(),
                "InterDiffusionCoefficient::check_variables_consistency: error while "
                "getting auxiliary variables. Additional informations must be given for each "
                "auxiliary variables.");
    const int vsize = aux_infos.size();
    const std::string_view symbol(toLowerCase(aux_infos.back()));

    if (symbol == "x") {
      has_x = true;
      MFEM_VERIFY(vsize == 2,
                  "InterDiffusionCoefficient::check_variables_consistency: error while "
                  "getting molar fractions.Expected [element, 'x']");
    } else if (symbol == "mob") {
      has_mob = true;
      MFEM_VERIFY(vsize == 3,
                  "InterDiffusionCoefficient::check_variables_consistency: error while "
                  "getting mobilities.Expected at least [phase, element, 'mob']");
    } else if (symbol == "phi") {
      MFEM_VERIFY(vsize == 2,
                  "InterDiffusionCoefficient::check_variables_consistency: error while "
                  "getting phase indicator.Expected at least [phase, 'phi']");

      this->single_phase_ = false;
    }
  }
  MFEM_VERIFY(has_mob,
              "InterDiffusionCoefficient::check_variables_consistency: error while "
              "getting mob auxiliary variables");

  MFEM_VERIFY(has_x,
              "InterDiffusionCoefficient::check_variables_consistency: error while "
              "getting x auxiliary variables");
}

/**
 * @brief Get the values of the property
 *
 * The values of the property are stored in the outputs (primary variables of the property
 * problem). They are calculated from the inputs (auxiliary variables of the property problem).
 *
 * @param output_system The outputs property problem (primary variables).
 * @param input_system  The inputs of the property problem (auxiliary variables).
 */
void InterDiffusionCoefficient::get_property(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
        output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  primary_mob_gf_.clear();
  secondary_mob_gf_.clear();

  x_gf_.clear();
  phi_gf_.clear();
  for (const auto& [aux_infos, aux_gf] : input_system) {
    if (!aux_infos.empty()) {
      const int size = aux_infos.size();
      const std::string_view type_info_view(aux_infos.back());

      if (!type_info_view.empty() && type_info_view.starts_with("mob")) {
        const std::string phase_info(aux_infos[size - 3]);
        const std::string elem_info(aux_infos[size - 2]);
        if (phase_info == this->primary_phase_) {
          this->primary_mob_gf_.emplace(toUpperCase(elem_info), aux_gf);

        } else if (phase_info == this->secondary_phase_) {
          this->secondary_mob_gf_.emplace(toUpperCase(elem_info), aux_gf);
        }
      }
      if (!type_info_view.empty() && type_info_view == "x") {
        const std::string elem_info(aux_infos[size - 2]);
        this->x_gf_.emplace(toUpperCase(elem_info), aux_gf);
      }
      if (!type_info_view.empty() && type_info_view == "phi") {
        const std::string elem_info(aux_infos[size - 2]);
        this->phi_gf_.emplace_back(aux_gf);
      }
    }
  }
  const int size_gf = this->primary_mob_gf_.at(this->first_component_).Size();

  auto mult_vectors = [](const std::vector<double>& v1, const std::vector<double>& v2,
                         std::vector<double>& out) {
    std::ranges::transform(v1, v2, out.begin(), [](double x, double y) { return x * y; });
  };

  auto sum_vectors = [](const std::vector<double>& v1, const std::vector<double>& v2,
                        std::vector<double>& out) {
    std::ranges::transform(v1, v2, out.begin(), [](double x, double y) { return x + y; });
  };

  auto scale_vector = [](std::vector<double>& v, double scalar) {
    std::ranges::transform(v, v.begin(), [=](double x) { return x * scalar; });
  };

  auto add_scalar = [](std::vector<double>& v, double scalar) {
    std::ranges::transform(v, v.begin(), [=](double x) { return x + scalar; });
  };

  std::vector<std::vector<double>> vmob;
  vmob.reserve(this->primary_mob_gf_.size());

  //
  std::vector<double> phi(size_gf, 1.0);
  if (!this->single_phase_) {
    for (int k = 0; k < size_gf; k++) {
      phi[k] = std::clamp(this->phi_gf_[0](k), 0.0, 1.0);
    }
  }
  std::vector<double> xo(this->x_gf_.at(this->first_component_).Size(), 0.);
  for (int k = 0; k < this->x_gf_.at(this->first_component_).Size(); k++) {
    xo[k] = this->x_gf_.at(this->first_component_)(k);
  }

  std::vector<double> mobo(this->primary_mob_gf_.at(this->first_component_).Size(), 0.);
  if (this->single_phase_) {
    for (int k = 0; k < this->primary_mob_gf_.at(this->first_component_).Size(); k++) {
      mobo[k] = this->primary_mob_gf_.at(this->first_component_)(k);
    }
  } else {
    for (int k = 0; k < this->primary_mob_gf_.at(this->first_component_).Size(); k++) {
      mobo[k] = phi[k] * this->primary_mob_gf_.at(this->first_component_)(k);

      mobo[k] += (1. - phi[k]) * this->secondary_mob_gf_.at(this->first_component_)(k);
    }
  }

  std::vector<double> mob(size_gf, 0.0);
  std::vector<double> sum_mob(size_gf, 0.0);

  for (const auto& [compo, mobi] : this->primary_mob_gf_) {
    if (compo == this->first_component_) continue;

    std::vector<double> xj(size_gf, 0.0);
    std::vector<double> mobj(size_gf, 0.0);

    if (compo == this->last_component_) {
      add_scalar(xj, -1.0);
      for (const auto& [other_compo, _] : this->primary_mob_gf_) {
        if (other_compo != this->last_component_) {
          for (int k = 0; k < this->x_gf_.at(other_compo).Size(); k++) {
            xj[k] += this->x_gf_.at(other_compo)(k);
          }
        }
      }
      scale_vector(xj, -1.0);
    } else {
      for (int k = 0; k < this->x_gf_.at(compo).Size(); k++) {
        xj[k] += this->x_gf_.at(compo)(k);
      }
    }

    mult_vectors(mobo, xj, mob);
    if (this->single_phase_) {
      for (int k = 0; k < mobi.Size(); k++) {
        mobj[k] = mobi(k) * xo[k];
      }
    } else {
      for (int k = 0; k < mobi.Size(); k++) {
        mobj[k] = (mobi(k) * phi[k] + (1. - phi[k]) * this->secondary_mob_gf_.at(compo)(k)) * xo[k];
      }
    }
    sum_vectors(mob, mobj, mob);
    mult_vectors(mob, xj, mob);
    mult_vectors(mob, xo, mob);
    sum_vectors(sum_mob, mob, sum_mob);
    scale_vector(mob, -1.0);
    vmob.emplace_back(mob);
  }
  vmob.insert(vmob.begin(), sum_mob);
  ///
  int j = 0;
  for (auto& [output_infos, output_value] : output_system) {
    mfem::Vector vv(vmob[j].size());
    for (int k = 0; k < vmob[j].size(); k++) {
      vv(k) = vmob[j][k];
      std::string msg = "Error for j " + std::to_string(j) + "  at node " + std::to_string(k) +
                        " with mob " + std::to_string(vv(k)) + " phi " + std::to_string(phi[k]);
      MFEM_VERIFY(vv(k) < 1., msg.c_str());
    }
    output_value.get() = vv;
    j++;
  }
}

/**
 * @brief Destroy the InterDiffusionCoefficient::InterDiffusionCoefficient object
 *
 */
InterDiffusionCoefficient::~InterDiffusionCoefficient() {}
