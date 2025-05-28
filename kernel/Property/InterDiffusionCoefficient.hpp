
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
#include "Profiling/Profiling.hpp"
#include "Property/PropertyBase.hpp"

#pragma once

class InterDiffusionCoefficient : public PropertyBase {
 private:
  std::set<std::string> list_components_;
  std::map<std::string, int> mob_index_;
  std::map<std::string, int> x_index_;
  void check_nan(mfem::Vector vv, std::string name);

 protected:
  std::string first_component_;
  std::string last_component_;
  std::vector<mfem::Vector> phi_gf_;
  std::map<std::string, mfem::Vector> x_gf_;
  std::map<std::string, mfem::Vector> mob_gf_;
  void check_variables_consistency(
      const std::vector<std::vector<std::string>>& unks_info,
      const std::vector<std::vector<std::string>>& vect_aux_infos) override;
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
}

/**
 * @briefCheck consistency of variables (primary and auxiliary)
 *
 */
void InterDiffusionCoefficient::check_variables_consistency(
    const std::vector<std::vector<std::string>>& unks_info,
    const std::vector<std::vector<std::string>>& vect_aux_infos) {
  // // Check variables
  // for (size_t i = 0; i < unks_info.size(); ++i) {
  //   const auto &variable_info = unks_info[i];
  //   MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
  //   const std::string &symbol = toLowerCase(variable_info.back());
  // }
  // Check auxiliary variables
  for (size_t i = 0; i < vect_aux_infos.size(); ++i) {
    const auto& variable_info = vect_aux_infos[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");

    size_t vsize = variable_info.size();
    const std::string& symbol = toLowerCase(variable_info[vsize - 2]);
    if (symbol == "x") {
      MFEM_VERIFY(variable_info.size() == 3,
                  "InterDiffusionCoefficient::check_variables_consistency: error while "
                  "getting molar fractions.Expected [element, 'x']");
      const std::string elem_info = variable_info[variable_info.size() - 3];
      this->x_index_.emplace(toUpperCase(elem_info), i);
    } else if (symbol == "mob") {
      MFEM_VERIFY(variable_info.size() >= 2,
                  "InterDiffusionCoefficient::check_variables_consistency: error while "
                  "getting mobilities.Expected at least [element, 'mob']");
      const std::string elem_info = variable_info[variable_info.size() - 3];
      this->list_components_.insert(toUpperCase(elem_info));
      this->mob_index_.emplace(toUpperCase(elem_info), i);
    }
  }
}

/**
 * @brief Compute the variables (properties) as functions of auxialiary variables
 *
 * @param vect_unk
 * @param unks_info
 * @param vect_aux_gf
 * @param vect_aux_infos
 */
void InterDiffusionCoefficient::get_property(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>&
        output_system,
    std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system) {
  //
  // M12+M13 grad mu1 - M12 grad mu2 -M13 grad mu3
  // -M12 grad mu1 + M12 + M23 grad mu2 -M23 grad mu3
  // -M13 grad mu1 - M23 grad mu2 +M13+M23 grad mu3
  //
  mob_gf_.clear();
  x_gf_.clear();
  phi_gf_.clear();
  for (const auto& [aux_infos, aux_gf] : input_system) {
    if (!aux_infos.empty()) {
      const int size = aux_infos.size();
      const std::string_view type_info_view(aux_infos.back());

      if (!type_info_view.empty() && type_info_view.starts_with("mob")) {
        const std::string elem_info(aux_infos[size - 2]);
        this->mob_gf_.emplace(toUpperCase(elem_info), aux_gf);
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
  const int size_gf = this->mob_gf_.at(this->first_component_).Size();

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
  vmob.reserve(this->mob_gf_.size());

  std::vector<double> xo(this->x_gf_.at(this->first_component_).Size(), 0.);
  for (int k = 0; k < this->x_gf_.at(this->first_component_).Size(); k++) {
    xo[k] = this->x_gf_.at(this->first_component_)(k);
  }
  std::vector<double> mobo(this->mob_gf_.at(this->first_component_).Size(), 0.);
  for (int k = 0; k < this->mob_gf_.at(this->first_component_).Size(); k++) {
    mobo[k] = this->mob_gf_.at(this->first_component_)(k);
  }
  std::vector<double> mob(size_gf, 0.0);
  std::vector<double> sum_mob(size_gf, 0.0);

  // Two phase case (temporary)
  std::vector<double> phi(size_gf, 0.);
  if (!this->phi_gf_.empty()) {
    for (int k = 0; k < size_gf; k++) {
      phi[k] = this->phi_gf_[0](k);
    }
    for (auto& [compo, mobi] : this->mob_gf_) {
      const auto mobi_max = mobi.Max();
      for (int k = 0; k < mobi.Size(); k++) {
        if (mobi(k) < 0.) {
          // No solid detected
          mobi(k) = mobi_max;
        } else {
          std::cout << " Two phase detected " << k << " compo " << compo << " : before " << mobi(k)
                    << " phi " << phi[k] << " mob max " << mobi_max << std::endl;
          mobi(k) = mobi(k) * phi[k] + (1. - phi[k]) * mobi_max;
          std::cout << " Two phase detected : after " << mobi(k) << std::endl;
        }
      }
    }
  }

  for (const auto& [compo, mobi] : this->mob_gf_) {
    if (compo == this->first_component_) continue;

    std::vector<double> xj(size_gf, 0.0);
    std::vector<double> mobj(size_gf, 0.0);

    if (compo == this->last_component_) {
      add_scalar(xj, -1.0);
      for (const auto& [other_compo, _] : this->mob_gf_) {
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
    for (int k = 0; k < mobi.Size(); k++) {
      mobj[k] = mobi(k) * xo[k];
    }
    sum_vectors(mob, mobj, mob);
    mult_vectors(mob, xj, mob);
    mult_vectors(mob, xo, mob);
    sum_vectors(sum_mob, mob, sum_mob);
    scale_vector(mob, -1.0);
    vmob.emplace_back(mob);
  }
  vmob.insert(vmob.begin(), sum_mob);

  int j = 0;
  for (auto& [output_infos, output_value] : output_system) {
    mfem::Vector vv(vmob[j].size());
    for (int k = 0; k < vmob[j].size(); k++) {
      vv(k) = vmob[j][k];
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
