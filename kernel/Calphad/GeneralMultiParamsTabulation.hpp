/**
 * @file GeneralMultiParamsTabulation.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Tabulated thermodynamic description for an binary isothermal system
 * @version 0.1
 * @date 2025-03-17
 *
 * Copyright CEA (c) 2025
 *
 */

// TO DO : PASSER EN INTERPOLATION TRILIN2AIRE, C4EST LONG MAIS LARGEMENT FAISABLE 
#include <H5Cpp.h>

#include <algorithm>
#include <boost/multi_array.hpp>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Inputs/HDF54Sloth.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

template <typename T>
class GeneralMultiParamsTabulation : CalphadBase<T> {
 private:
  std::vector<std::vector<double>> tab_p;

  boost::multi_array<double, 2> array_g;
  boost::multi_array<double, 2> array_mu;
  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(
      const size_t nb_nodes, const std::vector<T>& tp_gf,
      const std::vector<std::tuple<std::string, std::string>>& chemical_system,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  explicit GeneralMultiParamsTabulation(const Parameters& params);

  void initialize() override;

  void execute(const int dt, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>&
                   output_system) override;
  void get_1D_data_from_HDF5(const H5std_string& file_name, const std::string& dataset_name,
                             std::vector<double>& output_vector);
  void get_2D_data_from_HDF5(const H5std_string& file_name, const std::string& dataset_name,
                             boost::multi_array<double, 2>& output_multi_array);
  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~GeneralMultiParamsTabulation();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the GeneralMultiParamsTabulation object
 *
 * @tparam T
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description",
      "Analytical thermodynamic description for an ideal solution using tabulated data ");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new GeneralMultiParamsTabulation::GeneralMultiParamsTabulation object
 *
 * @param params
 */
template <typename T>
GeneralMultiParamsTabulation<T>::GeneralMultiParamsTabulation(const Parameters& params)
    : CalphadBase<T>(params) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::initialize() {
  // const H5std_string FILE_NAME("thermodata.h5"); "thermodatafromTAFID.h5"
  std::string filename = this->params_.template get_param_value_or_default<std::string>(
      "data filename", "no input file");
  const H5std_string FILE_NAME(filename);

  H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

  H5std_string DATASET_NAME;
  HDF54Sloth<2> hdf5_for_2D_table;
  DATASET_NAME = "g";
  hdf5_for_2D_table.get_data_from_HDF5(FILE_NAME, DATASET_NAME, array_g);
  DATASET_NAME = "mu";
  hdf5_for_2D_table.get_data_from_HDF5(FILE_NAME, DATASET_NAME, array_mu);
  HDF54Sloth<1> hdf5_for_vector;
  std::vector<double> tab_T;
  std::vector<double> tab_x;
  DATASET_NAME = "x";
  hdf5_for_vector.get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_x);
  DATASET_NAME = "T";
  hdf5_for_vector.get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_T);

  tab_p.emplace_back(std::move(tab_x));
  tab_p.emplace_back(std::move(tab_T));
}

/**
 * @brief Main method to calculate equilibrium states
 *
 * @tparam T
 * @param dt
 * @param aux_gf
 * @param chemical_system
 * @param output_system
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::execute(
    const int dt, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  // Clear containers and recalculation of the numbers of nodes
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  this->clear_containers();
  if (dt == 1) {
    this->check_variables_consistency(output_system);
  }
  // Thermodynamic Calculations
  this->compute(nb_nodes, tp_gf, chemical_system, output_system);

  // Use containers to update output_system
  this->update_outputs(nb_nodes, output_system);
}

/**
 * @brief Compute the CALPHAD contributions
 *
 * @tparam T
 * @param tp_gf
 * @param chemical_system
 * @param output_system
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::compute(
    const size_t nb_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  // Let us assume an ideal mixing solution
  const std::string& phase = "SOLUTION";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  std::vector<int> sorted_n_t_p = this->CU_->sort_nodes(tp_gf[0], tp_gf[1], "No", "No");

  // Process CALPHAD calculations for each node
  for (const auto& id : sorted_n_t_p) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });
    const auto loc_temperature = tp_gf_at_node[0];
    const auto& epsilon = 1.e-10;
    int i;
    int j;
    const auto& elem = std::get<0>(chemical_system[0]);
    const auto x = tp_gf_at_node[2];

    int nbtab = 2;
    std::vector<int> index(nbtab);

    std::vector<double> val;
    std::vector<double> lambda;
    val.emplace_back(std::move(x));
    val.emplace_back(std::move(loc_temperature));

    for (size_t k = 0; k < nbtab; k++) {
      auto lower_x = std::lower_bound(tab_p[k].begin(), tab_p[k].end(), val[k]);

      if (lower_x != tab_p[k].begin() && (lower_x + 1) != tab_p[k].end()) {
        index[k] = lower_x - tab_p[k].begin() - 1;
      } else if ((lower_x + 1) == tab_p[k].end()) {
        index[k] = tab_p[k].size() - 1;
      }
      double v = (val[k] - tab_p[k][index[k]]) / (tab_p[k][index[k] + 1] - tab_p[k][index[k]]);

      lambda.emplace_back(std::move(v));
    }
    double t = lambda[1];
    double u = lambda[0];

    // Energy
    for (const auto& energy_name : energy_names) {
      const auto& key = std::make_tuple(id, phase, energy_name);
      i = index[1];
      j = index[0];
      this->energies_of_phases_[key] =
          (1 - t) * (1 - u) * array_g[i][j] + t * (1 - u) * array_g[i + 1][j] +
          t * u * array_g[i + 1][j + 1] + (1 - t) * u * array_g[i][j + 1];
    }

    // Chemical potential
    this->chemical_potentials_[std::make_tuple(id, elem)] =
        (1 - t) * (1 - u) * array_mu[i][j] + t * (1 - u) * array_mu[i + 1][j] +
        t * u * array_mu[i + 1][j + 1] + (1 - t) * u * array_mu[i][j + 1];

    // Molar fraction
    const auto& key = std::make_tuple(id, phase, elem);
    this->elem_mole_fraction_by_phase_[key] = val[0];
  }
}
/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 *
 * @tparam T
 * @param output_system
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(false, "GeneralMultiParamsTabulation is only built for mu, x and g.\n");

        SlothInfo::debug("Output not available for this Calphad problem: ", output_type);
        break;
      }
    }
  }
}
/**
 * @brief Finalization actions (free memory)
 *
 * @tparam T
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
GeneralMultiParamsTabulation<T>::~GeneralMultiParamsTabulation() {}
