/**
 * @file AnalyticalSolidBinary.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Analytical thermodynamic description for binary system in solid phase (ideal case)
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
 *
 */
#include <functional>
#include <string>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

template <typename T>
class AnalyticalSolidBinary : CalphadBase<T> {
 public:
  explicit AnalyticalSolidBinary(const Parameters& params);

  void initialize() override;

  void execute(const int dt, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>&
                   output_system) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~AnalyticalSolidBinary();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the AnalyticalSolidBinary object
 *
 * @tparam T
 */
template <typename T>
void AnalyticalSolidBinary<T>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description", "Analytical thermodynamic description for a binary system in a solid phase. ");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new AnalyticalSolidBinary::AnalyticalSolidBinary object
 *
 * @param params
 */
template <typename T>
AnalyticalSolidBinary<T>::AnalyticalSolidBinary(const Parameters& params) : CalphadBase<T>(params) {
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void AnalyticalSolidBinary<T>::initialize() {}

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
void AnalyticalSolidBinary<T>::execute(
    const int dt, const std::vector<T>& aux_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {}

/**
 * @brief Finalization actions (free memory)
 *
 * @tparam T
 */
template <typename T>
void AnalyticalSolidBinary<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
AnalyticalSolidBinary<T>::~AnalyticalSolidBinary() {}
