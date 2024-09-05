/**
 * @file MobilityFunctions.hpp
 * @author  ci230846  (clement.introini@cea.fr)
 * @brief Mobility functions
 * @version 0.1
 * @date 2024-09-04
 *
 * @copyright Copyright (c) 2024
 *
 */
#pragma once
#include <algorithm>
#include <functional>
#include <string>
#include <vector>

#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"

template <int ORDER>
struct Mobility_function {};

///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct Mobility_function<0> {
  /**
   * @brief Get the Constant Mobility fonction
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    auto a0 = params.get_param_value<double>("mob");
    return FType([a0](double x) { return a0; });
  }
  /**
   * @brief Get the Degenerated Mobility function
   *
   * @param params
   * @return FType
   */
  FType getDegenerated(const Parameters& params) {
    auto a0 = params.get_param_value<double>("mob");
    auto a1 = params.get_param_value<double>("power");

    return FType([a0, a1](double x) {
      const double eps = 1.e-6;
      const auto xx = std::min(1. - eps, std::max(eps, x));
      return a0 * std::pow(xx, a1) * std::pow(1. - xx, a1);
    });
  }
};
///////////////////////
// ORDER  = 1
///////////////////////
template <>
struct Mobility_function<1> {
  /**
   * @brief Get the Constant Mobility function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    return FType([](double x) { return 0.; });
  }
  /**
   * @brief Get the Degenerated Mobility function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getDegenerated(const Parameters& params) {
    auto a0 = params.get_param_value<double>("mob");
    auto a1 = params.get_param_value<double>("power");

    return FType([a0, a1](double x) {
      const double eps = 1.e-6;
      const auto xx = std::min(1. - eps, std::max(eps, x));
      return a0 * a1 * std::pow(xx, a1 - 1.) * std::pow(1. - xx, a1 - 1.) * (1. - 2. * xx);
    });
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <int ORDER, Mobility NAME>
class MobilityFunctions {
 private:
  FType getDegenerated(const Parameters& params) {
    Mobility_function<ORDER> func;
    return func.getDegenerated(params);
  }
  FType getConstant(const Parameters& params) {
    Mobility_function<ORDER> func;
    return func.getConstant(params);
  }

 public:
  MobilityFunctions();

  FType getFunction(const Parameters& params);

  ~MobilityFunctions();
};

/**
 * @brief Construct a new MobilityFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Mobility NAME>
MobilityFunctions<ORDER, NAME>::MobilityFunctions() {}

/**
 * @brief Return the property function depending on the order of derivative and its name
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param args
 * @return FType
 */
template <int ORDER, Mobility NAME>
FType MobilityFunctions<ORDER, NAME>::getFunction(const Parameters& params) {
  switch (NAME) {
    case Mobility::Constant:
      return this->getConstant(params);
    case Mobility::Degenerated:
      return this->getDegenerated(params);
    default:
      mfem::mfem_error("MobilityFunctions::getFunction: Linear are available");
      break;
  }
}

/**
 * @brief Destroy the Mobility  Functions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Mobility NAME>
MobilityFunctions<ORDER, NAME>::~MobilityFunctions() {}
