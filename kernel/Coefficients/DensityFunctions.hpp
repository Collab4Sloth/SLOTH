/**
 * @file DensityFunctions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Density functions
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
struct density_function {};

///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct density_function<0> {
  /**
   * @brief Get the Constant Density fonction
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    auto a0 = params.get_param_value<double>("rho");
    return FType([a0](double x) { return a0; });
  }
  /**
   * @brief Get the Linear Density function
   *
   * @param params
   * @return FType
   */
  FType getLinear(const Parameters& params) {
    auto a0 = params.get_param_value<double>("rho_0");
    auto a1 = params.get_param_value<double>("rho_1");
    return FType([a0, a1](double x) { return a0 + a1 * x; });
  }
};
///////////////////////
// ORDER  = 1
///////////////////////
template <>
struct density_function<1> {
  /**
   * @brief Get the Constant Density function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    return FType([](double x) { return 0.; });
  }
  /**
   * @brief Get the Linear Density function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getLinear(const Parameters& params) {
    auto a1 = params.get_param_value<double>("rho_1");
    return FType([a1](double x) { return a1; });
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <int ORDER, Density NAME>
class DensityFunctions {
 private:
  FType getLinear(const Parameters& params) {
    density_function<ORDER> func;
    return func.getLinear(params);
  }
  FType getConstant(const Parameters& params) {
    density_function<ORDER> func;
    return func.getConstant(params);
  }

 public:
  DensityFunctions();

  FType getFunction(const Parameters& params);

  ~DensityFunctions();
};

/**
 * @brief Construct a new DensityFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Density NAME>
DensityFunctions<ORDER, NAME>::DensityFunctions() {}

/**
 * @brief Return the property function depending on the order of derivative and its name
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param args
 * @return FType
 */
template <int ORDER, Density NAME>
FType DensityFunctions<ORDER, NAME>::getFunction(const Parameters& params) {
  switch (NAME) {
    case Density::Constant:
      return this->getConstant(params);
    case Density::Linear:
      return this->getLinear(params);
    default:
      mfem::mfem_error("DensityFunctions::getFunction: Linear are available");
      break;
  }
}

/**
 * @brief Destroy the DensityFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Density NAME>
DensityFunctions<ORDER, NAME>::~DensityFunctions() {}
