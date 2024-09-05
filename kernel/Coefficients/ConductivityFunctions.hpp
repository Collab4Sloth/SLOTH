/**
 * @file ConductivityFunctions.hpp
 * @author  ci230846  (clement.introini@cea.fr)
 * @brief Conductivity functions
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
struct Conductivity_function {};

///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct Conductivity_function<0> {
  /**
   * @brief Get the Constant Conductivity fonction
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    auto a0 = params.get_param_value<double>("lambda");
    return FType([a0](double x) { return a0; });
  }
  /**
   * @brief Get the Linear Conductivity function
   *
   * @param params
   * @return FType
   */
  FType getLinear(const Parameters& params) {
    auto a0 = params.get_param_value<double>("lambda_0");
    auto a1 = params.get_param_value<double>("lambda_1");
    return FType([a0, a1](double x) { return a0 + a1 * x; });
  }
};
///////////////////////
// ORDER  = 1
///////////////////////
template <>
struct Conductivity_function<1> {
  /**
   * @brief Get the Constant Conductivity function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    return FType([](double x) { return 0.; });
  }
  /**
   * @brief Get the Linear Conductivity function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getLinear(const Parameters& params) {
    auto a1 = params.get_param_value<double>("cp_1");
    return FType([a1](double x) { return a1; });
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <int ORDER, Conductivity NAME>
class ConductivityFunctions {
 private:
  FType getLinear(const Parameters& params) {
    Conductivity_function<ORDER> func;
    return func.getLinear(params);
  }
  FType getConstant(const Parameters& params) {
    Conductivity_function<ORDER> func;
    return func.getConstant(params);
  }

 public:
  ConductivityFunctions();

  FType getFunction(const Parameters& params);

  ~ConductivityFunctions();
};

/**
 * @brief Construct a new ConductivityFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Conductivity NAME>
ConductivityFunctions<ORDER, NAME>::ConductivityFunctions() {}

/**
 * @brief Return the property function depending on the order of derivative and its name
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param args
 * @return FType
 */
template <int ORDER, Conductivity NAME>
FType ConductivityFunctions<ORDER, NAME>::getFunction(const Parameters& params) {
  switch (NAME) {
    case Conductivity::Constant:
      return this->getConstant(params);
    case Conductivity::Linear:
      return this->getLinear(params);
    default:
      mfem::mfem_error("ConductivityFunctions::getFunction: Linear are available");
      break;
  }
}

/**
 * @brief Destroy the ConductivityFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Conductivity NAME>
ConductivityFunctions<ORDER, NAME>::~ConductivityFunctions() {}
