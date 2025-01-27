/**
 * @file OmegaFunctions.hpp
 * @author  ci230846  (clement.introini@cea.fr)
 * @brief Omega functions
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

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"

template <int ORDER>
struct Omega_function {};

///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct Omega_function<0> {
  /**
   * @brief Get the Constant Omega fonction
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    auto a0 = params.get_param_value<double>("omega");
    return FType([a0](double x) { return a0; });
  }
};
///////////////////////
// ORDER  = 1
///////////////////////
template <>
struct Omega_function<1> {
  /**
   * @brief Get the Constant Omega function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    return FType([](double x) { return 0.; });
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <int ORDER, Omega NAME>
class OmegaFunctions {
 private:
  FType getConstant(const Parameters& params) {
    Omega_function<ORDER> func;
    return func.getConstant(params);
  }

 public:
  OmegaFunctions();

  FType getFunction(const Parameters& params);

  ~OmegaFunctions();
};

/**
 * @brief Construct a new OmegaFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Omega NAME>
OmegaFunctions<ORDER, NAME>::OmegaFunctions() {}

/**
 * @brief Return the property function depending on the order of derivative and its name
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param args
 * @return FType
 */
template <int ORDER, Omega NAME>
FType OmegaFunctions<ORDER, NAME>::getFunction(const Parameters& params) {
  switch (NAME) {
    case Omega::Constant:
      return this->getConstant(params);
    default:
      mfem::mfem_error("OmegaFunctions::getFunction: Constant is available");
      break;
  }
}

/**
 * @brief Destroy the Omega  Functions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Omega NAME>
OmegaFunctions<ORDER, NAME>::~OmegaFunctions() {}
