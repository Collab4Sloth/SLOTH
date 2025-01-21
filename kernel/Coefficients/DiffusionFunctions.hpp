/**
 * @file DiffusionFunctions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Diffusion functions
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
struct Diffusion_function {};

///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct Diffusion_function<0> {
  /**
   * @brief Get the Constant Diffusion fonction
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    auto a0 = params.get_param_value<double>("D");
    return FType([a0](double x) { return a0; });
  }
  /**
   * @brief Get the Linear Diffusion function
   *
   * @param params
   * @return FType
   */
  FType getLinear(const Parameters& params) {
    auto a0 = params.get_param_value<double>("D_0");
    auto a1 = params.get_param_value<double>("D_1");
    return FType([a0, a1](double x) { return a0 + a1 * x; });
  }
  /**
   * @brief Get the Constant Diffusion fonction
   *
   * @param params
   * @return FType
   */
  FType getLog(const Parameters& params) {
    auto a0 = params.get_param_value<double>("D");
    return FType([a0](double x) { return (a0 * (std::log(std::max(x, 1e-10)) + 1.)); });
  }
};
///////////////////////
// ORDER  = 1
///////////////////////
template <>
struct Diffusion_function<1> {
  /**
   * @brief Get the Constant Diffusion function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    return FType([](double x) { return 0.; });
  }
  /**
   * @brief Get the Linear Diffusion function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getLinear(const Parameters& params) {
    auto a1 = params.get_param_value<double>("D_1");
    return FType([a1](double x) { return a1; });
  }
  /**
   * @brief Get the Log Diffusion function (1st derivative)
   *
   * @param params
   * @return FType
   */
  FType getLog(const Parameters& params) {
    auto a1 = params.get_param_value<double>("D");
    return FType([a1](double x) { return std::min(a1 / x, 1e-5); });
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <int ORDER, Diffusion NAME>
class DiffusionFunctions {
 private:
  FType getLinear(const Parameters& params) {
    Diffusion_function<ORDER> func;
    return func.getLinear(params);
  }
  FType getConstant(const Parameters& params) {
    Diffusion_function<ORDER> func;
    return func.getConstant(params);
  }
  FType getLog(const Parameters& params) {
    Diffusion_function<ORDER> func;
    return func.getLog(params);
  }

 public:
  DiffusionFunctions();

  FType getFunction(const Parameters& params);

  ~DiffusionFunctions();
};

/**
 * @brief Construct a new DiffusionFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Diffusion NAME>
DiffusionFunctions<ORDER, NAME>::DiffusionFunctions() {}

/**
 * @brief Return the property function depending on the order of derivative and its name
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param args
 * @return FType
 */
template <int ORDER, Diffusion NAME>
FType DiffusionFunctions<ORDER, NAME>::getFunction(const Parameters& params) {
  switch (NAME) {
    case Diffusion::Constant:
      return this->getConstant(params);
    case Diffusion::Linear:
      return this->getLinear(params);
    case Diffusion::Log:
      return this->getLog(params);
    default:
      mfem::mfem_error("DiffusionFunctions::getFunction: Linear are available");
      break;
  }
}

/**
 * @brief Destroy the DiffusionFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Diffusion NAME>
DiffusionFunctions<ORDER, NAME>::~DiffusionFunctions() {}
