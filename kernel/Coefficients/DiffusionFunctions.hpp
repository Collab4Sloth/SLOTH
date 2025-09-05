/**
 * @file DiffusionFunctions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Diffusion functions
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
 * 
 * This file is part of SLOTH.
 * 
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
