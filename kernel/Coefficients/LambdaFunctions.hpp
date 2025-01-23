/**
 * @file LambdaFunctions.hpp
 * @author  ci230846  (clement.introini@cea.fr)
 * @brief Lambda functions
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
struct Lambda_function {};

///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct Lambda_function<0> {
  /**
   * @brief Get the Constant Lambda fonction
   *
   * @param params
   * @return FType
   */
  FType getConstant(const Parameters& params) {
    auto a0 = params.get_param_value<double>("lambda");
    return FType([a0](double x) { return a0; });
  }
};
///////////////////////
// ORDER  = 1
///////////////////////
template <>
struct Lambda_function<1> {
  /**
   * @brief Get the Constant Lambda function (1st derivative)
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

template <int ORDER, Lambda NAME>
class LambdaFunctions {
 private:
  FType getConstant(const Parameters& params) {
    Lambda_function<ORDER> func;
    return func.getConstant(params);
  }

 public:
  LambdaFunctions();

  FType getFunction(const Parameters& params);

  ~LambdaFunctions();
};

/**
 * @brief Construct a new LambdaFunctions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Lambda NAME>
LambdaFunctions<ORDER, NAME>::LambdaFunctions() {}

/**
 * @brief Return the property function depending on the order of derivative and its name
 *
 * @tparam ORDER
 * @tparam NAME
 * @tparam Args
 * @param args
 * @return FType
 */
template <int ORDER, Lambda NAME>
FType LambdaFunctions<ORDER, NAME>::getFunction(const Parameters& params) {
  switch (NAME) {
    case Lambda::Constant:
      return this->getConstant(params);
    default:
      mfem::mfem_error("LambdaFunctions::getFunction: Constant is available");
      break;
  }
}

/**
 * @brief Destroy the Lambda  Functions object
 *
 * @tparam ORDER
 * @tparam NAME
 */
template <int ORDER, Lambda NAME>
LambdaFunctions<ORDER, NAME>::~LambdaFunctions() {}
