/**
 * @file DiffusionCoeffients.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Diffusion coefficient and derivatives used by diffusion models
 * @version 0.1
 * @date 2024-06-28
 *
 * @copyright Copyright (c) 2024
 *
 */
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "Utils/PhaseFieldOptions.hpp"
// TODO(cci) faire en sorte qu'il y ait une convention (xn tjs puis des parametres)

template <int ORDER, DiffusionCoefficientDiscretization SCHEME>
struct diffusion_function {};

template <int ORDER, DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
class DiffusionFunctions {
 private:
  template <class... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    diffusion_function<ORDER, SCHEME> func;
    return func.getLinear(args...);
  }

 public:
  DiffusionFunctions();

  template <class... Args>
  std::function<double(const double&)> getFunction(Args... args);

  ~DiffusionFunctions();
};

//////////////////////////////////////////////////////
// IMPLICIT, EXPLICIT, SEMI_IMPLICIT
//////////////////////////////////////////////////////
///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct diffusion_function<0, DiffusionCoefficientDiscretization::Implicit> {
  /**
   * @brief  Get the  Implicit Linear coefficient
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() > 2) {
      const auto xn = v[0];
      const auto alpha = v[1];
      const auto kappa = v[2];
      return std::function<double(const double&)>([alpha, kappa](double x) {
        const auto coeff = alpha * x + kappa;
        return coeff;
      });
    } else {
      throw std::runtime_error(
          "diffusion_function::getLinear: two arguments are expected for implicit scheme");
    }
  }
};

template <>
struct diffusion_function<0, DiffusionCoefficientDiscretization::Explicit> {
  /**
   * @brief  Get the  Explicit Linear coefficient
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() > 2) {
      const auto xn = v[0];
      const auto alpha = v[1];
      const auto kappa = v[2];
      return std::function<double(const double&)>([xn, alpha, kappa](double x) {
        const auto coeff = alpha * xn + kappa;
        return coeff;
      });
    } else {
      throw std::runtime_error(
          "diffusion_function::getLinear: three arguments are expected for explicit scheme");
    }
  }
};

template <>
struct diffusion_function<0, DiffusionCoefficientDiscretization::SemiImplicit> {
  /**
   * @brief  Get the  SemiExplicit Linear coefficient
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    throw std::runtime_error("diffusion_function::getLinear: semi-explicit scheme not available");
  }
};

///////////////////////
// ORDER = 1
///////////////////////
template <>
struct diffusion_function<1, DiffusionCoefficientDiscretization::Implicit> {
  /**
   * @brief  Get the First derivative of Linear coefficient
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() > 2) {
      const auto xn = v[0];
      const auto alpha = v[1];
      const auto kappa = v[2];
      return std::function<double(const double&)>([alpha](double x) { return alpha; });
    } else {
      throw std::runtime_error(
          "diffusion_function::getLinear: only one argument is expected for implicit/explicit "
          "schemes");
    }
  }
};
template <>
struct diffusion_function<1, DiffusionCoefficientDiscretization::Explicit> {
  /**
   * @brief  Get the First derivative of Linear coefficient
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() > 2) {
      const auto xn = v[0];
      const auto alpha = v[1];
      const auto kappa = v[2];
      return std::function<double(const double&)>([alpha](double x) { return 0.; });
    } else {
      throw std::runtime_error(
          "diffusion_function::getLinear: only one argument is expected for implicit/explicit "
          "schemes");
    }
  }
};
///////////////////////
// ORDER  2
///////////////////////
template <DiffusionCoefficientDiscretization SCHEME>
struct diffusion_function<2, SCHEME> {
  /**
   * @brief Get the Second derivative of Linear coefficient
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLinear(Args... args) {
    return std::function<double(const double&)>([](double x) { return 0.; });
  }
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Diffusion Functions<ORDER, SCHEME, COEFFICIENT>:: Diffusion Functions
 * object
 *
 * @tparam ORDER
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */
template <int ORDER, DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
DiffusionFunctions<ORDER, SCHEME, COEFFICIENT>::DiffusionFunctions() {}

/**
 * @brief return the function associated with the COEFFICIENT, its ORDER of
 * derivative and its time discretization scheme
 *
 * @tparam ORDER
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @tparam Args
 * @param args
 * @return std::function<double(const double&)>
 */
template <int ORDER, DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
template <class... Args>
std::function<double(const double&)> DiffusionFunctions<ORDER, SCHEME, COEFFICIENT>::getFunction(
    Args... args) {
  switch (COEFFICIENT) {
    case DiffusionCoefficients::Linear:
      return this->getLinear(args...);
    default:
      throw std::runtime_error("DiffusionFunctions::getFunctions: linear function is available");
      break;
  }
}

/**
 * @brief Destroy the Diffusion Functions<ORDER, SCHEME, COEFFICIENT>::Diffusion Functions object
 *
 * @tparam ORDER
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */
template <int ORDER, DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
DiffusionFunctions<ORDER, SCHEME, COEFFICIENT>::~DiffusionFunctions() {}
