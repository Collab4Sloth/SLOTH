/**
 * @file AnalyticalFunctions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief List of Analytical functions used by phase-field models
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#pragma once
#include <any>
#include <functional>
#include <string>
#include <utility>  // std::forward
#include <vector>

#include "Options/Options.hpp"
#include "Utils/Utils.hpp"

template <int DIM>
struct multidimension_function {};

template <int DIM>
class AnalyticalFunctions {
 private:
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHeaviside(Args... args) {
    multidimension_function<DIM> func;
    return func.getHeaviside(args...);
  }

  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHyperbolicTangent(Args... args) {
    multidimension_function<DIM> func;
    return func.getHyperbolicTangent(args...);
  }

  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getParabolic(Args... args) {
    multidimension_function<DIM> func;
    return func.getParabolic(args...);
  }

  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    multidimension_function<DIM> func;
    return func.getUniform(args...);
  }
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide(Args... args) {
    multidimension_function<DIM> func;
    return func.getSinusoide(args...);
  }
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide2(Args... args) {
    multidimension_function<DIM> func;
    return func.getSinusoide2(args...);
  }

 public:
  template <class... Args>
  AnalyticalFunctions(AnalyticalFunctionsType::value function_name, Args... args_func);

  explicit AnalyticalFunctions(std::function<double(const mfem::Vector &, double)> user_function);

  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getAnalyticalFunctions(
      AnalyticalFunctionsType::value function_name, Args... args);

  std::function<double(const mfem::Vector &, double)> analytical_function_;

  std::function<double(const mfem::Vector &, double)> getFunction() const;

  ~AnalyticalFunctions();
};

///////////////////////
// DIM = 1
///////////////////////
template <>
struct multidimension_function<1> {
  // Heaviside
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHeaviside(Args... args) {
    return std::function<double(const mfem::Vector &, double)>([](mfem::Vector x, double time) {
      if (x[0] < 0.5) {
        return 0.0;
      } else {
        return 1.0;
      }
    });
  }
  // Sinusoide
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide(Args... args) {
    mfem::mfem_error("Not implemented");
  }
  // Sinusoide2
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide2(Args... args) {
    mfem::mfem_error("Not implemented");
  }

  // TANH
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHyperbolicTangent(Args... args) {
    auto v = std::vector<double>{args...};
    if (v.size() == 4) {
      const auto center_x = v[0];
      const auto a_x = v[1];
      const auto thickness = v[2];
      const auto radius = v[3];

      return std::function<double(const mfem::Vector &, double)>(
          [center_x, a_x, radius, thickness](const mfem::Vector &x, double time) {
            const auto r = a_x * (x[0] - center_x);
            const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
            return func;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getHyperbolicTangent: four arguments are expected (1- "
          "center_x, 2-  a_x, 3- thickness, 4- radius");
    }
  }

  // PARABOLIC
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getParabolic(Args... args) {
    auto v = std::vector<double>{args...};
    if (v.size() == 4) {
      const auto rmax = v[0];
      const auto fo = v[1];
      const auto lin_pow = v[2];
      const auto cond = v[3];

      return std::function<double(const mfem::Vector &, double)>(
          [rmax, fo, cond, lin_pow](const mfem::Vector &x, double time) {
            const auto r2 = x[0] * x[0];
            const auto rmax2 = rmax * rmax;

            const auto func = fo + lin_pow * (rmax2 - r2) / (4.0 * M_PI * cond * rmax2);
            return func;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getParabolic: three arguments are expected (1- "
          "maximum radius, 2-  initial temperature, 3- lineic power, 4- conductivity");
    }
  }

  // Uniform
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto value = v[0];
      return std::function<double(const mfem::Vector &, double)>(
          [value](const mfem::Vector &x, double time) { return value; });
    } else {
      mfem::mfem_error("multidimension_function::getUniform: only one argument is expected");
    }
  }
};
///////////////////////
// DIM = 2
///////////////////////
template <>
struct multidimension_function<2> {
  // Heaviside
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHeaviside(Args... args) {
    return std::function<double(const mfem::Vector &, double)>(
        [](const mfem::Vector &x, double time) {
          if (x[0] < 0.5) {
            return 0.0;
          } else {
            return 1.0;
          }
        });
  }
  // Sinusoide
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide(Args... args) {
    auto v = std::vector<double>{args...};
    if (v.size() == 1) {
      const auto mult_fact = v[0];
      return std::function<double(const mfem::Vector &, double)>(
          [mult_fact](const mfem::Vector &x, double time) {
            const auto sinusoide = std::exp(-2. * time) * std::sin(x[0] + x[1]);
            return mult_fact * sinusoide;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getSinusoide: one argument is expected : multiplicative "
          "factor");
    }
  }

  // Sinusoide2
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide2(Args... args) {
    auto v = std::vector<double>{args...};
    if (v.size() == 1) {
      const auto mult_fact = v[0];
      return std::function<double(const mfem::Vector &, double)>(
          [mult_fact](const mfem::Vector &x, double time) {
            const auto u = std::exp(-2. * time) * std::sin(x[0] + x[1]);
            const auto sinusoide = u * u * u - u;
            return mult_fact * sinusoide;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getSinusoide2: one argument is expected : multiplicative "
          "factor");
    }
  }

  // TANH
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHyperbolicTangent(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 6) {
      const auto center_x = v[0];
      const auto center_y = v[1];
      const auto a_x = v[2];
      const auto a_y = v[3];
      const auto thickness = v[4];
      const auto radius = v[5];

      return std::function<double(const mfem::Vector &, double)>(
          [center_x, center_y, a_x, a_y, radius, thickness](const mfem::Vector &x, double time) {
            const auto xx = a_x * (x[0] - center_x);
            const auto yy = a_y * (x[1] - center_y);
            const auto r = std::sqrt(xx * xx + yy * yy);
            const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
            return func;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getHyperbolicTangent: six arguments are expected (1- center_x, "
          "2- center_y, 3- a_x, 4- a_y, 5- thickness, 6- radius");
    }
  }
  // PARABOLIC
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getParabolic(Args... args) {
    // mfem::mfem_error("Not implemented");
    auto v = std::vector<double>{args...};
    if (v.size() == 4) {
      const auto rmax = v[0];
      const auto fo = v[1];
      const auto lin_pow = v[2];
      const auto cond = v[3];

      return std::function<double(const mfem::Vector &, double)>(
          [rmax, fo, cond, lin_pow](const mfem::Vector &x, double time) {
            const auto r2 = x[0] * x[0];
            const auto rmax2 = rmax * rmax;

            const auto func = fo + lin_pow * (rmax2 - r2) / (4.0 * M_PI * cond * rmax2);
            return func;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getParabolic: three arguments are expected (1- "
          "maximum radius, 2-  initial temperature, 3- lineic power, 4- conductivity");
    }
  }

  // Uniform
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto value = v[0];
      return std::function<double(const mfem::Vector &, double)>(
          [value](const mfem::Vector &x, double time) { return value; });
    } else {
      mfem::mfem_error("multidimension_function::getUniform: only one argument is expected");
    }
  }
};
///////////////////////
// DIM = 3
///////////////////////
template <>
struct multidimension_function<3> {
  // Heaviside
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHeaviside(Args... args) {
    return std::function<double(const mfem::Vector &, double)>(
        [](const mfem::Vector &x, double time) {
          if (x[0] < 0.5) {
            return 0.0;
          } else {
            return 1.0;
          }
        });
  }
  // Sinusoide
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide(Args... args) {
    mfem::mfem_error("Not implemented");
  }

  // Sinusoide2
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide2(Args... args) {
    mfem::mfem_error("Not implemented");
  }

  // TANH
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getHyperbolicTangent(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 8) {
      const auto center_x = v[0];
      const auto center_y = v[1];
      const auto center_z = v[2];
      const auto a_x = v[3];
      const auto a_y = v[4];
      const auto a_z = v[5];
      const auto thickness = v[6];
      const auto radius = v[7];

      return std::function<double(const mfem::Vector &, double)>(
          [center_x, center_y, center_z, a_x, a_y, a_z, radius, thickness](const mfem::Vector &x,
                                                                           double time) {
            const auto xx = a_x * (x[0] - center_x);
            const auto yy = a_y * (x[1] - center_y);
            const auto zz = a_z * (x[2] - center_z);
            const auto r = std::sqrt(xx * xx + yy * yy + zz * zz);
            const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
            return func;
          });
    } else {
      mfem::mfem_error(
          "multidimension_function::getHyperbolicTangent: height arguments are expected (1- "
          "center_x, "
          "2- center_y, 3- center_z, 4- a_x, 5- a_y, 6- a_z, 7- thickness, 8- radius");
    }
  }
  // PARABOLIC
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getParabolic(Args... args) {
    mfem::mfem_error("Not implemented");
  }
  // Uniform
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto value = v[0];
      return std::function<double(const mfem::Vector &, double)>(
          [value](const mfem::Vector &x, double time) { return value; });
    } else {
      mfem::mfem_error("multidimension_function::getUniform: only one argument is expected");
    }
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new analytical function:: analytical function object
 *
 */
template <int DIM>
template <class... Args>
AnalyticalFunctions<DIM>::AnalyticalFunctions(AnalyticalFunctionsType::value function_name,
                                              Args... args_func) {
  // std::apply(
  //     [function_name, this](Args... args) {
  //   this->analytical_function_ = this->getAnalyticalFunctions(function_name, args...);
  // },
  // args_func);
  this->analytical_function_ = this->getAnalyticalFunctions(function_name, args_func...);
}

template <int DIM>
AnalyticalFunctions<DIM>::AnalyticalFunctions(
    std::function<double(const mfem::Vector &, double)> user_function) {
  this->analytical_function_ = user_function;
}

/**
 * @brief Return the analytical function
 *
 * @tparam DIM
 * @tparam Args
 * @return std::function<double(const mfem::Vector &, double)>
 */
template <int DIM>
std::function<double(const mfem::Vector &, double)> AnalyticalFunctions<DIM>::getFunction() const {
  return this->analytical_function_;
}

/**
 * @brief return the function associated with the analytical_function_name
 *
 * @tparam DIM
 * @tparam Args
 * @param function_name
 * @param args
 * @return std::function<double(const mfem::Vector &, double)>
 */
template <int DIM>
template <class... Args>
std::function<double(const mfem::Vector &, double)>
AnalyticalFunctions<DIM>::getAnalyticalFunctions(AnalyticalFunctionsType::value function_name,
                                                 Args... args) {
  switch (function_name) {
    case AnalyticalFunctionsType::Heaviside:
      return this->getHeaviside(args...);
    case AnalyticalFunctionsType::Sinusoide:
      return this->getSinusoide(args...);
    case AnalyticalFunctionsType::Sinusoide2:
      return this->getSinusoide2(args...);
    case AnalyticalFunctionsType::HyperbolicTangent:
      return this->getHyperbolicTangent(args...);
    case AnalyticalFunctionsType::Parabolic:
      return this->getParabolic(args...);
    case AnalyticalFunctionsType::Uniform:
      return this->getUniform(args...);
    default:
      mfem::mfem_error(
          "AnalyticalFunctions::getAnalyticalFunctions: Heaviside, Sinusoide, HyperbolicTangent, "
          "Parabolic, Uniform "
          "analytical function  are available");
      break;
  }
}

/**
 * @brief Destroy the analytical function :: analytical function  object
 *
 */
template <int DIM>
AnalyticalFunctions<DIM>::~AnalyticalFunctions() {}
