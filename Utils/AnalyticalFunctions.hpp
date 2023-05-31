/*
 * Copyright © CEA 2023
 *
 * \brief Analytical functions used by phase-field models
 *
 * \file AnalyticalFunctions.hpp
 * \author ci230846
 * \date 20/03/2023
 */

#pragma once
#include <any>
#include <functional>
#include <string>
#include <tuple>
#include <utility>  // std::forward
#include <vector>
#include "../Utils/PhaseFieldOptions.hpp"

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
  AnalyticalFunctions();

  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getAnalyticalFunctions(
      const std::string &analytical_function_name, Args... args);

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
    throw std::runtime_error("Not implemented");
  }
  // Sinusoide2
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide2(Args... args) {
    throw std::runtime_error("Not implemented");
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
          [center_x, a_x, radius, thickness](mfem::Vector x, double time) {
            const auto xx = a_x * (x[0] - center_x);
            const auto r = xx;
            const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
            return func;
          });
    } else {
      throw std::runtime_error(
          "multidimension_function::getHyperbolicTangent: four arguments are expected (1- "
          "center_x, "
          "2-  a_x, 3- thickness, 4- radius");
    }
  }
  // Uniform
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    const auto value = 0.;
    return std::function<double(const mfem::Vector &, double)>(
        [value](mfem::Vector x, double time) { return value; });
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
    auto v = std::vector<double>{args...};
    if (v.size() == 1) {
      const auto mult_fact = v[0];
      return std::function<double(const mfem::Vector &, double)>(
          [mult_fact](mfem::Vector x, double time) {
            const auto sinusoide = std::exp(-2. * time) * std::sin(x[0] + x[1]);
            return mult_fact * sinusoide;
          });
    } else {
      throw std::runtime_error(
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
          [mult_fact](mfem::Vector x, double time) {
            const auto u = std::exp(-2. * time) * std::sin(x[0] + x[1]);
            const auto sinusoide = u * u * u - u;
            return mult_fact * sinusoide;
          });
    } else {
      throw std::runtime_error(
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
          [center_x, center_y, a_x, a_y, radius, thickness](mfem::Vector x, double time) {
            const auto xx = a_x * (x[0] - center_x);
            const auto yy = a_y * (x[1] - center_y);
            const auto r = std::sqrt(xx * xx + yy * yy);
            const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
            return func;
          });
    } else {
      throw std::runtime_error(
          "multidimension_function::getHyperbolicTangent: six arguments are expected (1- center_x, "
          "2- center_y, 3- a_x, 4- a_y, 5- thickness, 6- radius");
    }
  }
  // Uniform
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    const auto value = 0.;
    return std::function<double(const mfem::Vector &, double)>(
        [value](mfem::Vector x, double time) { return value; });
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
    throw std::runtime_error("Not implemented");
  }

  // Sinusoide2
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getSinusoide2(Args... args) {
    throw std::runtime_error("Not implemented");
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
          [center_x, center_y, center_z, a_x, a_y, a_z, radius, thickness](mfem::Vector x,
                                                                           double time) {
            const auto xx = a_x * (x[0] - center_x);
            const auto yy = a_y * (x[1] - center_y);
            const auto zz = a_z * (x[2] - center_z);
            const auto r = std::sqrt(xx * xx + yy * yy + zz * zz);
            const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
            return func;
          });
    } else {
      throw std::runtime_error(
          "multidimension_function::getHyperbolicTangent: height arguments are expected (1- "
          "center_x, "
          "2- center_y, 3- center_z, 4- a_x, 5- a_y, 6- a_z, 7- thickness, 8- radius");
    }
  }
  // Uniform
  template <typename... Args>
  std::function<double(const mfem::Vector &, double)> getUniform(Args... args) {
    const auto value = 0.;
    return std::function<double(const mfem::Vector &, double)>(
        [value](mfem::Vector x, double time) { return value; });
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

//  explicit instantiation of AnalyticalFunctions class template/
template class AnalyticalFunctions<1>;
template class AnalyticalFunctions<2>;
template class AnalyticalFunctions<3>;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new analytical function:: analytical function object
 *
 */
template <int DIM>
AnalyticalFunctions<DIM>::AnalyticalFunctions() {}

/**
 * @brief return the function associated with the analytical_function_name
 *
 * @param analytical_function_name
 * @return const double
 */
template <int DIM>
template <class... Args>
std::function<double(const mfem::Vector &, double)>
AnalyticalFunctions<DIM>::getAnalyticalFunctions(const std::string &analytical_function_name,
                                                 Args... args) {
  switch (AnalyticalFunctionsType::from(analytical_function_name)) {
    case AnalyticalFunctionsType::Heaviside:
      return this->getHeaviside(args...);
    case AnalyticalFunctionsType::Sinusoide:
      return this->getSinusoide(args...);
    case AnalyticalFunctionsType::Sinusoide2:
      return this->getSinusoide2(args...);
    case AnalyticalFunctionsType::HyperbolicTangent:
      return this->getHyperbolicTangent(args...);
    case AnalyticalFunctionsType::Uniform:
      return this->getUniform(args...);
    default:
      throw std::runtime_error(
          "AnalyticalFunctions::getAnalyticalFunctions: Heaviside, Sinusoide, HyperbolicTangent "
          "and Uniform "
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
