/**
 *
 * Copyright CEA (C) 2026
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"
#include "Options/PhysicalPropertiesOptions.hpp"

#pragma once

/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = c**3*dvol*(6*c**2 - 15*c + 10) + c*dsurf*(1 - c) + dgb*(-eta1**2 - eta2**2 + (eta1 +
 * eta2)**2) + dvap*(-c**3*(6*c**2 - 15*c + 10) + 1)
 */
class D : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  D() : prefactor_(1.0) {}
  explicit D(const double prefactor) : prefactor_(prefactor) {}
  virtual ~D() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const
 * std::span<const double>&, const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
D::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = std::clamp(input_vector[0], 0.0, 1.0);
    double eta1 = auxiliary_vector[0];
    double eta2 = auxiliary_vector[1];
    double dvol = 0.04;
    double dvap = 0.002;
    double dsurf = 16.0;
    double dgb = 1.6;
    double F = std::pow(c, 3) * dvol * (6 * std::pow(c, 2) - 15 * c + 10) + c * dsurf * (1 - c) +
               dgb * (-std::pow(eta1, 2) - std::pow(eta2, 2) + std::pow(eta1 + eta2, 2)) +
               dvap * (-std::pow(c, 3) * (6 * std::pow(c, 2) - 15 * c + 10) + 1);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
D::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = auxiliary_vector[0];
    double eta2 = auxiliary_vector[1];
    double dvol = 0.04;
    double dvap = 0.002;
    double dsurf = 16.0;
    double dgb = 1.6;
    std::vector<double> gradient(1);
    gradient[0] =
        this->prefactor_ * (std::pow(c, 3) * dvol * (12 * c - 15) +
                            3 * std::pow(c, 2) * dvol * (6 * std::pow(c, 2) - 15 * c + 10) -
                            c * dsurf + dsurf * (1 - c) +
                            dvap * (-std::pow(c, 3) * (12 * c - 15) -
                                    3 * std::pow(c, 2) * (6 * std::pow(c, 2) - 15 * c + 10)));
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const std::span<const double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
D::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = std::clamp(input_vector[0], 0.0, 1.0);
    double eta1 = auxiliary_vector[0];
    double eta2 = auxiliary_vector[1];
    double dvol = 0.04;
    double dvap = 0.002;
    double dsurf = 16.0;
    double dgb = 1.6;
    std::vector<double> hessian(1);
    hessian[0] =
        this->prefactor_ * (12 * std::pow(c, 3) * dvol + 6 * std::pow(c, 2) * dvol * (12 * c - 15) +
                            6 * c * dvol * (6 * std::pow(c, 2) - 15 * c + 10) - 2 * dsurf +
                            dvap * (-12 * std::pow(c, 3) - 6 * std::pow(c, 2) * (12 * c - 15) -
                                    6 * c * (6 * std::pow(c, 2) - 15 * c + 10)));
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = A*c**2*(1 - c)**2 + B*(c**2 + (6 - 6*c)*(eta1**2 + eta2**2) - (8 - 4*c)*(eta1**3 +
 * eta2**3) + 3*(eta1**2 + eta2**2)**2)
 */
class Gimp : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  Gimp() : prefactor_(1.0) {}
  explicit Gimp(const double prefactor) : prefactor_(prefactor) {}
  virtual ~Gimp() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const
 * std::span<const double>&, const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Gimp::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = input_vector[1];
    double eta2 = input_vector[2];
    double A = 16.0;
    double B = 1.0;
    double F = A * std::pow(c, 2) * std::pow(1 - c, 2) +
               B * (std::pow(c, 2) + (6 - 6 * c) * (std::pow(eta1, 2) + std::pow(eta2, 2)) -
                    (8 - 4 * c) * (std::pow(eta1, 3) + std::pow(eta2, 3)) +
                    3 * std::pow(std::pow(eta1, 2) + std::pow(eta2, 2), 2));
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gimp::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = input_vector[1];
    double eta2 = input_vector[2];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> gradient(3);
    gradient[0] =
        this->prefactor_ * (A * std::pow(c, 2) * (2 * c - 2) + 2 * A * c * std::pow(1 - c, 2) +
                            B * (2 * c + 4 * std::pow(eta1, 3) - 6 * std::pow(eta1, 2) +
                                 4 * std::pow(eta2, 3) - 6 * std::pow(eta2, 2)));
    gradient[1] =
        this->prefactor_ * (B * (-3 * std::pow(eta1, 2) * (8 - 4 * c) + 2 * eta1 * (6 - 6 * c) +
                                 12 * eta1 * (std::pow(eta1, 2) + std::pow(eta2, 2))));
    gradient[2] =
        this->prefactor_ * (B * (-3 * std::pow(eta2, 2) * (8 - 4 * c) + 2 * eta2 * (6 - 6 * c) +
                                 12 * eta2 * (std::pow(eta1, 2) + std::pow(eta2, 2))));
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const std::span<const double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gimp::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = input_vector[1];
    double eta2 = input_vector[2];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> hessian(9);
    hessian[0] = this->prefactor_ * (2 * A * std::pow(c, 2) + 4 * A * c * (2 * c - 2) +
                                     2 * A * std::pow(1 - c, 2) + 2 * B);
    hessian[1] = this->prefactor_ * (B * (12 * std::pow(eta1, 2) - 12 * eta1));
    hessian[2] = this->prefactor_ * (B * (12 * std::pow(eta2, 2) - 12 * eta2));
    hessian[3] = this->prefactor_ * (B * (12 * std::pow(eta1, 2) - 12 * eta1));
    hessian[4] = this->prefactor_ * (B * (-12 * c + 36 * std::pow(eta1, 2) -
                                          6 * eta1 * (8 - 4 * c) + 12 * std::pow(eta2, 2) + 12));
    hessian[5] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[6] = this->prefactor_ * (B * (12 * std::pow(eta2, 2) - 12 * eta2));
    hessian[7] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[8] = this->prefactor_ * (B * (-12 * c + 12 * std::pow(eta1, 2) +
                                          36 * std::pow(eta2, 2) - 6 * eta2 * (8 - 4 * c) + 12));
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = A*c**2*(1 - c)**2 + B*(c**2 + (6 - 6*c)*(eta1**2 + eta2**2) - (8 - 4*c)*(eta1**3 +
 * eta2**3) + 3*(eta1**2 + eta2**2)**2)
 */
class Gc : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  Gc() : prefactor_(1.0) {}
  explicit Gc(const double prefactor) : prefactor_(prefactor) {}
  virtual ~Gc() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const
 * std::span<const double>&, const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Gc::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = auxiliary_vector[0];
    double eta2 = auxiliary_vector[1];
    double A = 16.0;
    double B = 1.0;
    double F = A * std::pow(c, 2) * std::pow(1 - c, 2) +
               B * (std::pow(c, 2) + (6 - 6 * c) * (std::pow(eta1, 2) + std::pow(eta2, 2)) -
                    (8 - 4 * c) * (std::pow(eta1, 3) + std::pow(eta2, 3)) +
                    3 * std::pow(std::pow(eta1, 2) + std::pow(eta2, 2), 2));
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gc::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = auxiliary_vector[0];
    double eta2 = auxiliary_vector[1];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> gradient(1);
    gradient[0] =
        this->prefactor_ * (A * std::pow(c, 2) * (2 * c - 2) + 2 * A * c * std::pow(1 - c, 2) +
                            B * (2 * c + 4 * std::pow(eta1, 3) - 6 * std::pow(eta1, 2) +
                                 4 * std::pow(eta2, 3) - 6 * std::pow(eta2, 2)));
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const std::span<const double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gc::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double c = input_vector[0];
    double eta1 = auxiliary_vector[0];
    double eta2 = auxiliary_vector[1];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (2 * A * std::pow(c, 2) + 4 * A * c * (2 * c - 2) +
                                     2 * A * std::pow(1 - c, 2) + 2 * B);
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = A*c**2*(1 - c)**2 + B*(c**2 + (6 - 6*c)*(eta1**2 + eta2**2) - (8 - 4*c)*(eta1**3 +
 * eta2**3) + 3*(eta1**2 + eta2**2)**2)
 */
class Geta : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  Geta() : prefactor_(1.0) {}
  explicit Geta(const double prefactor) : prefactor_(prefactor) {}
  virtual ~Geta() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const
 * std::span<const double>&, const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Geta::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double eta1 = input_vector[0];
    double eta2 = input_vector[1];
    double c = auxiliary_vector[0];
    double A = 16.0;
    double B = 1.0;
    double F = A * std::pow(c, 2) * std::pow(1 - c, 2) +
               B * (std::pow(c, 2) + (6 - 6 * c) * (std::pow(eta1, 2) + std::pow(eta2, 2)) -
                    (8 - 4 * c) * (std::pow(eta1, 3) + std::pow(eta2, 3)) +
                    3 * std::pow(std::pow(eta1, 2) + std::pow(eta2, 2), 2));
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Geta::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double eta1 = input_vector[0];
    double eta2 = input_vector[1];
    double c = auxiliary_vector[0];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> gradient(2);
    gradient[0] =
        this->prefactor_ * (B * (-3 * std::pow(eta1, 2) * (8 - 4 * c) + 2 * eta1 * (6 - 6 * c) +
                                 12 * eta1 * (std::pow(eta1, 2) + std::pow(eta2, 2))));
    gradient[1] =
        this->prefactor_ * (B * (-3 * std::pow(eta2, 2) * (8 - 4 * c) + 2 * eta2 * (6 - 6 * c) +
                                 12 * eta2 * (std::pow(eta1, 2) + std::pow(eta2, 2))));
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const std::span<const double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Geta::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double eta1 = input_vector[0];
    double eta2 = input_vector[1];
    double c = auxiliary_vector[0];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> hessian(4);
    hessian[0] = this->prefactor_ * (B * (-12 * c + 36 * std::pow(eta1, 2) -
                                          6 * eta1 * (8 - 4 * c) + 12 * std::pow(eta2, 2) + 12));
    hessian[1] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[2] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[3] = this->prefactor_ * (B * (-12 * c + 12 * std::pow(eta1, 2) +
                                          36 * std::pow(eta2, 2) - 6 * eta2 * (8 - 4 * c) + 12));
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = dot(eta1, eta1) + dot(eta2, eta2) + 5*dot(x, x)/2
 */
class GradC : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  GradC() : prefactor_(1.0) {}
  explicit GradC(const double prefactor) : prefactor_(prefactor) {}
  virtual ~GradC() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const
 * std::span<const double>&, const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
GradC::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> x;
    for (unsigned int i = 0; i < dimension; i++) x.push_back(input_vector[0 * dimension + i]);
    std::vector<double> eta1;
    for (unsigned int i = 0; i < dimension; i++)
      eta1.push_back(auxiliary_vector[0 * dimension + i]);
    std::vector<double> eta2;
    for (unsigned int i = 0; i < dimension; i++)
      eta2.push_back(auxiliary_vector[1 * dimension + i]);
    double F = std::inner_product(eta1.begin(), eta1.end(), eta1.begin(), 0.0) +
               std::inner_product(eta2.begin(), eta2.end(), eta2.begin(), 0.0) +
               2.5 * std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
GradC::GradientF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1, 0.0);
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const std::span<const double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
GradC::HessianF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1, 0.0);
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = dot(eta1, eta1) + dot(eta2, eta2) + 5*dot(x, x)/2
 */
class Gradeta : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  Gradeta() : prefactor_(1.0) {}
  explicit Gradeta(const double prefactor) : prefactor_(prefactor) {}
  virtual ~Gradeta() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const
 * std::span<const double>&, const unsigned int dimension)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Gradeta::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> eta1;
    for (unsigned int i = 0; i < dimension; i++) eta1.push_back(input_vector[0 * dimension + i]);
    std::vector<double> eta2;
    for (unsigned int i = 0; i < dimension; i++) eta2.push_back(input_vector[1 * dimension + i]);
    std::vector<double> x;
    for (unsigned int i = 0; i < dimension; i++) x.push_back(auxiliary_vector[0 * dimension + i]);
    double F = std::inner_product(eta1.begin(), eta1.end(), eta1.begin(), 0.0) +
               std::inner_product(eta2.begin(), eta2.end(), eta2.begin(), 0.0) +
               2.5 * std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gradeta::GradientF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(2, 0.0);
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const std::span<const double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gradeta::HessianF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(4, 0.0);
    return hessian;
  };
  return func;
}
