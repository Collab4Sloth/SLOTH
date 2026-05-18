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

#include "Options/PhysicalPropertiesOptions.hpp"

#include "Coefficients/FunctionCoefficient.hpp"

#pragma once

/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = c**3*dvol*(6*c**2 - 15*c + 10) + c*dsurf*(1 - c) + dgb*(-eta1**2 - eta2**2 - eta3**2 -
 * eta4**2 - eta5**2 - eta6**2 - eta7**2 - eta8**2 - eta9**2 + (eta1 + eta2 + eta3 + eta4 + eta5 +
 * eta6 + eta7 + eta8 + eta9)**2) + dvap*(-c**3*(6*c**2 - 15*c + 10) + 1)
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
    double eta3 = auxiliary_vector[2];
    double eta4 = auxiliary_vector[3];
    double eta5 = auxiliary_vector[4];
    double eta6 = auxiliary_vector[5];
    double eta7 = auxiliary_vector[6];
    double eta8 = auxiliary_vector[7];
    double eta9 = auxiliary_vector[8];
    double dvol = 0.04;
    double dvap = 0.002;
    double dsurf = 16.0;
    double dgb = 1.6;
    double F = std::pow(c, 3) * dvol * (6 * std::pow(c, 2) - 15 * c + 10) + c * dsurf * (1 - c) +
               dgb * (-std::pow(eta1, 2) - std::pow(eta2, 2) - std::pow(eta3, 2) -
                      std::pow(eta4, 2) - std::pow(eta5, 2) - std::pow(eta6, 2) -
                      std::pow(eta7, 2) - std::pow(eta8, 2) - std::pow(eta9, 2) +
                      std::pow(eta1 + eta2 + eta3 + eta4 + eta5 + eta6 + eta7 + eta8 + eta9, 2)) +
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
    double eta3 = auxiliary_vector[2];
    double eta4 = auxiliary_vector[3];
    double eta5 = auxiliary_vector[4];
    double eta6 = auxiliary_vector[5];
    double eta7 = auxiliary_vector[6];
    double eta8 = auxiliary_vector[7];
    double eta9 = auxiliary_vector[8];
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
    double eta3 = auxiliary_vector[2];
    double eta4 = auxiliary_vector[3];
    double eta5 = auxiliary_vector[4];
    double eta6 = auxiliary_vector[5];
    double eta7 = auxiliary_vector[6];
    double eta8 = auxiliary_vector[7];
    double eta9 = auxiliary_vector[8];
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
 *       F = A*c**2*(1 - c)**2 + B*(c**2 + (6 - 6*c)*(eta1**2 + eta2**2 + eta3**2 + eta4**2 +
 * eta5**2 + eta6**2 + eta7**2 + eta8**2 + eta9**2) - (8 - 4*c)*(eta1**3 + eta2**3 + eta3**3 +
 * eta4**3 + eta5**3 + eta6**3 + eta7**3 + eta8**3 + eta9**3) + 3*(eta1**2 + eta2**2 + eta3**2 +
 * eta4**2 + eta5**2 + eta6**2 + eta7**2 + eta8**2 + eta9**2)**2)
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
    double eta3 = input_vector[3];
    double eta4 = input_vector[4];
    double eta5 = input_vector[5];
    double eta6 = input_vector[6];
    double eta7 = input_vector[7];
    double eta8 = input_vector[8];
    double eta9 = input_vector[9];
    double A = 16.0;
    double B = 1.0;
    double F = A * std::pow(c, 2) * std::pow(1 - c, 2) +
               B * (std::pow(c, 2) +
                    (6 - 6 * c) * (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                   std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                   std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2)) -
                    (8 - 4 * c) * (std::pow(eta1, 3) + std::pow(eta2, 3) + std::pow(eta3, 3) +
                                   std::pow(eta4, 3) + std::pow(eta5, 3) + std::pow(eta6, 3) +
                                   std::pow(eta7, 3) + std::pow(eta8, 3) + std::pow(eta9, 3)) +
                    3 * std::pow(std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                     std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                     std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2),
                                 2));
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
    double eta3 = input_vector[3];
    double eta4 = input_vector[4];
    double eta5 = input_vector[5];
    double eta6 = input_vector[6];
    double eta7 = input_vector[7];
    double eta8 = input_vector[8];
    double eta9 = input_vector[9];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> gradient(10);
    gradient[0] =
        this->prefactor_ *
        (A * std::pow(c, 2) * (2 * c - 2) + 2 * A * c * std::pow(1 - c, 2) +
         B * (2 * c + 4 * std::pow(eta1, 3) - 6 * std::pow(eta1, 2) + 4 * std::pow(eta2, 3) -
              6 * std::pow(eta2, 2) + 4 * std::pow(eta3, 3) - 6 * std::pow(eta3, 2) +
              4 * std::pow(eta4, 3) - 6 * std::pow(eta4, 2) + 4 * std::pow(eta5, 3) -
              6 * std::pow(eta5, 2) + 4 * std::pow(eta6, 3) - 6 * std::pow(eta6, 2) +
              4 * std::pow(eta7, 3) - 6 * std::pow(eta7, 2) + 4 * std::pow(eta8, 3) -
              6 * std::pow(eta8, 2) + 4 * std::pow(eta9, 3) - 6 * std::pow(eta9, 2)));
    gradient[1] =
        this->prefactor_ * (B * (-3 * std::pow(eta1, 2) * (8 - 4 * c) + 2 * eta1 * (6 - 6 * c) +
                                 12 * eta1 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[2] =
        this->prefactor_ * (B * (-3 * std::pow(eta2, 2) * (8 - 4 * c) + 2 * eta2 * (6 - 6 * c) +
                                 12 * eta2 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[3] =
        this->prefactor_ * (B * (-3 * std::pow(eta3, 2) * (8 - 4 * c) + 2 * eta3 * (6 - 6 * c) +
                                 12 * eta3 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[4] =
        this->prefactor_ * (B * (-3 * std::pow(eta4, 2) * (8 - 4 * c) + 2 * eta4 * (6 - 6 * c) +
                                 12 * eta4 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[5] =
        this->prefactor_ * (B * (-3 * std::pow(eta5, 2) * (8 - 4 * c) + 2 * eta5 * (6 - 6 * c) +
                                 12 * eta5 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[6] =
        this->prefactor_ * (B * (-3 * std::pow(eta6, 2) * (8 - 4 * c) + 2 * eta6 * (6 - 6 * c) +
                                 12 * eta6 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[7] =
        this->prefactor_ * (B * (-3 * std::pow(eta7, 2) * (8 - 4 * c) + 2 * eta7 * (6 - 6 * c) +
                                 12 * eta7 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[8] =
        this->prefactor_ * (B * (-3 * std::pow(eta8, 2) * (8 - 4 * c) + 2 * eta8 * (6 - 6 * c) +
                                 12 * eta8 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[9] =
        this->prefactor_ * (B * (-3 * std::pow(eta9, 2) * (8 - 4 * c) + 2 * eta9 * (6 - 6 * c) +
                                 12 * eta9 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
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
    double eta3 = input_vector[3];
    double eta4 = input_vector[4];
    double eta5 = input_vector[5];
    double eta6 = input_vector[6];
    double eta7 = input_vector[7];
    double eta8 = input_vector[8];
    double eta9 = input_vector[9];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> hessian(100);
    hessian[0] = this->prefactor_ * (2 * A * std::pow(c, 2) + 4 * A * c * (2 * c - 2) +
                                     2 * A * std::pow(1 - c, 2) + 2 * B);
    hessian[1] = this->prefactor_ * (B * (12 * std::pow(eta1, 2) - 12 * eta1));
    hessian[2] = this->prefactor_ * (B * (12 * std::pow(eta2, 2) - 12 * eta2));
    hessian[3] = this->prefactor_ * (B * (12 * std::pow(eta3, 2) - 12 * eta3));
    hessian[4] = this->prefactor_ * (B * (12 * std::pow(eta4, 2) - 12 * eta4));
    hessian[5] = this->prefactor_ * (B * (12 * std::pow(eta5, 2) - 12 * eta5));
    hessian[6] = this->prefactor_ * (B * (12 * std::pow(eta6, 2) - 12 * eta6));
    hessian[7] = this->prefactor_ * (B * (12 * std::pow(eta7, 2) - 12 * eta7));
    hessian[8] = this->prefactor_ * (B * (12 * std::pow(eta8, 2) - 12 * eta8));
    hessian[9] = this->prefactor_ * (B * (12 * std::pow(eta9, 2) - 12 * eta9));
    hessian[10] = this->prefactor_ * (B * (12 * std::pow(eta1, 2) - 12 * eta1));
    hessian[11] = this->prefactor_ *
                  (B * (-12 * c + 36 * std::pow(eta1, 2) - 6 * eta1 * (8 - 4 * c) +
                        12 * std::pow(eta2, 2) + 12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[12] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[13] = this->prefactor_ * (24 * B * eta1 * eta3);
    hessian[14] = this->prefactor_ * (24 * B * eta1 * eta4);
    hessian[15] = this->prefactor_ * (24 * B * eta1 * eta5);
    hessian[16] = this->prefactor_ * (24 * B * eta1 * eta6);
    hessian[17] = this->prefactor_ * (24 * B * eta1 * eta7);
    hessian[18] = this->prefactor_ * (24 * B * eta1 * eta8);
    hessian[19] = this->prefactor_ * (24 * B * eta1 * eta9);
    hessian[20] = this->prefactor_ * (B * (12 * std::pow(eta2, 2) - 12 * eta2));
    hessian[21] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[22] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 36 * std::pow(eta2, 2) -
                        6 * eta2 * (8 - 4 * c) + 12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[23] = this->prefactor_ * (24 * B * eta2 * eta3);
    hessian[24] = this->prefactor_ * (24 * B * eta2 * eta4);
    hessian[25] = this->prefactor_ * (24 * B * eta2 * eta5);
    hessian[26] = this->prefactor_ * (24 * B * eta2 * eta6);
    hessian[27] = this->prefactor_ * (24 * B * eta2 * eta7);
    hessian[28] = this->prefactor_ * (24 * B * eta2 * eta8);
    hessian[29] = this->prefactor_ * (24 * B * eta2 * eta9);
    hessian[30] = this->prefactor_ * (B * (12 * std::pow(eta3, 2) - 12 * eta3));
    hessian[31] = this->prefactor_ * (24 * B * eta1 * eta3);
    hessian[32] = this->prefactor_ * (24 * B * eta2 * eta3);
    hessian[33] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        36 * std::pow(eta3, 2) - 6 * eta3 * (8 - 4 * c) + 12 * std::pow(eta4, 2) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[34] = this->prefactor_ * (24 * B * eta3 * eta4);
    hessian[35] = this->prefactor_ * (24 * B * eta3 * eta5);
    hessian[36] = this->prefactor_ * (24 * B * eta3 * eta6);
    hessian[37] = this->prefactor_ * (24 * B * eta3 * eta7);
    hessian[38] = this->prefactor_ * (24 * B * eta3 * eta8);
    hessian[39] = this->prefactor_ * (24 * B * eta3 * eta9);
    hessian[40] = this->prefactor_ * (B * (12 * std::pow(eta4, 2) - 12 * eta4));
    hessian[41] = this->prefactor_ * (24 * B * eta1 * eta4);
    hessian[42] = this->prefactor_ * (24 * B * eta2 * eta4);
    hessian[43] = this->prefactor_ * (24 * B * eta3 * eta4);
    hessian[44] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 36 * std::pow(eta4, 2) - 6 * eta4 * (8 - 4 * c) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[45] = this->prefactor_ * (24 * B * eta4 * eta5);
    hessian[46] = this->prefactor_ * (24 * B * eta4 * eta6);
    hessian[47] = this->prefactor_ * (24 * B * eta4 * eta7);
    hessian[48] = this->prefactor_ * (24 * B * eta4 * eta8);
    hessian[49] = this->prefactor_ * (24 * B * eta4 * eta9);
    hessian[50] = this->prefactor_ * (B * (12 * std::pow(eta5, 2) - 12 * eta5));
    hessian[51] = this->prefactor_ * (24 * B * eta1 * eta5);
    hessian[52] = this->prefactor_ * (24 * B * eta2 * eta5);
    hessian[53] = this->prefactor_ * (24 * B * eta3 * eta5);
    hessian[54] = this->prefactor_ * (24 * B * eta4 * eta5);
    hessian[55] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 36 * std::pow(eta5, 2) -
                        6 * eta5 * (8 - 4 * c) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[56] = this->prefactor_ * (24 * B * eta5 * eta6);
    hessian[57] = this->prefactor_ * (24 * B * eta5 * eta7);
    hessian[58] = this->prefactor_ * (24 * B * eta5 * eta8);
    hessian[59] = this->prefactor_ * (24 * B * eta5 * eta9);
    hessian[60] = this->prefactor_ * (B * (12 * std::pow(eta6, 2) - 12 * eta6));
    hessian[61] = this->prefactor_ * (24 * B * eta1 * eta6);
    hessian[62] = this->prefactor_ * (24 * B * eta2 * eta6);
    hessian[63] = this->prefactor_ * (24 * B * eta3 * eta6);
    hessian[64] = this->prefactor_ * (24 * B * eta4 * eta6);
    hessian[65] = this->prefactor_ * (24 * B * eta5 * eta6);
    hessian[66] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        36 * std::pow(eta6, 2) - 6 * eta6 * (8 - 4 * c) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[67] = this->prefactor_ * (24 * B * eta6 * eta7);
    hessian[68] = this->prefactor_ * (24 * B * eta6 * eta8);
    hessian[69] = this->prefactor_ * (24 * B * eta6 * eta9);
    hessian[70] = this->prefactor_ * (B * (12 * std::pow(eta7, 2) - 12 * eta7));
    hessian[71] = this->prefactor_ * (24 * B * eta1 * eta7);
    hessian[72] = this->prefactor_ * (24 * B * eta2 * eta7);
    hessian[73] = this->prefactor_ * (24 * B * eta3 * eta7);
    hessian[74] = this->prefactor_ * (24 * B * eta4 * eta7);
    hessian[75] = this->prefactor_ * (24 * B * eta5 * eta7);
    hessian[76] = this->prefactor_ * (24 * B * eta6 * eta7);
    hessian[77] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        12 * std::pow(eta6, 2) + 36 * std::pow(eta7, 2) - 6 * eta7 * (8 - 4 * c) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[78] = this->prefactor_ * (24 * B * eta7 * eta8);
    hessian[79] = this->prefactor_ * (24 * B * eta7 * eta9);
    hessian[80] = this->prefactor_ * (B * (12 * std::pow(eta8, 2) - 12 * eta8));
    hessian[81] = this->prefactor_ * (24 * B * eta1 * eta8);
    hessian[82] = this->prefactor_ * (24 * B * eta2 * eta8);
    hessian[83] = this->prefactor_ * (24 * B * eta3 * eta8);
    hessian[84] = this->prefactor_ * (24 * B * eta4 * eta8);
    hessian[85] = this->prefactor_ * (24 * B * eta5 * eta8);
    hessian[86] = this->prefactor_ * (24 * B * eta6 * eta8);
    hessian[87] = this->prefactor_ * (24 * B * eta7 * eta8);
    hessian[88] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) + 36 * std::pow(eta8, 2) -
                        6 * eta8 * (8 - 4 * c) + 12 * std::pow(eta9, 2) + 12));
    hessian[89] = this->prefactor_ * (24 * B * eta8 * eta9);
    hessian[90] = this->prefactor_ * (B * (12 * std::pow(eta9, 2) - 12 * eta9));
    hessian[91] = this->prefactor_ * (24 * B * eta1 * eta9);
    hessian[92] = this->prefactor_ * (24 * B * eta2 * eta9);
    hessian[93] = this->prefactor_ * (24 * B * eta3 * eta9);
    hessian[94] = this->prefactor_ * (24 * B * eta4 * eta9);
    hessian[95] = this->prefactor_ * (24 * B * eta5 * eta9);
    hessian[96] = this->prefactor_ * (24 * B * eta6 * eta9);
    hessian[97] = this->prefactor_ * (24 * B * eta7 * eta9);
    hessian[98] = this->prefactor_ * (24 * B * eta8 * eta9);
    hessian[99] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) + 12 * std::pow(eta8, 2) +
                        36 * std::pow(eta9, 2) - 6 * eta9 * (8 - 4 * c) + 12));
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = A*c**2*(1 - c)**2 + B*(c**2 + (6 - 6*c)*(eta1**2 + eta2**2 + eta3**2 + eta4**2 +
 * eta5**2 + eta6**2 + eta7**2 + eta8**2 + eta9**2) - (8 - 4*c)*(eta1**3 + eta2**3 + eta3**3 +
 * eta4**3 + eta5**3 + eta6**3 + eta7**3 + eta8**3 + eta9**3) + 3*(eta1**2 + eta2**2 + eta3**2 +
 * eta4**2 + eta5**2 + eta6**2 + eta7**2 + eta8**2 + eta9**2)**2)
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
    double eta3 = auxiliary_vector[2];
    double eta4 = auxiliary_vector[3];
    double eta5 = auxiliary_vector[4];
    double eta6 = auxiliary_vector[5];
    double eta7 = auxiliary_vector[6];
    double eta8 = auxiliary_vector[7];
    double eta9 = auxiliary_vector[8];
    double A = 16.0;
    double B = 1.0;
    double F = A * std::pow(c, 2) * std::pow(1 - c, 2) +
               B * (std::pow(c, 2) +
                    (6 - 6 * c) * (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                   std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                   std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2)) -
                    (8 - 4 * c) * (std::pow(eta1, 3) + std::pow(eta2, 3) + std::pow(eta3, 3) +
                                   std::pow(eta4, 3) + std::pow(eta5, 3) + std::pow(eta6, 3) +
                                   std::pow(eta7, 3) + std::pow(eta8, 3) + std::pow(eta9, 3)) +
                    3 * std::pow(std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                     std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                     std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2),
                                 2));
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
    double eta3 = auxiliary_vector[2];
    double eta4 = auxiliary_vector[3];
    double eta5 = auxiliary_vector[4];
    double eta6 = auxiliary_vector[5];
    double eta7 = auxiliary_vector[6];
    double eta8 = auxiliary_vector[7];
    double eta9 = auxiliary_vector[8];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> gradient(1);
    gradient[0] =
        this->prefactor_ *
        (A * std::pow(c, 2) * (2 * c - 2) + 2 * A * c * std::pow(1 - c, 2) +
         B * (2 * c + 4 * std::pow(eta1, 3) - 6 * std::pow(eta1, 2) + 4 * std::pow(eta2, 3) -
              6 * std::pow(eta2, 2) + 4 * std::pow(eta3, 3) - 6 * std::pow(eta3, 2) +
              4 * std::pow(eta4, 3) - 6 * std::pow(eta4, 2) + 4 * std::pow(eta5, 3) -
              6 * std::pow(eta5, 2) + 4 * std::pow(eta6, 3) - 6 * std::pow(eta6, 2) +
              4 * std::pow(eta7, 3) - 6 * std::pow(eta7, 2) + 4 * std::pow(eta8, 3) -
              6 * std::pow(eta8, 2) + 4 * std::pow(eta9, 3) - 6 * std::pow(eta9, 2)));
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
    double eta3 = auxiliary_vector[2];
    double eta4 = auxiliary_vector[3];
    double eta5 = auxiliary_vector[4];
    double eta6 = auxiliary_vector[5];
    double eta7 = auxiliary_vector[6];
    double eta8 = auxiliary_vector[7];
    double eta9 = auxiliary_vector[8];
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
 *       F = A*c**2*(1 - c)**2 + B*(c**2 + (6 - 6*c)*(eta1**2 + eta2**2 + eta3**2 + eta4**2 +
 * eta5**2 + eta6**2 + eta7**2 + eta8**2 + eta9**2) - (8 - 4*c)*(eta1**3 + eta2**3 + eta3**3 +
 * eta4**3 + eta5**3 + eta6**3 + eta7**3 + eta8**3 + eta9**3) + 3*(eta1**2 + eta2**2 + eta3**2 +
 * eta4**2 + eta5**2 + eta6**2 + eta7**2 + eta8**2 + eta9**2)**2)
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
    double eta3 = input_vector[2];
    double eta4 = input_vector[3];
    double eta5 = input_vector[4];
    double eta6 = input_vector[5];
    double eta7 = input_vector[6];
    double eta8 = input_vector[7];
    double eta9 = input_vector[8];
    double c = auxiliary_vector[0];
    double A = 16.0;
    double B = 1.0;
    double F = A * std::pow(c, 2) * std::pow(1 - c, 2) +
               B * (std::pow(c, 2) +
                    (6 - 6 * c) * (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                   std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                   std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2)) -
                    (8 - 4 * c) * (std::pow(eta1, 3) + std::pow(eta2, 3) + std::pow(eta3, 3) +
                                   std::pow(eta4, 3) + std::pow(eta5, 3) + std::pow(eta6, 3) +
                                   std::pow(eta7, 3) + std::pow(eta8, 3) + std::pow(eta9, 3)) +
                    3 * std::pow(std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                     std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                     std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2),
                                 2));
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
    double eta3 = input_vector[2];
    double eta4 = input_vector[3];
    double eta5 = input_vector[4];
    double eta6 = input_vector[5];
    double eta7 = input_vector[6];
    double eta8 = input_vector[7];
    double eta9 = input_vector[8];
    double c = auxiliary_vector[0];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> gradient(9);
    gradient[0] =
        this->prefactor_ * (B * (-3 * std::pow(eta1, 2) * (8 - 4 * c) + 2 * eta1 * (6 - 6 * c) +
                                 12 * eta1 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[1] =
        this->prefactor_ * (B * (-3 * std::pow(eta2, 2) * (8 - 4 * c) + 2 * eta2 * (6 - 6 * c) +
                                 12 * eta2 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[2] =
        this->prefactor_ * (B * (-3 * std::pow(eta3, 2) * (8 - 4 * c) + 2 * eta3 * (6 - 6 * c) +
                                 12 * eta3 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[3] =
        this->prefactor_ * (B * (-3 * std::pow(eta4, 2) * (8 - 4 * c) + 2 * eta4 * (6 - 6 * c) +
                                 12 * eta4 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[4] =
        this->prefactor_ * (B * (-3 * std::pow(eta5, 2) * (8 - 4 * c) + 2 * eta5 * (6 - 6 * c) +
                                 12 * eta5 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[5] =
        this->prefactor_ * (B * (-3 * std::pow(eta6, 2) * (8 - 4 * c) + 2 * eta6 * (6 - 6 * c) +
                                 12 * eta6 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[6] =
        this->prefactor_ * (B * (-3 * std::pow(eta7, 2) * (8 - 4 * c) + 2 * eta7 * (6 - 6 * c) +
                                 12 * eta7 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[7] =
        this->prefactor_ * (B * (-3 * std::pow(eta8, 2) * (8 - 4 * c) + 2 * eta8 * (6 - 6 * c) +
                                 12 * eta8 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
    gradient[8] =
        this->prefactor_ * (B * (-3 * std::pow(eta9, 2) * (8 - 4 * c) + 2 * eta9 * (6 - 6 * c) +
                                 12 * eta9 *
                                     (std::pow(eta1, 2) + std::pow(eta2, 2) + std::pow(eta3, 2) +
                                      std::pow(eta4, 2) + std::pow(eta5, 2) + std::pow(eta6, 2) +
                                      std::pow(eta7, 2) + std::pow(eta8, 2) + std::pow(eta9, 2))));
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
    double eta3 = input_vector[2];
    double eta4 = input_vector[3];
    double eta5 = input_vector[4];
    double eta6 = input_vector[5];
    double eta7 = input_vector[6];
    double eta8 = input_vector[7];
    double eta9 = input_vector[8];
    double c = auxiliary_vector[0];
    double A = 16.0;
    double B = 1.0;
    std::vector<double> hessian(81);
    hessian[0] = this->prefactor_ *
                 (B * (-12 * c + 36 * std::pow(eta1, 2) - 6 * eta1 * (8 - 4 * c) +
                       12 * std::pow(eta2, 2) + 12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) +
                       12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                       12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[1] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[2] = this->prefactor_ * (24 * B * eta1 * eta3);
    hessian[3] = this->prefactor_ * (24 * B * eta1 * eta4);
    hessian[4] = this->prefactor_ * (24 * B * eta1 * eta5);
    hessian[5] = this->prefactor_ * (24 * B * eta1 * eta6);
    hessian[6] = this->prefactor_ * (24 * B * eta1 * eta7);
    hessian[7] = this->prefactor_ * (24 * B * eta1 * eta8);
    hessian[8] = this->prefactor_ * (24 * B * eta1 * eta9);
    hessian[9] = this->prefactor_ * (24 * B * eta1 * eta2);
    hessian[10] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 36 * std::pow(eta2, 2) -
                        6 * eta2 * (8 - 4 * c) + 12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[11] = this->prefactor_ * (24 * B * eta2 * eta3);
    hessian[12] = this->prefactor_ * (24 * B * eta2 * eta4);
    hessian[13] = this->prefactor_ * (24 * B * eta2 * eta5);
    hessian[14] = this->prefactor_ * (24 * B * eta2 * eta6);
    hessian[15] = this->prefactor_ * (24 * B * eta2 * eta7);
    hessian[16] = this->prefactor_ * (24 * B * eta2 * eta8);
    hessian[17] = this->prefactor_ * (24 * B * eta2 * eta9);
    hessian[18] = this->prefactor_ * (24 * B * eta1 * eta3);
    hessian[19] = this->prefactor_ * (24 * B * eta2 * eta3);
    hessian[20] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        36 * std::pow(eta3, 2) - 6 * eta3 * (8 - 4 * c) + 12 * std::pow(eta4, 2) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[21] = this->prefactor_ * (24 * B * eta3 * eta4);
    hessian[22] = this->prefactor_ * (24 * B * eta3 * eta5);
    hessian[23] = this->prefactor_ * (24 * B * eta3 * eta6);
    hessian[24] = this->prefactor_ * (24 * B * eta3 * eta7);
    hessian[25] = this->prefactor_ * (24 * B * eta3 * eta8);
    hessian[26] = this->prefactor_ * (24 * B * eta3 * eta9);
    hessian[27] = this->prefactor_ * (24 * B * eta1 * eta4);
    hessian[28] = this->prefactor_ * (24 * B * eta2 * eta4);
    hessian[29] = this->prefactor_ * (24 * B * eta3 * eta4);
    hessian[30] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 36 * std::pow(eta4, 2) - 6 * eta4 * (8 - 4 * c) +
                        12 * std::pow(eta5, 2) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[31] = this->prefactor_ * (24 * B * eta4 * eta5);
    hessian[32] = this->prefactor_ * (24 * B * eta4 * eta6);
    hessian[33] = this->prefactor_ * (24 * B * eta4 * eta7);
    hessian[34] = this->prefactor_ * (24 * B * eta4 * eta8);
    hessian[35] = this->prefactor_ * (24 * B * eta4 * eta9);
    hessian[36] = this->prefactor_ * (24 * B * eta1 * eta5);
    hessian[37] = this->prefactor_ * (24 * B * eta2 * eta5);
    hessian[38] = this->prefactor_ * (24 * B * eta3 * eta5);
    hessian[39] = this->prefactor_ * (24 * B * eta4 * eta5);
    hessian[40] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 36 * std::pow(eta5, 2) -
                        6 * eta5 * (8 - 4 * c) + 12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[41] = this->prefactor_ * (24 * B * eta5 * eta6);
    hessian[42] = this->prefactor_ * (24 * B * eta5 * eta7);
    hessian[43] = this->prefactor_ * (24 * B * eta5 * eta8);
    hessian[44] = this->prefactor_ * (24 * B * eta5 * eta9);
    hessian[45] = this->prefactor_ * (24 * B * eta1 * eta6);
    hessian[46] = this->prefactor_ * (24 * B * eta2 * eta6);
    hessian[47] = this->prefactor_ * (24 * B * eta3 * eta6);
    hessian[48] = this->prefactor_ * (24 * B * eta4 * eta6);
    hessian[49] = this->prefactor_ * (24 * B * eta5 * eta6);
    hessian[50] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        36 * std::pow(eta6, 2) - 6 * eta6 * (8 - 4 * c) + 12 * std::pow(eta7, 2) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[51] = this->prefactor_ * (24 * B * eta6 * eta7);
    hessian[52] = this->prefactor_ * (24 * B * eta6 * eta8);
    hessian[53] = this->prefactor_ * (24 * B * eta6 * eta9);
    hessian[54] = this->prefactor_ * (24 * B * eta1 * eta7);
    hessian[55] = this->prefactor_ * (24 * B * eta2 * eta7);
    hessian[56] = this->prefactor_ * (24 * B * eta3 * eta7);
    hessian[57] = this->prefactor_ * (24 * B * eta4 * eta7);
    hessian[58] = this->prefactor_ * (24 * B * eta5 * eta7);
    hessian[59] = this->prefactor_ * (24 * B * eta6 * eta7);
    hessian[60] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        12 * std::pow(eta6, 2) + 36 * std::pow(eta7, 2) - 6 * eta7 * (8 - 4 * c) +
                        12 * std::pow(eta8, 2) + 12 * std::pow(eta9, 2) + 12));
    hessian[61] = this->prefactor_ * (24 * B * eta7 * eta8);
    hessian[62] = this->prefactor_ * (24 * B * eta7 * eta9);
    hessian[63] = this->prefactor_ * (24 * B * eta1 * eta8);
    hessian[64] = this->prefactor_ * (24 * B * eta2 * eta8);
    hessian[65] = this->prefactor_ * (24 * B * eta3 * eta8);
    hessian[66] = this->prefactor_ * (24 * B * eta4 * eta8);
    hessian[67] = this->prefactor_ * (24 * B * eta5 * eta8);
    hessian[68] = this->prefactor_ * (24 * B * eta6 * eta8);
    hessian[69] = this->prefactor_ * (24 * B * eta7 * eta8);
    hessian[70] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) + 36 * std::pow(eta8, 2) -
                        6 * eta8 * (8 - 4 * c) + 12 * std::pow(eta9, 2) + 12));
    hessian[71] = this->prefactor_ * (24 * B * eta8 * eta9);
    hessian[72] = this->prefactor_ * (24 * B * eta1 * eta9);
    hessian[73] = this->prefactor_ * (24 * B * eta2 * eta9);
    hessian[74] = this->prefactor_ * (24 * B * eta3 * eta9);
    hessian[75] = this->prefactor_ * (24 * B * eta4 * eta9);
    hessian[76] = this->prefactor_ * (24 * B * eta5 * eta9);
    hessian[77] = this->prefactor_ * (24 * B * eta6 * eta9);
    hessian[78] = this->prefactor_ * (24 * B * eta7 * eta9);
    hessian[79] = this->prefactor_ * (24 * B * eta8 * eta9);
    hessian[80] = this->prefactor_ *
                  (B * (-12 * c + 12 * std::pow(eta1, 2) + 12 * std::pow(eta2, 2) +
                        12 * std::pow(eta3, 2) + 12 * std::pow(eta4, 2) + 12 * std::pow(eta5, 2) +
                        12 * std::pow(eta6, 2) + 12 * std::pow(eta7, 2) + 12 * std::pow(eta8, 2) +
                        36 * std::pow(eta9, 2) - 6 * eta9 * (8 - 4 * c) + 12));
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = dot(eta1, eta1) + dot(eta2, eta2) + dot(eta3, eta3) + dot(eta4, eta4) + dot(eta5, eta5)
 * + dot(eta6, eta6) + dot(eta7, eta7) + dot(eta8, eta8) + dot(eta9, eta9) + 5*dot(x, x)/2
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
    std::vector<double> eta3;
    for (unsigned int i = 0; i < dimension; i++)
      eta3.push_back(auxiliary_vector[2 * dimension + i]);
    std::vector<double> eta4;
    for (unsigned int i = 0; i < dimension; i++)
      eta4.push_back(auxiliary_vector[3 * dimension + i]);
    std::vector<double> eta5;
    for (unsigned int i = 0; i < dimension; i++)
      eta5.push_back(auxiliary_vector[4 * dimension + i]);
    std::vector<double> eta6;
    for (unsigned int i = 0; i < dimension; i++)
      eta6.push_back(auxiliary_vector[5 * dimension + i]);
    std::vector<double> eta7;
    for (unsigned int i = 0; i < dimension; i++)
      eta7.push_back(auxiliary_vector[6 * dimension + i]);
    std::vector<double> eta8;
    for (unsigned int i = 0; i < dimension; i++)
      eta8.push_back(auxiliary_vector[7 * dimension + i]);
    std::vector<double> eta9;
    for (unsigned int i = 0; i < dimension; i++)
      eta9.push_back(auxiliary_vector[8 * dimension + i]);
    double F = std::inner_product(eta1.begin(), eta1.end(), eta1.begin(), 0.0) +
               std::inner_product(eta2.begin(), eta2.end(), eta2.begin(), 0.0) +
               std::inner_product(eta3.begin(), eta3.end(), eta3.begin(), 0.0) +
               std::inner_product(eta4.begin(), eta4.end(), eta4.begin(), 0.0) +
               std::inner_product(eta5.begin(), eta5.end(), eta5.begin(), 0.0) +
               std::inner_product(eta6.begin(), eta6.end(), eta6.begin(), 0.0) +
               std::inner_product(eta7.begin(), eta7.end(), eta7.begin(), 0.0) +
               std::inner_product(eta8.begin(), eta8.end(), eta8.begin(), 0.0) +
               std::inner_product(eta9.begin(), eta9.end(), eta9.begin(), 0.0) +
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
 *       F = dot(eta1, eta1) + dot(eta2, eta2) + dot(eta3, eta3) + dot(eta4, eta4) + dot(eta5, eta5)
 * + dot(eta6, eta6) + dot(eta7, eta7) + dot(eta8, eta8) + dot(eta9, eta9) + 5*dot(x, x)/2
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
    std::vector<double> eta3;
    for (unsigned int i = 0; i < dimension; i++) eta3.push_back(input_vector[2 * dimension + i]);
    std::vector<double> eta4;
    for (unsigned int i = 0; i < dimension; i++) eta4.push_back(input_vector[3 * dimension + i]);
    std::vector<double> eta5;
    for (unsigned int i = 0; i < dimension; i++) eta5.push_back(input_vector[4 * dimension + i]);
    std::vector<double> eta6;
    for (unsigned int i = 0; i < dimension; i++) eta6.push_back(input_vector[5 * dimension + i]);
    std::vector<double> eta7;
    for (unsigned int i = 0; i < dimension; i++) eta7.push_back(input_vector[6 * dimension + i]);
    std::vector<double> eta8;
    for (unsigned int i = 0; i < dimension; i++) eta8.push_back(input_vector[7 * dimension + i]);
    std::vector<double> eta9;
    for (unsigned int i = 0; i < dimension; i++) eta9.push_back(input_vector[8 * dimension + i]);
    std::vector<double> x;
    for (unsigned int i = 0; i < dimension; i++) x.push_back(auxiliary_vector[0 * dimension + i]);
    double F = std::inner_product(eta1.begin(), eta1.end(), eta1.begin(), 0.0) +
               std::inner_product(eta2.begin(), eta2.end(), eta2.begin(), 0.0) +
               std::inner_product(eta3.begin(), eta3.end(), eta3.begin(), 0.0) +
               std::inner_product(eta4.begin(), eta4.end(), eta4.begin(), 0.0) +
               std::inner_product(eta5.begin(), eta5.end(), eta5.begin(), 0.0) +
               std::inner_product(eta6.begin(), eta6.end(), eta6.begin(), 0.0) +
               std::inner_product(eta7.begin(), eta7.end(), eta7.begin(), 0.0) +
               std::inner_product(eta8.begin(), eta8.end(), eta8.begin(), 0.0) +
               std::inner_product(eta9.begin(), eta9.end(), eta9.begin(), 0.0) +
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
    std::vector<double> gradient(9, 0.0);
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
    std::vector<double> hessian(81, 0.0);
    return hessian;
  };
  return func;
}
