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
#include <stdexcept>
#include <vector>

#include "Options/PhysicalPropertiesOptions.hpp"

#include "Coefficients/FunctionCoefficient.hpp"

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Evaluation of D at one point, on the host or on the device
 *
 * Same expressions as D (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 *
 * Bounded variables (F, gradient and Hessian evaluated at the clamped values):
 *   c in [eps, (1.0 - eps)]
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct DEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 9;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;
  double dvol = 0.04;
  double dvap = 0.002;
  double dsurf = 16.0;
  double dgb = 1.6;
  double eps = 1.e-10;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    // Bounded variables: clamped to their bounds
    [[maybe_unused]] const double c = fmin((1.0 - eps), fmax(eps, values[0]));
    [[maybe_unused]] const double eta1 = auxiliary_values[0];
    [[maybe_unused]] const double eta2 = auxiliary_values[1];
    [[maybe_unused]] const double eta3 = auxiliary_values[2];
    [[maybe_unused]] const double eta4 = auxiliary_values[3];
    [[maybe_unused]] const double eta5 = auxiliary_values[4];
    [[maybe_unused]] const double eta6 = auxiliary_values[5];
    [[maybe_unused]] const double eta7 = auxiliary_values[6];
    [[maybe_unused]] const double eta8 = auxiliary_values[7];
    [[maybe_unused]] const double eta9 = auxiliary_values[8];
    if (F != nullptr) {
      *F = prefactor_ *
           (pow(c, 3) * dvol * (6 * pow(c, 2) - 15 * c + 10) + c * dsurf * (1 - c) +
            dgb * (-pow(eta1, 2) - pow(eta2, 2) - pow(eta3, 2) - pow(eta4, 2) - pow(eta5, 2) -
                   pow(eta6, 2) - pow(eta7, 2) - pow(eta8, 2) - pow(eta9, 2) +
                   pow(eta1 + eta2 + eta3 + eta4 + eta5 + eta6 + eta7 + eta8 + eta9, 2)) +
            dvap * (-pow(c, 3) * (6 * pow(c, 2) - 15 * c + 10) + 1));
    }
    if (gradient != nullptr) {
      gradient[0] =
          prefactor_ *
          (pow(c, 3) * dvol * (12 * c - 15) + 3 * pow(c, 2) * dvol * (6 * pow(c, 2) - 15 * c + 10) -
           c * dsurf + dsurf * (1 - c) +
           dvap * (-pow(c, 3) * (12 * c - 15) - 3 * pow(c, 2) * (6 * pow(c, 2) - 15 * c + 10)));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ * (12 * pow(c, 3) * dvol + 6 * pow(c, 2) * dvol * (12 * c - 15) +
                                 6 * c * dvol * (6 * pow(c, 2) - 15 * c + 10) - 2 * dsurf +
                                 dvap * (-12 * pow(c, 3) - 6 * pow(c, 2) * (12 * c - 15) -
                                         6 * c * (6 * pow(c, 2) - 15 * c + 10)));
    }
  }
};

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
  DEval eval_;  ///< F, gradient and Hessian (host and device)

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
  D() = default;
  explicit D(const double prefactor) { eval_.prefactor_ = prefactor; }
  virtual ~D() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const DEval& get_eval() const { return eval_; }
};

/**
 * @brief Value of D (see DEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
D::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of D (see DEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
D::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of D (see DEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
D::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
    return hessian;
  };
  return func;
}
/**
 * @brief Evaluation of Gimp at one point, on the host or on the device
 *
 * Same expressions as Gimp (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct GimpEval {
  static constexpr int nb_variables = 10;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 0;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;
  double A = 16.0;
  double B = 1.0;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    [[maybe_unused]] const double c = values[0];
    [[maybe_unused]] const double eta1 = values[1];
    [[maybe_unused]] const double eta2 = values[2];
    [[maybe_unused]] const double eta3 = values[3];
    [[maybe_unused]] const double eta4 = values[4];
    [[maybe_unused]] const double eta5 = values[5];
    [[maybe_unused]] const double eta6 = values[6];
    [[maybe_unused]] const double eta7 = values[7];
    [[maybe_unused]] const double eta8 = values[8];
    [[maybe_unused]] const double eta9 = values[9];
    if (F != nullptr) {
      *F = prefactor_ *
           (A * pow(c, 2) * pow(1 - c, 2) +
            B * (pow(c, 2) +
                 (6 - 6 * c) *
                     (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                      pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2)) -
                 (8 - 4 * c) *
                     (pow(eta1, 3) + pow(eta2, 3) + pow(eta3, 3) + pow(eta4, 3) + pow(eta5, 3) +
                      pow(eta6, 3) + pow(eta7, 3) + pow(eta8, 3) + pow(eta9, 3)) +
                 3 * pow(pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                             pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2),
                         2)));
    }
    if (gradient != nullptr) {
      gradient[0] =
          prefactor_ *
          (A * pow(c, 2) * (2 * c - 2) + 2 * A * c * pow(1 - c, 2) +
           B * (2 * c + 4 * pow(eta1, 3) - 6 * pow(eta1, 2) + 4 * pow(eta2, 3) - 6 * pow(eta2, 2) +
                4 * pow(eta3, 3) - 6 * pow(eta3, 2) + 4 * pow(eta4, 3) - 6 * pow(eta4, 2) +
                4 * pow(eta5, 3) - 6 * pow(eta5, 2) + 4 * pow(eta6, 3) - 6 * pow(eta6, 2) +
                4 * pow(eta7, 3) - 6 * pow(eta7, 2) + 4 * pow(eta8, 3) - 6 * pow(eta8, 2) +
                4 * pow(eta9, 3) - 6 * pow(eta9, 2)));
      gradient[1] =
          prefactor_ *
          (B * (-3 * pow(eta1, 2) * (8 - 4 * c) + 2 * eta1 * (6 - 6 * c) +
                12 * eta1 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[2] =
          prefactor_ *
          (B * (-3 * pow(eta2, 2) * (8 - 4 * c) + 2 * eta2 * (6 - 6 * c) +
                12 * eta2 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[3] =
          prefactor_ *
          (B * (-3 * pow(eta3, 2) * (8 - 4 * c) + 2 * eta3 * (6 - 6 * c) +
                12 * eta3 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[4] =
          prefactor_ *
          (B * (-3 * pow(eta4, 2) * (8 - 4 * c) + 2 * eta4 * (6 - 6 * c) +
                12 * eta4 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[5] =
          prefactor_ *
          (B * (-3 * pow(eta5, 2) * (8 - 4 * c) + 2 * eta5 * (6 - 6 * c) +
                12 * eta5 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[6] =
          prefactor_ *
          (B * (-3 * pow(eta6, 2) * (8 - 4 * c) + 2 * eta6 * (6 - 6 * c) +
                12 * eta6 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[7] =
          prefactor_ *
          (B * (-3 * pow(eta7, 2) * (8 - 4 * c) + 2 * eta7 * (6 - 6 * c) +
                12 * eta7 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[8] =
          prefactor_ *
          (B * (-3 * pow(eta8, 2) * (8 - 4 * c) + 2 * eta8 * (6 - 6 * c) +
                12 * eta8 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[9] =
          prefactor_ *
          (B * (-3 * pow(eta9, 2) * (8 - 4 * c) + 2 * eta9 * (6 - 6 * c) +
                12 * eta9 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ *
                   (2 * A * pow(c, 2) + 4 * A * c * (2 * c - 2) + 2 * A * pow(1 - c, 2) + 2 * B);
      hessian[1] = prefactor_ * (B * (12 * pow(eta1, 2) - 12 * eta1));
      hessian[2] = prefactor_ * (B * (12 * pow(eta2, 2) - 12 * eta2));
      hessian[3] = prefactor_ * (B * (12 * pow(eta3, 2) - 12 * eta3));
      hessian[4] = prefactor_ * (B * (12 * pow(eta4, 2) - 12 * eta4));
      hessian[5] = prefactor_ * (B * (12 * pow(eta5, 2) - 12 * eta5));
      hessian[6] = prefactor_ * (B * (12 * pow(eta6, 2) - 12 * eta6));
      hessian[7] = prefactor_ * (B * (12 * pow(eta7, 2) - 12 * eta7));
      hessian[8] = prefactor_ * (B * (12 * pow(eta8, 2) - 12 * eta8));
      hessian[9] = prefactor_ * (B * (12 * pow(eta9, 2) - 12 * eta9));
      hessian[10] = hessian[1];
      hessian[11] = prefactor_ * (B * (-12 * c + 36 * pow(eta1, 2) - 6 * eta1 * (8 - 4 * c) +
                                       12 * pow(eta2, 2) + 12 * pow(eta3, 2) + 12 * pow(eta4, 2) +
                                       12 * pow(eta5, 2) + 12 * pow(eta6, 2) + 12 * pow(eta7, 2) +
                                       12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[12] = prefactor_ * (24 * B * eta1 * eta2);
      hessian[13] = prefactor_ * (24 * B * eta1 * eta3);
      hessian[14] = prefactor_ * (24 * B * eta1 * eta4);
      hessian[15] = prefactor_ * (24 * B * eta1 * eta5);
      hessian[16] = prefactor_ * (24 * B * eta1 * eta6);
      hessian[17] = prefactor_ * (24 * B * eta1 * eta7);
      hessian[18] = prefactor_ * (24 * B * eta1 * eta8);
      hessian[19] = prefactor_ * (24 * B * eta1 * eta9);
      hessian[20] = hessian[2];
      hessian[21] = hessian[12];
      hessian[22] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 36 * pow(eta2, 2) - 6 * eta2 * (8 - 4 * c) +
                12 * pow(eta3, 2) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[23] = prefactor_ * (24 * B * eta2 * eta3);
      hessian[24] = prefactor_ * (24 * B * eta2 * eta4);
      hessian[25] = prefactor_ * (24 * B * eta2 * eta5);
      hessian[26] = prefactor_ * (24 * B * eta2 * eta6);
      hessian[27] = prefactor_ * (24 * B * eta2 * eta7);
      hessian[28] = prefactor_ * (24 * B * eta2 * eta8);
      hessian[29] = prefactor_ * (24 * B * eta2 * eta9);
      hessian[30] = hessian[3];
      hessian[31] = hessian[13];
      hessian[32] = hessian[23];
      hessian[33] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 36 * pow(eta3, 2) -
                6 * eta3 * (8 - 4 * c) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[34] = prefactor_ * (24 * B * eta3 * eta4);
      hessian[35] = prefactor_ * (24 * B * eta3 * eta5);
      hessian[36] = prefactor_ * (24 * B * eta3 * eta6);
      hessian[37] = prefactor_ * (24 * B * eta3 * eta7);
      hessian[38] = prefactor_ * (24 * B * eta3 * eta8);
      hessian[39] = prefactor_ * (24 * B * eta3 * eta9);
      hessian[40] = hessian[4];
      hessian[41] = hessian[14];
      hessian[42] = hessian[24];
      hessian[43] = hessian[34];
      hessian[44] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                36 * pow(eta4, 2) - 6 * eta4 * (8 - 4 * c) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[45] = prefactor_ * (24 * B * eta4 * eta5);
      hessian[46] = prefactor_ * (24 * B * eta4 * eta6);
      hessian[47] = prefactor_ * (24 * B * eta4 * eta7);
      hessian[48] = prefactor_ * (24 * B * eta4 * eta8);
      hessian[49] = prefactor_ * (24 * B * eta4 * eta9);
      hessian[50] = hessian[5];
      hessian[51] = hessian[15];
      hessian[52] = hessian[25];
      hessian[53] = hessian[35];
      hessian[54] = hessian[45];
      hessian[55] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                12 * pow(eta4, 2) + 36 * pow(eta5, 2) - 6 * eta5 * (8 - 4 * c) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[56] = prefactor_ * (24 * B * eta5 * eta6);
      hessian[57] = prefactor_ * (24 * B * eta5 * eta7);
      hessian[58] = prefactor_ * (24 * B * eta5 * eta8);
      hessian[59] = prefactor_ * (24 * B * eta5 * eta9);
      hessian[60] = hessian[6];
      hessian[61] = hessian[16];
      hessian[62] = hessian[26];
      hessian[63] = hessian[36];
      hessian[64] = hessian[46];
      hessian[65] = hessian[56];
      hessian[66] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 36 * pow(eta6, 2) - 6 * eta6 * (8 - 4 * c) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[67] = prefactor_ * (24 * B * eta6 * eta7);
      hessian[68] = prefactor_ * (24 * B * eta6 * eta8);
      hessian[69] = prefactor_ * (24 * B * eta6 * eta9);
      hessian[70] = hessian[7];
      hessian[71] = hessian[17];
      hessian[72] = hessian[27];
      hessian[73] = hessian[37];
      hessian[74] = hessian[47];
      hessian[75] = hessian[57];
      hessian[76] = hessian[67];
      hessian[77] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) + 36 * pow(eta7, 2) -
                6 * eta7 * (8 - 4 * c) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[78] = prefactor_ * (24 * B * eta7 * eta8);
      hessian[79] = prefactor_ * (24 * B * eta7 * eta9);
      hessian[80] = hessian[8];
      hessian[81] = hessian[18];
      hessian[82] = hessian[28];
      hessian[83] = hessian[38];
      hessian[84] = hessian[48];
      hessian[85] = hessian[58];
      hessian[86] = hessian[68];
      hessian[87] = hessian[78];
      hessian[88] = prefactor_ * (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) +
                                       12 * pow(eta3, 2) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) +
                                       12 * pow(eta6, 2) + 12 * pow(eta7, 2) + 36 * pow(eta8, 2) -
                                       6 * eta8 * (8 - 4 * c) + 12 * pow(eta9, 2) + 12));
      hessian[89] = prefactor_ * (24 * B * eta8 * eta9);
      hessian[90] = hessian[9];
      hessian[91] = hessian[19];
      hessian[92] = hessian[29];
      hessian[93] = hessian[39];
      hessian[94] = hessian[49];
      hessian[95] = hessian[59];
      hessian[96] = hessian[69];
      hessian[97] = hessian[79];
      hessian[98] = hessian[89];
      hessian[99] = prefactor_ * (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) +
                                       12 * pow(eta3, 2) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) +
                                       12 * pow(eta6, 2) + 12 * pow(eta7, 2) + 12 * pow(eta8, 2) +
                                       36 * pow(eta9, 2) - 6 * eta9 * (8 - 4 * c) + 12));
    }
  }
};

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
  GimpEval eval_;  ///< F, gradient and Hessian (host and device)

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
  Gimp() = default;
  explicit Gimp(const double prefactor) { eval_.prefactor_ = prefactor; }
  virtual ~Gimp() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const GimpEval& get_eval() const { return eval_; }
};

/**
 * @brief Value of Gimp (see GimpEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Gimp::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of Gimp (see GimpEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gimp::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(10);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of Gimp (see GimpEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gimp::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(100);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
    return hessian;
  };
  return func;
}
/**
 * @brief Evaluation of Gc at one point, on the host or on the device
 *
 * Same expressions as Gc (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct GcEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 9;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;
  double A = 16.0;
  double B = 1.0;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    [[maybe_unused]] const double c = values[0];
    [[maybe_unused]] const double eta1 = auxiliary_values[0];
    [[maybe_unused]] const double eta2 = auxiliary_values[1];
    [[maybe_unused]] const double eta3 = auxiliary_values[2];
    [[maybe_unused]] const double eta4 = auxiliary_values[3];
    [[maybe_unused]] const double eta5 = auxiliary_values[4];
    [[maybe_unused]] const double eta6 = auxiliary_values[5];
    [[maybe_unused]] const double eta7 = auxiliary_values[6];
    [[maybe_unused]] const double eta8 = auxiliary_values[7];
    [[maybe_unused]] const double eta9 = auxiliary_values[8];
    if (F != nullptr) {
      *F = prefactor_ *
           (A * pow(c, 2) * pow(1 - c, 2) +
            B * (pow(c, 2) +
                 (6 - 6 * c) *
                     (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                      pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2)) -
                 (8 - 4 * c) *
                     (pow(eta1, 3) + pow(eta2, 3) + pow(eta3, 3) + pow(eta4, 3) + pow(eta5, 3) +
                      pow(eta6, 3) + pow(eta7, 3) + pow(eta8, 3) + pow(eta9, 3)) +
                 3 * pow(pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                             pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2),
                         2)));
    }
    if (gradient != nullptr) {
      gradient[0] =
          prefactor_ *
          (A * pow(c, 2) * (2 * c - 2) + 2 * A * c * pow(1 - c, 2) +
           B * (2 * c + 4 * pow(eta1, 3) - 6 * pow(eta1, 2) + 4 * pow(eta2, 3) - 6 * pow(eta2, 2) +
                4 * pow(eta3, 3) - 6 * pow(eta3, 2) + 4 * pow(eta4, 3) - 6 * pow(eta4, 2) +
                4 * pow(eta5, 3) - 6 * pow(eta5, 2) + 4 * pow(eta6, 3) - 6 * pow(eta6, 2) +
                4 * pow(eta7, 3) - 6 * pow(eta7, 2) + 4 * pow(eta8, 3) - 6 * pow(eta8, 2) +
                4 * pow(eta9, 3) - 6 * pow(eta9, 2)));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ *
                   (2 * A * pow(c, 2) + 4 * A * c * (2 * c - 2) + 2 * A * pow(1 - c, 2) + 2 * B);
    }
  }
};

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
  GcEval eval_;  ///< F, gradient and Hessian (host and device)

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
  Gc() = default;
  explicit Gc(const double prefactor) { eval_.prefactor_ = prefactor; }
  virtual ~Gc() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const GcEval& get_eval() const { return eval_; }
};

/**
 * @brief Value of Gc (see GcEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Gc::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of Gc (see GcEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gc::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of Gc (see GcEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Gc::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
    return hessian;
  };
  return func;
}
/**
 * @brief Evaluation of Geta at one point, on the host or on the device
 *
 * Same expressions as Geta (F, gradient and Hessian), without memory
 * allocation. Each output is computed only if its pointer is not null.
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
struct GetaEval {
  static constexpr int nb_variables = 9;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 1;
  static constexpr int nb_coordinates = 0;

  double prefactor_ = 1.0;
  double A = 16.0;
  double B = 1.0;

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   [[maybe_unused]] const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    [[maybe_unused]] const double eta1 = values[0];
    [[maybe_unused]] const double eta2 = values[1];
    [[maybe_unused]] const double eta3 = values[2];
    [[maybe_unused]] const double eta4 = values[3];
    [[maybe_unused]] const double eta5 = values[4];
    [[maybe_unused]] const double eta6 = values[5];
    [[maybe_unused]] const double eta7 = values[6];
    [[maybe_unused]] const double eta8 = values[7];
    [[maybe_unused]] const double eta9 = values[8];
    [[maybe_unused]] const double c = auxiliary_values[0];
    if (F != nullptr) {
      *F = prefactor_ *
           (A * pow(c, 2) * pow(1 - c, 2) +
            B * (pow(c, 2) +
                 (6 - 6 * c) *
                     (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                      pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2)) -
                 (8 - 4 * c) *
                     (pow(eta1, 3) + pow(eta2, 3) + pow(eta3, 3) + pow(eta4, 3) + pow(eta5, 3) +
                      pow(eta6, 3) + pow(eta7, 3) + pow(eta8, 3) + pow(eta9, 3)) +
                 3 * pow(pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                             pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2),
                         2)));
    }
    if (gradient != nullptr) {
      gradient[0] =
          prefactor_ *
          (B * (-3 * pow(eta1, 2) * (8 - 4 * c) + 2 * eta1 * (6 - 6 * c) +
                12 * eta1 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[1] =
          prefactor_ *
          (B * (-3 * pow(eta2, 2) * (8 - 4 * c) + 2 * eta2 * (6 - 6 * c) +
                12 * eta2 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[2] =
          prefactor_ *
          (B * (-3 * pow(eta3, 2) * (8 - 4 * c) + 2 * eta3 * (6 - 6 * c) +
                12 * eta3 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[3] =
          prefactor_ *
          (B * (-3 * pow(eta4, 2) * (8 - 4 * c) + 2 * eta4 * (6 - 6 * c) +
                12 * eta4 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[4] =
          prefactor_ *
          (B * (-3 * pow(eta5, 2) * (8 - 4 * c) + 2 * eta5 * (6 - 6 * c) +
                12 * eta5 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[5] =
          prefactor_ *
          (B * (-3 * pow(eta6, 2) * (8 - 4 * c) + 2 * eta6 * (6 - 6 * c) +
                12 * eta6 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[6] =
          prefactor_ *
          (B * (-3 * pow(eta7, 2) * (8 - 4 * c) + 2 * eta7 * (6 - 6 * c) +
                12 * eta7 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[7] =
          prefactor_ *
          (B * (-3 * pow(eta8, 2) * (8 - 4 * c) + 2 * eta8 * (6 - 6 * c) +
                12 * eta8 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
      gradient[8] =
          prefactor_ *
          (B * (-3 * pow(eta9, 2) * (8 - 4 * c) + 2 * eta9 * (6 - 6 * c) +
                12 * eta9 *
                    (pow(eta1, 2) + pow(eta2, 2) + pow(eta3, 2) + pow(eta4, 2) + pow(eta5, 2) +
                     pow(eta6, 2) + pow(eta7, 2) + pow(eta8, 2) + pow(eta9, 2))));
    }
    if (hessian != nullptr) {
      hessian[0] = prefactor_ * (B * (-12 * c + 36 * pow(eta1, 2) - 6 * eta1 * (8 - 4 * c) +
                                      12 * pow(eta2, 2) + 12 * pow(eta3, 2) + 12 * pow(eta4, 2) +
                                      12 * pow(eta5, 2) + 12 * pow(eta6, 2) + 12 * pow(eta7, 2) +
                                      12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[1] = prefactor_ * (24 * B * eta1 * eta2);
      hessian[2] = prefactor_ * (24 * B * eta1 * eta3);
      hessian[3] = prefactor_ * (24 * B * eta1 * eta4);
      hessian[4] = prefactor_ * (24 * B * eta1 * eta5);
      hessian[5] = prefactor_ * (24 * B * eta1 * eta6);
      hessian[6] = prefactor_ * (24 * B * eta1 * eta7);
      hessian[7] = prefactor_ * (24 * B * eta1 * eta8);
      hessian[8] = prefactor_ * (24 * B * eta1 * eta9);
      hessian[9] = hessian[1];
      hessian[10] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 36 * pow(eta2, 2) - 6 * eta2 * (8 - 4 * c) +
                12 * pow(eta3, 2) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[11] = prefactor_ * (24 * B * eta2 * eta3);
      hessian[12] = prefactor_ * (24 * B * eta2 * eta4);
      hessian[13] = prefactor_ * (24 * B * eta2 * eta5);
      hessian[14] = prefactor_ * (24 * B * eta2 * eta6);
      hessian[15] = prefactor_ * (24 * B * eta2 * eta7);
      hessian[16] = prefactor_ * (24 * B * eta2 * eta8);
      hessian[17] = prefactor_ * (24 * B * eta2 * eta9);
      hessian[18] = hessian[2];
      hessian[19] = hessian[11];
      hessian[20] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 36 * pow(eta3, 2) -
                6 * eta3 * (8 - 4 * c) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[21] = prefactor_ * (24 * B * eta3 * eta4);
      hessian[22] = prefactor_ * (24 * B * eta3 * eta5);
      hessian[23] = prefactor_ * (24 * B * eta3 * eta6);
      hessian[24] = prefactor_ * (24 * B * eta3 * eta7);
      hessian[25] = prefactor_ * (24 * B * eta3 * eta8);
      hessian[26] = prefactor_ * (24 * B * eta3 * eta9);
      hessian[27] = hessian[3];
      hessian[28] = hessian[12];
      hessian[29] = hessian[21];
      hessian[30] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                36 * pow(eta4, 2) - 6 * eta4 * (8 - 4 * c) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[31] = prefactor_ * (24 * B * eta4 * eta5);
      hessian[32] = prefactor_ * (24 * B * eta4 * eta6);
      hessian[33] = prefactor_ * (24 * B * eta4 * eta7);
      hessian[34] = prefactor_ * (24 * B * eta4 * eta8);
      hessian[35] = prefactor_ * (24 * B * eta4 * eta9);
      hessian[36] = hessian[4];
      hessian[37] = hessian[13];
      hessian[38] = hessian[22];
      hessian[39] = hessian[31];
      hessian[40] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                12 * pow(eta4, 2) + 36 * pow(eta5, 2) - 6 * eta5 * (8 - 4 * c) + 12 * pow(eta6, 2) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[41] = prefactor_ * (24 * B * eta5 * eta6);
      hessian[42] = prefactor_ * (24 * B * eta5 * eta7);
      hessian[43] = prefactor_ * (24 * B * eta5 * eta8);
      hessian[44] = prefactor_ * (24 * B * eta5 * eta9);
      hessian[45] = hessian[5];
      hessian[46] = hessian[14];
      hessian[47] = hessian[23];
      hessian[48] = hessian[32];
      hessian[49] = hessian[41];
      hessian[50] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 36 * pow(eta6, 2) - 6 * eta6 * (8 - 4 * c) +
                12 * pow(eta7, 2) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[51] = prefactor_ * (24 * B * eta6 * eta7);
      hessian[52] = prefactor_ * (24 * B * eta6 * eta8);
      hessian[53] = prefactor_ * (24 * B * eta6 * eta9);
      hessian[54] = hessian[6];
      hessian[55] = hessian[15];
      hessian[56] = hessian[24];
      hessian[57] = hessian[33];
      hessian[58] = hessian[42];
      hessian[59] = hessian[51];
      hessian[60] =
          prefactor_ *
          (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) + 12 * pow(eta3, 2) +
                12 * pow(eta4, 2) + 12 * pow(eta5, 2) + 12 * pow(eta6, 2) + 36 * pow(eta7, 2) -
                6 * eta7 * (8 - 4 * c) + 12 * pow(eta8, 2) + 12 * pow(eta9, 2) + 12));
      hessian[61] = prefactor_ * (24 * B * eta7 * eta8);
      hessian[62] = prefactor_ * (24 * B * eta7 * eta9);
      hessian[63] = hessian[7];
      hessian[64] = hessian[16];
      hessian[65] = hessian[25];
      hessian[66] = hessian[34];
      hessian[67] = hessian[43];
      hessian[68] = hessian[52];
      hessian[69] = hessian[61];
      hessian[70] = prefactor_ * (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) +
                                       12 * pow(eta3, 2) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) +
                                       12 * pow(eta6, 2) + 12 * pow(eta7, 2) + 36 * pow(eta8, 2) -
                                       6 * eta8 * (8 - 4 * c) + 12 * pow(eta9, 2) + 12));
      hessian[71] = prefactor_ * (24 * B * eta8 * eta9);
      hessian[72] = hessian[8];
      hessian[73] = hessian[17];
      hessian[74] = hessian[26];
      hessian[75] = hessian[35];
      hessian[76] = hessian[44];
      hessian[77] = hessian[53];
      hessian[78] = hessian[62];
      hessian[79] = hessian[71];
      hessian[80] = prefactor_ * (B * (-12 * c + 12 * pow(eta1, 2) + 12 * pow(eta2, 2) +
                                       12 * pow(eta3, 2) + 12 * pow(eta4, 2) + 12 * pow(eta5, 2) +
                                       12 * pow(eta6, 2) + 12 * pow(eta7, 2) + 12 * pow(eta8, 2) +
                                       36 * pow(eta9, 2) - 6 * eta9 * (8 - 4 * c) + 12));
    }
  }
};

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
  GetaEval eval_;  ///< F, gradient and Hessian (host and device)

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
  Geta() = default;
  explicit Geta(const double prefactor) { eval_.prefactor_ = prefactor; }
  virtual ~Geta() = default;
  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points
  const GetaEval& get_eval() const { return eval_; }
};

/**
 * @brief Value of Geta (see GetaEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Geta::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of Geta (see GetaEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Geta::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(9);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of Geta (see GetaEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Geta::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(81);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
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
