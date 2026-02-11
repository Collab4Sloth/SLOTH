/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Covering test for SlothBaseCoefficients
 * @version 0.1
 * @date 2025-11-17
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
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "./Function.hpp"
#include "./OtherFunction.hpp"
#include "kernel/Coefficients/Coefficient.hpp"
#include "kernel/Coefficients/Coefficients.hpp"
#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

struct TestParameters {
  int tcase = 1;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(
      &p.tcase, "-tc", "--test-case",
      "0: evaluation of variables, 1: gradient and hessian evaluations, 2: check the "
      "use of Coefficients, 3: check the sum of two coefficients, 4: check the products of "
      "two coefficients");

  args.Parse();

  if (!args.Good()) {
    args.PrintUsage(mfem::out);
    std::exit(EXIT_FAILURE);
  }
  args.PrintOptions(mfem::out);
}
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Verbose);

  // ################ //
  // ################ //
  //   Read options   //
  // ################ //
  // ################ //
  constexpr double epsilon = 1.e-10;

  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);
  switch (p.tcase) {
    case 0: {
      // Evaluation of variables
      SlothInfo::print("Running test case 0: variable evaluation");
      Coefficient coeff(Glossary::Temperature, Scheme::Implicit, FunctionA());
      try {
        auto result = coeff.compute({2, 3, 1});
        if (std::abs(result - 94.1) > epsilon) {
          throw std::runtime_error("Wrong variable evaluation");
        }
        std::cout << "Successfull variable evaluation: " << result << " " << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }
    case 1: {
      // Evaluation of gradient
      SlothInfo::print("Running test case 1: gradient/hessian evaluation");
      Coefficient coeff(Glossary::Temperature, Scheme::Implicit, FunctionB());
      std::vector<double> gradient_solution{98.3, 60.2, 0.6};
      std::vector<double> hessian_solution{4.0, 30.1, 0.3, 30.1, 0.0, 0.2, 0.3, 0.2, 0.0};
      try {
        for (int i = 0; i < 3; i++) {
          const double gradient_i = coeff.compute_gradient(i, {2, 3, 1});

          if (std::abs(gradient_i - gradient_solution[i]) > epsilon) {
            throw std::runtime_error("Wrong gradient evaluation");
          }

          for (int j = 0; j < 3; j++) {
            const double hessian_ij = coeff.compute_hessian(i, j, {2, 3, 1});

            if (std::abs(hessian_ij - hessian_solution[i * 3 + j]) > epsilon) {
              throw std::runtime_error("Wrong hessian evaluation");
            }
          }
        }
        std::cout << "Successfull gradient and hessian evaluation " << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }
    case 2: {
      // Check Coefficients
      SlothInfo::print("Running test case 2: gradient/hessian evaluation");
      Coefficient coeffA(Glossary::Temperature, Scheme::Implicit, FunctionA());
      Coefficient coeffB(Glossary::Temperature, Scheme::Implicit, FunctionB());
      Coefficient coeffC(Glossary::Temperature, Scheme::Implicit, FunctionB());
      std::vector<double> gradient_solution{98.3, 60.2, 0.6};
      std::vector<double> hessian_solution{4.0, 30.1, 0.3, 30.1, 0.0, 0.2, 0.3, 0.2, 0.0};
      Coefficients coeffAB = Coefficients(coeffA, coeffB);
      coeffAB.add(coeffC);

      try {
        if (coeffAB.size() != 3) {
          throw std::runtime_error("Wrong size of coefficients");
        }
        auto coeff = coeffAB[1];
        for (int i = 0; i < 3; i++) {
          const double gradient_i = coeff.compute_gradient(i, {2, 3, 1});

          if (std::abs(gradient_i - gradient_solution[i]) > epsilon) {
            throw std::runtime_error("Wrong gradient evaluation");
          }

          for (int j = 0; j < 3; j++) {
            const double hessian_ij = coeff.compute_hessian(i, j, {2, 3, 1});

            if (std::abs(hessian_ij - hessian_solution[i * 3 + j]) > epsilon) {
              throw std::runtime_error("Wrong hessian evaluation");
            }
          }
        }
        std::cout << "Successfull gradient and hessian evaluation " << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }
    case 3: {
      // Check sum of Coefficients
      SlothInfo::print("Running test case 3: sum evaluation");
      Coefficient coeffA(Glossary::Temperature, Scheme::Implicit, FunctionB());
      auto A = FunctionB();

      Coefficient coeffB(Glossary::Temperature, Scheme::Implicit, SumCoefficient({&A, &A}));
      Coefficients coeffAB = Coefficients(coeffB);
      std::vector<double> gradient_solution{98.3, 60.2, 0.6};
      std::vector<double> hessian_solution{4.0, 30.1, 0.3, 30.1, 0.0, 0.2, 0.3, 0.2, 0.0};

      try {
        auto coeff = coeffAB[0];
        for (int i = 0; i < 3; i++) {
          const double gradient_i = coeff.compute_gradient(i, {2, 3, 1});

          if (std::abs(gradient_i - 2.0 * gradient_solution[i]) > epsilon) {
            throw std::runtime_error("Wrong gradient evaluation");
          }

          for (int j = 0; j < 3; j++) {
            const double hessian_ij = coeff.compute_hessian(i, j, {2, 3, 1});

            if (std::abs(hessian_ij - 2.0 * hessian_solution[i * 3 + j]) > epsilon) {
              throw std::runtime_error("Wrong hessian evaluation");
            }
          }
        }
        std::cout << "Successfull evaluation of sum of coefficients " << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }

    case 4: {
      // Check product of Coefficients
      SlothInfo::print("Running test case 4: product evaluation");
      auto A = FunctionA();
      auto B = FunctionB();

      Coefficient coeffB(Glossary::Temperature, Scheme::Implicit, ProductCoefficient({&A, &B}));
      Coefficients coeffAB = Coefficients(coeffB);
      std::vector<double> gradient_solution{9627.23, 11322.82, 75.32};
      std::vector<double> hessian_solution{769.6, 5901.81, 39.26, 5901.81, 3612,
                                           42.84, 39.26,   42.84, 0.12};
      try {
        auto coeff = coeffAB[0];
        const double value = coeff.compute({2.0, 3.0, 1.0});
        if (std::abs(value - 17747.26) > epsilon) {
          throw std::runtime_error("Wrong value evaluation");
        }
        for (int i = 0; i < 3; i++) {
          const double gradient_i = coeff.compute_gradient(i, {2.0, 3.0, 1.0});

          if (std::abs(gradient_i - gradient_solution[i]) > epsilon) {
            throw std::runtime_error("Wrong gradient evaluation " +
                                     std::format("{:.14f}", gradient_i - gradient_solution[i]));
          }

          for (int j = 0; j < 3; j++) {
            const double hessian_ij = coeff.compute_hessian(i, j, {2, 3, 1});

            if (std::abs(hessian_ij - hessian_solution[i * 3 + j]) > epsilon) {
              throw std::runtime_error(
                  "Wrong hessian evaluation " +
                  std::format("{:.14f}", hessian_ij - hessian_solution[i * 3 + j]));
            }
          }
        }
        std::cout << "Successfull evaluation of product of coefficients" << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }

    case 5: {
      // Check sum of two product of Coefficients
      SlothInfo::print("Running test case 5:  sum of product evaluation");
      auto A = FunctionA();
      auto B = FunctionB();

      auto P1 = ProductCoefficient({&A, &B});
      auto S1 = SumCoefficient({&P1, &P1});
      Coefficient coeffB(Glossary::Temperature, Scheme::Implicit, S1);
      Coefficients coeffAB = Coefficients(coeffB);
      std::vector<double> gradient_solution{9627.23, 11322.82, 75.32};
      std::vector<double> hessian_solution{769.6, 5901.81, 39.26, 5901.81, 3612,
                                           42.84, 39.26,   42.84, 0.12};
      try {
        auto coeff = coeffAB[0];
        const double value = coeff.compute({2.0, 3.0, 1.0});
        if (std::abs(value - 2.0 * 17747.26) > epsilon) {
          throw std::runtime_error("Wrong value evaluation " + std::format("{:.14f}", value));
        }
        for (int i = 0; i < 3; i++) {
          const double gradient_i = coeff.compute_gradient(i, {2.0, 3.0, 1.0});

          if (std::abs(gradient_i - 2.0 * gradient_solution[i]) > epsilon) {
            throw std::runtime_error(
                "Wrong gradient evaluation " +
                std::format("{:.14f}", gradient_i - 2.0 * gradient_solution[i]));
          }

          for (int j = 0; j < 3; j++) {
            const double hessian_ij = coeff.compute_hessian(i, j, {2, 3, 1});

            if (std::abs(hessian_ij - 2.0 * hessian_solution[i * 3 + j]) > epsilon) {
              throw std::runtime_error(
                  "Wrong hessian evaluation " +
                  std::format("{:.14f}", hessian_ij - 2.0 * hessian_solution[i * 3 + j]));
            }
          }
        }
        std::cout << "Successfull evaluation of sum of product of coefficients" << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }
  }

  return 0;
}
