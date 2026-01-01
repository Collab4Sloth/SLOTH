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

#include "kernel/Coefficients/Coefficient.hpp"
#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

struct TestParameters {
  int tcase = 5;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(
      &p.tcase, "-tc", "--test-case",
      "0: evaluation of variables, 1: wrong number of arguments, 2: wrong expression, "
      "3: check GlossaryType, 4: constant evaluation, 5: gradient and hessian evaluations");

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
  constexpr double epsilon = 1.e-12;

  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);
  switch (p.tcase) {
    case 0: {
      // Evaluation of variables
      SlothInfo::print("Running test case 0: variable evaluation");
      Coefficient coeff(Glossary::Temperature, Scheme::Implicit, "2*x+30*y+0.1*z", "x", "y", "z");
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
      // Wrong number of variables
      SlothInfo::print("Running test case 2: wrong number of arguments");
      Coefficient coeff(Glossary::Temperature, Scheme::Implicit, "2*x+30*y+0.1*z", "x", "y", "z");
      try {
        const double bad_result = coeff.compute({1.0, 2.0});
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }

      break;
    }
    case 2: {
      // Wrong expression
      SlothInfo::print("Running test case 3: wrong expression");

      try {
        Coefficient coeff(Glossary::Temperature, Scheme::Implicit, "2*x+30*y+0.1*z*", "x", "y",
                          "z");
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }

      break;
    }
    case 3: {
      // Wrong expression
      SlothInfo::print("Running test case 4: check GlossaryType");

      Coefficient coeff(Glossary::Temperature, Scheme::Implicit, "2*x+30*y+0.1*z", "x", "y", "z");
      try {
        if (coeff.get_type() != GlossaryType::Temperature) {
          throw std::runtime_error("Wrong GlossaryType");
        }
        std::cout << "Successfull checking of the GlossaryType" << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }
      break;
    }
    case 4: {
      // Evaluation of a constant
      SlothInfo::print("Running test case 4: constant evaluation");
      Coefficient coeff(Glossary::Temperature, 1);
      try {
        const double result = coeff.compute();
        if (std::abs(result - 1) > epsilon) {
          throw std::runtime_error("Wrong constant evaluation");
        }
        std::cout << "Successfull constant evaluation: " << result << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }

      break;
    }
    case 5: {
      // Evaluation of gradient
      SlothInfo::print("Running test case 5: gradient/hessian evaluation");
      std::string function = "2*x*x+30*x*y+0.1*z*x*y";
      std::string gx = "4*x+30*y+0.1*z*y";
      std::string gy = "30*x+0.1*z*x";
      std::string gz = "0.1*x*y";
      std::vector<std::string> grad_functions{gx, gy, gz};
      std::string gxx = "4.0";
      std::string gxy = "30+0.1*z";
      std::string gxz = "0.1*y";
      std::string gyx = "30+0.1*z";
      std::string gyy = "0.0";
      std::string gyz = "0.1*x";
      std::string gzx = "0.1*y";
      std::string gzy = "0.1*x";
      std::string gzz = "0.0";
      std::vector<std::string> hess_functions{gxx, gxy, gxz, gyx, gyy, gyz, gzx, gyz, gzz};

      Coefficient coeff(Glossary::Temperature, Scheme::Implicit, hess_functions, grad_functions,
                        function, "x", "y", "z");
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
  }

  return 0;
}
