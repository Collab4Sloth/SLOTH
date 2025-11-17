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

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

struct TestParameters {
  int tcase = 0;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.tcase, "-tc", "--test-case",
                 "0: evaluation of variables, 1: wrong number of arguments, 2: wrong expression, "
                 "3: check GlossaryType, 4: constant evaluation");

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
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);
  switch (p.tcase) {
    case 0: {
      // Evaluation of variables
      SlothInfo::print("Running test case 0: variable evaluation");
      SlothBaseCoefficient<double> coeff(Glossary::Temperature, "2*x+30*y+0.1*z", "x", "y", "z");
      auto result = coeff.evaluate(2, 3, 1);
      try {
        if (result != 94.1) {
          throw std::runtime_error("Wrong variable evaluation");
        }
        std::cout << "Successfull variable evaluation: " << result << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }

      break;
    }
    case 1: {
      // Wrong number of variables
      SlothInfo::print("Running test case 2: wrong number of arguments");
      std::vector<std::string> variable_names = {"x", "y", "z"};
      SlothBaseCoefficient<double> coeff(Glossary::Temperature, "2*x+30*y+0.1*z", "x", "y", "z");
      try {
        auto bad_result = coeff.evaluate(1.0, 2.0);
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }

      break;
    }
    case 2: {
      // Wrong expression
      SlothInfo::print("Running test case 3: wrong expression");

      std::vector<std::string> variable_names = {"x", "y", "z"};
      try {
        SlothBaseCoefficient<double> bad_coeff(Glossary::Temperature, "2*x+30*y+0.1*z*", "x", "y",
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

      std::vector<std::string> variable_names = {"x", "y", "z"};
      SlothBaseCoefficient<double> coeff(Glossary::Temperature, "2*x+30*y+0.1*z", "x", "y", "z");
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
      SlothBaseCoefficient<double> coeff(Glossary::Temperature, "1");
      auto result = coeff.evaluate();
      try {
        if (result != 1) {
          throw std::runtime_error("Wrong constant evaluation");
        }
        std::cout << "Successfull constant evaluation: " << result << std::endl;
      } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
      }

      break;
    }
  }

  return 0;
}
