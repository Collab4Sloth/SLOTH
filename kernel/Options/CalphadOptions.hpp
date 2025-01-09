/**
 * @file CalphadOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for Calphad problems
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <limits>
#include <string>
#include <vector>

#include "Utils/Utils.hpp"

#pragma once

/**
 * @brief Available outputs for Calphad problems
 */
struct calphad_outputs {
  enum value { mu, g, gm, h, hm, x };
  static value from(const std::string&);
};

calphad_outputs::value calphad_outputs::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<calphad_outputs::value> m{
      {"mu", calphad_outputs::mu}, {"x", calphad_outputs::x}, {"g", calphad_outputs::g},
      {"gm", calphad_outputs::gm}, {"h", calphad_outputs::h}, {"hm", calphad_outputs::hm}};
  return m.find("calphad_outputs", v);
}
