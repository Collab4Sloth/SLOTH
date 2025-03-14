/**
 * @file ProblemsOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for Problem objects
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <string>

#include "Utils/Utils.hpp"

#pragma once

///////////////////////////////////////////////////
//////// Problems
///////////////////////////////////////////////////

struct Problems {
  enum value { Diffusion, AllenCahn, Calphad };
  static value from(const std::string&);
};

inline Problems::value Problems::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<Problems::value> m{{"Diffusion", Problems::Diffusion},
                                                    {"AllenCahn", Problems::AllenCahn},
                                                    {"Calphad", Problems::Calphad}};
  return m.find("Problems", v);
}
