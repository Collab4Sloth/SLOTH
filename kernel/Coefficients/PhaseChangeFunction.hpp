/**
 * @file PhaseChangeFunction.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief PhaseChange Functions
 * @version 0.1
 * @date 2025-09-05
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

#ifndef PhaseChangeFunction_HPP_
#define PhaseChangeFunction_HPP_
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "Utils/Utils.hpp"

template <PhaseChange PHASECHANGE>
class PhaseChangeFunction {
 private:
  template <class... Args>
  std::function<double(const double&)> getConstant(Args... args);

  template <class... Args>
  std::function<double(const double&)> getCalphad(Args... args);

 public:
  PhaseChangeFunction();

  template <class... Args>
  std::function<double(const double&)> getPhaseChangeFunction(Args... args);

  ~PhaseChangeFunction();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Phase Change Functions<PHASECHANGE>:: Phase Change Functions object
 *
 * @tparam PHASECHANGE
 */
template <PhaseChange PHASECHANGE>
PhaseChangeFunction<PHASECHANGE>::PhaseChangeFunction() {}

/**
 * @brief Get a constant PhaseChange
 *
 * @tparam PHASECHANGE
 * @tparam Args
 * @param args
 * @return std::function<double(const double&)>
 */
template <PhaseChange PHASECHANGE>
template <typename... Args>
std::function<double(const double&)> PhaseChangeFunction<PHASECHANGE>::getConstant(Args... args) {
  return std::function<double(const double&)>([](double x) { return x; });
}

/**
 * @brief Get the PhaseChange from CALPHAD
 *
 * @tparam PHASECHANGE
 * @tparam Args
 */
template <PhaseChange PHASECHANGE>
template <typename... Args>
std::function<double(const double&)> PhaseChangeFunction<PHASECHANGE>::getCalphad(Args... args) {
  mfem::mfem_error("PhaseChangeFunction::getCalphad: not implemented yet");
}

/**
 * @brief Get the PhaseChange function
 *
 * @tparam PHASECHANGE
 * @tparam Args
 * @param args
 * @return std::function<double(const double&)>
 */
template <PhaseChange PHASECHANGE>
template <class... Args>
std::function<double(const double&)> PhaseChangeFunction<PHASECHANGE>::getPhaseChangeFunction(
    Args... args) {
  switch (PHASECHANGE) {
    case PhaseChange::Constant:
      return this->getConstant(args...);
    case PhaseChange::Calphad:
      return this->getCalphad(args...);
    default:
      mfem::mfem_error(
          "PhaseChangeFunction::getPhaseChangeFunction: only constant and calphad mobilities "
          " are available");
      break;
  }
}

/**
 * @brief Destroy the Phase Change Functions< PHASECHANGE>:: Phase Change Functions object
 *
 * @tparam PHASECHANGE
 */
template <PhaseChange PHASECHANGE>
PhaseChangeFunction<PHASECHANGE>::~PhaseChangeFunction() {}

#endif
