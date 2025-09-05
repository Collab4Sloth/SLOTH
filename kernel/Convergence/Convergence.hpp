/**
 * @file Convergence.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to build and manage PhysicalConvergence
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

#pragma once
#include <algorithm>
#include <limits>
#include <tuple>
#include <vector>

#include "Convergence/PhysicalConvergence.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
class Convergence {
 private:
  std::vector<PhysicalConvergence> vect_convergence_;

 public:
  template <class... Args>
  explicit Convergence(Args... args);

  Convergence();

  std::vector<PhysicalConvergence> getPhysicalConvergence() const;

  ~Convergence();
};

/**
 * @brief Construct a new Convergence:: Convergence object
 *
 */
Convergence::Convergence() {}

/**
 * @brief Construct a new Convergence:: Convergence object
 *
 * @tparam Args
 * @param args
 */
template <class... Args>
Convergence::Convergence(Args... args) {
  this->vect_convergence_ = std::vector<PhysicalConvergence>{args...};
}

/**
 * @brief Return the vector of PhysicalConvergence Objects
 *
 * @return std::vector<PhysicalConvergence>
 */
std::vector<PhysicalConvergence> Convergence::getPhysicalConvergence() const {
  return this->vect_convergence_;
}

/**
 * @brief Destroy the Convergence:: Convergence object
 *
 */
Convergence::~Convergence() {}
