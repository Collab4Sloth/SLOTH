/**
 * @file Convergence.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to build and manage PhysicalConvergence
 * @version 0.1
 * @date 2025-07-25
 *
 * Copyright CEA (c) 2025
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
