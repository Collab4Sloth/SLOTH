/**
 * @file PhysicalConvergence.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to define a PhysicalConvergence criterion (physical increment allowed per
 * time-step)
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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

#include "Convergence/PhysicalConvergence.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Physical Convergence:: Physical Convergence object
 *
 * @param convergence_criterion_type
 * @param given_criterion
 */
PhysicalConvergence::PhysicalConvergence(ConvergenceType::value convergence_criterion_type,
                                         const double& given_criterion) {
  this->physical_convergence_ =
      this->getPhysicalConvergenceCriterion(convergence_criterion_type, given_criterion);
}

std::tuple<bool, double> PhysicalConvergence::getPhysicalConvergence(const mfem::Vector& x,
                                                                     const mfem::Vector& y) {
  return this->physical_convergence_(x, y);
}

/**
 * @brief Return the function used to assess the convergence by comparing two Vector with given
 * convergence function and criterion
 *
 * @param convergence_criterion_type
 * @param given_criterion
 * @return std::function<double(const mfem::Vector &, const mfem::Vector &)>
 */
std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>
PhysicalConvergence::getPhysicalConvergenceCriterion(
    ConvergenceType::value convergence_criterion_type, const double& given_criterion) {
  switch (convergence_criterion_type) {
    case ConvergenceType::RELATIVE_MAX: {
      return this->getRelativeMax(given_criterion);
    }
    case ConvergenceType::ABSOLUTE_MAX: {
      return this->getAbsoluteMax(given_criterion);
    }
    default:
      mfem::mfem_error(
          "PhysicalConvergence::getPhysicalConvergenceCriterion:  RELATIVE_MAX, ABSOLUTE_MAX "
          "criterion are  available");
      break;
  }
}

/**
 * @brief Return the boolean that indicates if convergence is reached or not, with the calculated
 * criterion using RELATIVE_MAX method
 *
 * @param given_criterion
 * @return std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
 */
std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>
PhysicalConvergence::getRelativeMax(const double& given_criterion) {
  return std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>(
      [given_criterion](const mfem::Vector& x, const mfem::Vector& xn) {
        // Maximum of the relative differences, element per element (0 for empty vectors,
        // e.g. a rank without dof: the criterion is then reduced over the ranks)
        constexpr double epsilon = 1.e-12;
        double criterion = 0.;
        for (int i = 0; i < x.Size(); ++i) {
          criterion =
              std::max(criterion, std::abs(x(i) - xn(i)) / std::max(epsilon, std::abs(xn(i))));
        }
        return std::make_tuple(criterion < given_criterion, criterion);
      });
}

/**
 * @brief  Return the boolean that indicates if convergence is reached or not, with the calculated
 * criterion using ABSOLUTE_MAX method
 *
 * @param given_criterion
 * @return std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
 */
std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>
PhysicalConvergence::getAbsoluteMax(const double& given_criterion) {
  return std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>(
      [given_criterion](const mfem::Vector& x, const mfem::Vector& xn) {
        // Maximum of the absolute differences, element per element (0 for empty vectors,
        // e.g. a rank without dof: the criterion is then reduced over the ranks)
        double criterion = 0.;
        for (int i = 0; i < x.Size(); ++i) {
          criterion = std::max(criterion, std::abs(x(i) - xn(i)));
        }
        return std::make_tuple(criterion < given_criterion, criterion);
      });
}

/**
 * @brief Destroy the Physical Convergence:: Physical Convergence object
 *
 */
PhysicalConvergence::~PhysicalConvergence() {}
