/**
 * @file PhysicalConvergence.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to define a PhysicalConvergence criterion (physical increment allowed per
 * time-step)
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

#include "Options/Options.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

class PhysicalConvergence {
 private:
  std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>
  getPhysicalConvergenceCriterion(ConvergenceType::value convergence_criterion_type,
                                  const double& given_criterion);
  std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)> getAbsoluteMax(
      const double& given_criterion);
  std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)> getRelativeMax(
      const double& given_criterion);

  std::function<std::tuple<bool, double>(const mfem::Vector&, const mfem::Vector&)>
      physical_convergence_;

 public:
  PhysicalConvergence(ConvergenceType::value convergence_criterion_type,
                      const double& given_criterion);
  std::tuple<bool, double> getPhysicalConvergence(const mfem::Vector& x, const mfem::Vector& y);
  ~PhysicalConvergence();
};
