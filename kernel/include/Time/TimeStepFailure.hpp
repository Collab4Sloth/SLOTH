/**
 * @file TimeStepFailure.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Exception thrown when a time-step cannot be completed
 * @version 0.1
 * @date 2026-10-27
 *
 * @copyright CEA (C) 2026
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

#include <stdexcept>
#include <string>

/**
 * @brief Failure of a time-step (e.g. Newton solver not converged)
 *
 * @details Thrown on every MPI rank: the failure is detected from collective quantities
 *          (residual norms reduced over all the ranks). Caught by TimeDiscretization, which
 *          either rejects the time-step and retries it with a smaller time-step (adaptive
 *          time-stepping) or stops the computation.
 */
class TimeStepFailure : public std::runtime_error {
 public:
  explicit TimeStepFailure(const std::string& message) : std::runtime_error(message) {}
};
