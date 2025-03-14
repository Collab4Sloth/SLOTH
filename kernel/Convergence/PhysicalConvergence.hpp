/**
 * @file PhysicalConvergence.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to define a PhysicalConvergence criterion (physical increment allowed per
 * time-step)
 * @version 0.1
 * @date 2024-05-27
 *
 * Copyright CEA (c) 2024
 *
 */

#pragma once
#include <algorithm>
#include <limits>
#include <tuple>
#include <vector>

#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

class PhysicalConvergence {
 private:
  std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
  getPhysicalConvergenceCriterion(ConvergenceType::value convergence_criterion_type,
                                  const double &given_criterion);
  std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
  getAbsoluteMax(const double &given_criterion);
  std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
  getRelativeMax(const double &given_criterion);

  std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
      physical_convergence_;

 public:
  PhysicalConvergence(ConvergenceType::value convergence_criterion_type,
                      const double &given_criterion);
  std::tuple<bool, double> getPhysicalConvergence(const mfem::Vector &x, const mfem::Vector &y);
  ~PhysicalConvergence() = default;
};

/**
 * @brief Construct a new Physical Convergence:: Physical Convergence object
 *
 * @param convergence_criterion_type
 * @param given_criterion
 */
DEBILE_INLINE PhysicalConvergence::PhysicalConvergence(ConvergenceType::value convergence_criterion_type,
                                         const double &given_criterion) {
  this->physical_convergence_ =
      this->getPhysicalConvergenceCriterion(convergence_criterion_type, given_criterion);
}

DEBILE_INLINE std::tuple<bool, double> PhysicalConvergence::getPhysicalConvergence(const mfem::Vector &x,
                                                                     const mfem::Vector &y) {
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
DEBILE_INLINE std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
PhysicalConvergence::getPhysicalConvergenceCriterion(
    ConvergenceType::value convergence_criterion_type, const double &given_criterion) {
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
DEBILE_INLINE std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
PhysicalConvergence::getRelativeMax(const double &given_criterion) {
  return std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>(
      [given_criterion](const mfem::Vector &x, const mfem::Vector &xn) {
        std::vector<double> abs_error;
        abs_error.resize(x.Size());
        const auto epsilon = 1.e-12;
        // Absolute difference element per element
        std::transform(x.begin(), x.end(), xn.begin(), abs_error.begin(),
                       [epsilon](double a, double b) {
                         return std::abs(a - b) / std::max(epsilon, std::abs(b));
                       });
        // Maximum among the absolute differences

        const auto &criterion = *(std::max_element(abs_error.begin(), abs_error.end()));
        // Checking if criterion is satisfied
        bool is_cvg = (criterion < given_criterion);
        abs_error.clear();
        return std::make_tuple(is_cvg, criterion);
      });
}

/**
 * @brief  Return the boolean that indicates if convergence is reached or not, with the calculated
 * criterion using ABSOLUTE_MAX method
 *
 * @param given_criterion
 * @return std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
 */
DEBILE_INLINE std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>
PhysicalConvergence::getAbsoluteMax(const double &given_criterion) {
  return std::function<std::tuple<bool, double>(const mfem::Vector &, const mfem::Vector &)>(
      [given_criterion](const mfem::Vector &x, const mfem::Vector &xn) {
        std::vector<double> abs_error;
        abs_error.resize(x.Size());
        // Absolute difference element per element
        std::transform(x.begin(), x.end(), xn.begin(), abs_error.begin(),
                       [](double a, double b) { return std::abs(a - b); });
        // Maximum among the absolute differences
        const auto &criterion = *(std::max_element(abs_error.begin(), abs_error.end()));
        // Checking if criterion is satisfied
        bool is_cvg = (criterion < given_criterion);
        abs_error.clear();
        return std::make_tuple(is_cvg, criterion);
      });
}
