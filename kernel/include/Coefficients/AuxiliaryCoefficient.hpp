/**
 * @file AuxiliaryCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Coefficient depending only on auxiliary variables, built by an integrator
 * @version 0.1
 * @date 2026-11-09
 *
 * Copyright CEA (C) 2026
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

#include <functional>
#include <span>
#include <string>
#include <vector>

#include "Coefficients/EvalAtPoints.hpp"
#include "Coefficients/FunctionCoefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Coefficient of one variable whose value depends only on auxiliary variables.
 *
 * Built by an integrator from its parameters (not given in the JSON file of the coefficients).
 * Its derivatives with respect to the unknown are zero. EVAL is an evaluation structure, as
 * those written by GenerateCoefficient.py (host and device), with in addition:
 *  - max_auxiliary_index(): largest index of the auxiliary variables read by the structure
 *    (indices of the auxiliary variables of the integrator, known at run time).
 *
 * @tparam EVAL Evaluation structure
 */
template <class EVAL>
class AuxiliaryCoefficient : public FunctionCoefficient {
 private:
  EVAL eval_;

  /// Checks that the auxiliary values hold the variables read by the structure
  void check_auxiliary_values(const std::size_t nb_auxiliary_values) const {
    MFEM_VERIFY(static_cast<int>(nb_auxiliary_values) > this->eval_.max_auxiliary_index(),
                "AuxiliaryCoefficient: auxiliary variable " +
                    std::to_string(this->eval_.max_auxiliary_index()) + " required, but only " +
                    std::to_string(nb_auxiliary_values) + " given");
  }

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final {
    return [this](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
      this->check_auxiliary_values(auxiliary_vector.size());
      double value = 0.0;
      this->eval_(input_vector.data(), nullptr, auxiliary_vector.data(), nullptr, this->time_,
                  &value, nullptr, nullptr);
      return value;
    };
  }

  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final {
    return [](const std::span<const double>&, const std::span<const double>&,
              const std::span<const double>&,
              const unsigned int) { return std::vector<double>(EVAL::nb_variables, 0.0); };
  }

  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final {
    return [](const std::span<const double>&, const std::span<const double>&,
              const std::span<const double>&, const unsigned int) {
      return std::vector<double>(EVAL::nb_variables * EVAL::nb_variables, 0.0);
    };
  }

 public:
  explicit AuxiliaryCoefficient(const EVAL& eval) : eval_(eval) {}
  virtual ~AuxiliaryCoefficient() = default;

  const EVAL& get_eval() const { return eval_; }

  bool eval_batch(const int npts, const mfem::Vector& values, const mfem::Vector& explicit_values,
                  const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                  mfem::Vector* F, mfem::Vector* gradient, mfem::Vector* hessian) const final {
    this->check_auxiliary_values(npts > 0 ? static_cast<std::size_t>(auxiliary_values.Size() / npts)
                                          : 0u);
    eval_at_points(eval_, npts, this->time_, values, explicit_values, auxiliary_values, coordinates,
                   F, gradient, hessian);
    return true;
  }
  int batch_nb_variables() const final { return EVAL::nb_variables; }
  bool is_single_variable() const final { return EVAL::nb_variables == 1; }
};
