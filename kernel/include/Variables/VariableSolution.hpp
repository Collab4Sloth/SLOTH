/**
 * @file VariableSolution.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Initial condition and analytical solution of a Variable, whatever the way they are given
 * @version 0.1
 * @date 2026-11-26
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
#include <concepts>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <variant>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "Coefficients/Coefficient.hpp"
#include "Coefficients/MfemFromSlothCoefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Common description of a field attached to a Variable (initial condition or
 *        analytical solution), whatever the way it is given by the user.
 *
 * @details All the constructors are implicit on purpose: they let the Variable
 *          constructors accept a double, an AnalyticalFunctions, an
 *          mfem::FunctionCoefficient or a SLOTH ::Coefficient exactly as before.
 *          SlothToMfemCoefficient() converts any of them into an mfem::Coefficient
 *          that can be projected onto a grid function or used to compute errors.
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
class VariableSolution {
 protected:
  using SolutionData = std::variant<std::monostate, double, AnalyticalFunctions<DIM>,
                                    mfem::FunctionCoefficient, ::Coefficient>;
  SolutionData data_;

  VariableSolution() = default;

 public:
  VariableSolution(double value);                                         // NOLINT [implicit]
  VariableSolution(const AnalyticalFunctions<DIM>& function);             // NOLINT [implicit]
  VariableSolution(const mfem::FunctionCoefficient& function);            // NOLINT [implicit]
  VariableSolution(std::function<double(const mfem::Vector&)> function);  // NOLINT [implicit]
  VariableSolution(                                                       // NOLINT [implicit]
      std::function<double(const mfem::Vector&, double)> function);
  VariableSolution(const ::Coefficient& coefficient);  // NOLINT [implicit]

  bool is_defined() const;
  bool is_sloth_coefficient() const;

  std::shared_ptr<mfem::Coefficient> SlothToMfemCoefficient() const;
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> get_function() const;

  virtual ~VariableSolution() = default;
};

/**
 * @brief Initial condition of a Variable.
 *
 * @details In addition to the cases handled by VariableSolution, an initial
 *          condition can be read from a MFEM grid function file (one file per
 *          MPI rank). Such an initial condition cannot be reprojected after a
 *          mesh change.
 *
 *          An initial condition given by a SLOTH ::Coefficient that is neither a
 *          constant nor file-based may depend on other variables (see
 *          Variable::set_dependencies()).
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
class InitialCondition : public VariableSolution<DIM> {
 private:
  std::optional<std::tuple<std::string, std::string>> input_file_dir_;

 public:
  using VariableSolution<DIM>::VariableSolution;

  InitialCondition(const std::tuple<std::string, std::string>& input_file_dir);  // NOLINT
  template <class A, class B>
    requires std::constructible_from<std::string, const A&> &&
             std::constructible_from<std::string, const B&>
  InitialCondition(const std::tuple<A, B>& input_file_dir);  // NOLINT [implicit]
  InitialCondition(const std::string& input_file, const std::string& input_dir);

  bool is_input_file() const;
  const std::tuple<std::string, std::string>& get_input_file_dir() const;
  bool requires_dependencies() const;
};

/**
 * @brief Analytical solution of a Variable, used to compute errors.
 *
 * @details Default-constructed, it means "no analytical solution". It never
 *          depends on other variables: a SLOTH ::Coefficient that is neither a
 *          constant nor file-based is evaluated with an empty list of values (and
 *          the physical point, if it uses coordinates).
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
class AnalyticalSolution : public VariableSolution<DIM> {
 public:
  using VariableSolution<DIM>::VariableSolution;
  AnalyticalSolution() = default;
};

#include "Variables/VariableSolution.tpp"
