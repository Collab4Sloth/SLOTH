/**
 * @file VariableSolution.tpp
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

////////////////////////////////
// VariableSolution
////////////////////////////////

/**
 * @brief Construct a VariableSolution from a constant value
 *
 * @tparam DIM Spatial dimension.
 * @param value
 */
template <int DIM>
VariableSolution<DIM>::VariableSolution(double value) : data_(value) {}

/**
 * @brief Construct a VariableSolution from an AnalyticalFunctions
 *
 * @tparam DIM Spatial dimension.
 * @param function
 */
template <int DIM>
VariableSolution<DIM>::VariableSolution(const AnalyticalFunctions<DIM>& function)
    : data_(function) {}

/**
 * @brief Construct a VariableSolution from an mfem::FunctionCoefficient
 *
 * @tparam DIM Spatial dimension.
 * @param function
 */
template <int DIM>
VariableSolution<DIM>::VariableSolution(const mfem::FunctionCoefficient& function)
    : data_(function) {}

/**
 * @brief Construct a VariableSolution from a time-independent std::function
 *        (same as passing an mfem::FunctionCoefficient built from it)
 *
 * @tparam DIM Spatial dimension.
 * @param function
 */
template <int DIM>
VariableSolution<DIM>::VariableSolution(std::function<double(const mfem::Vector&)> function)
    : data_(mfem::FunctionCoefficient(std::move(function))) {}

/**
 * @brief Construct a VariableSolution from a time-dependent std::function
 *        (same as passing an mfem::FunctionCoefficient built from it)
 *
 * @tparam DIM Spatial dimension.
 * @param function
 */
template <int DIM>
VariableSolution<DIM>::VariableSolution(std::function<double(const mfem::Vector&, double)> function)
    : data_(mfem::FunctionCoefficient(std::move(function))) {}

/**
 * @brief Construct a VariableSolution from a SLOTH ::Coefficient
 *
 * @tparam DIM Spatial dimension.
 * @param coefficient
 */
template <int DIM>
VariableSolution<DIM>::VariableSolution(const ::Coefficient& coefficient) : data_(coefficient) {}

/**
 * @brief True if a solution has been given (false for a default-constructed one)
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
bool VariableSolution<DIM>::is_defined() const {
  return !std::holds_alternative<std::monostate>(this->data_);
}

/**
 * @brief True if the solution is given by a SLOTH ::Coefficient
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
bool VariableSolution<DIM>::is_sloth_coefficient() const {
  return std::holds_alternative<::Coefficient>(this->data_);
}

/**
 * @brief Convert the solution into an mfem::Coefficient
 *
 * @details A SLOTH ::Coefficient is wrapped into an MfemFromSlothCoefficient.
 *
 * @tparam DIM Spatial dimension.
 * @return std::shared_ptr<mfem::Coefficient> nullptr if no solution has been given
 */
template <int DIM>
std::shared_ptr<mfem::Coefficient> VariableSolution<DIM>::SlothToMfemCoefficient() const {
  return std::visit(
      [](const auto& data) -> std::shared_ptr<mfem::Coefficient> {
        using D = std::decay_t<decltype(data)>;
        if constexpr (std::is_same_v<D, std::monostate>) {
          return nullptr;
        } else if constexpr (std::is_same_v<D, double>) {
          return std::make_shared<mfem::ConstantCoefficient>(data);
        } else if constexpr (std::is_same_v<D, AnalyticalFunctions<DIM>>) {
          return std::make_shared<mfem::FunctionCoefficient>(data.getFunction());
        } else if constexpr (std::is_same_v<D, mfem::FunctionCoefficient>) {
          return std::make_shared<mfem::FunctionCoefficient>(data);
        } else {
          return std::make_shared<MfemFromSlothCoefficient>(data);
        }
      },
      this->data_);
}

/**
 * @brief Return the solution as a std::function, when it is given by an AnalyticalFunctions
 *
 * @tparam DIM Spatial dimension.
 * @return std::shared_ptr<std::function<double(const mfem::Vector&, double)>> nullptr otherwise
 */
template <int DIM>
std::shared_ptr<std::function<double(const mfem::Vector&, double)>>
VariableSolution<DIM>::get_function() const {
  if (const auto* function = std::get_if<AnalyticalFunctions<DIM>>(&this->data_)) {
    return std::make_shared<std::function<double(const mfem::Vector&, double)>>(
        function->getFunction());
  }
  return nullptr;
}

////////////////////////////////
// InitialCondition
////////////////////////////////

/**
 * @brief Construct an InitialCondition read from a MFEM grid function file
 *
 * @tparam DIM Spatial dimension.
 * @param input_file_dir (file name, directory). The rank number is appended to the file name.
 */
template <int DIM>
InitialCondition<DIM>::InitialCondition(const std::tuple<std::string, std::string>& input_file_dir)
    : input_file_dir_(input_file_dir) {}

/**
 * @brief Construct an InitialCondition read from a MFEM grid function file
 *
 * @tparam DIM Spatial dimension.
 * @tparam A Type convertible to std::string
 * @tparam B Type convertible to std::string
 * @param input_file_dir (file name, directory). The rank number is appended to the file name.
 */
template <int DIM>
template <class A, class B>
  requires std::constructible_from<std::string, const A&> &&
           std::constructible_from<std::string, const B&>
InitialCondition<DIM>::InitialCondition(const std::tuple<A, B>& input_file_dir)
    : input_file_dir_(std::make_tuple(std::string(std::get<0>(input_file_dir)),
                                      std::string(std::get<1>(input_file_dir)))) {}

/**
 * @brief Construct an InitialCondition read from a MFEM grid function file
 *
 * @tparam DIM Spatial dimension.
 * @param input_file File name. The rank number is appended to it.
 * @param input_dir Directory
 */
template <int DIM>
InitialCondition<DIM>::InitialCondition(const std::string& input_file, const std::string& input_dir)
    : input_file_dir_(std::make_tuple(input_file, input_dir)) {}

/**
 * @brief True if the initial condition is read from a MFEM grid function file
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
bool InitialCondition<DIM>::is_input_file() const {
  return this->input_file_dir_.has_value();
}

/**
 * @brief Return the (file name, directory) of an initial condition read from file
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
const std::tuple<std::string, std::string>& InitialCondition<DIM>::get_input_file_dir() const {
  MFEM_VERIFY(this->is_input_file(),
              "InitialCondition::get_input_file_dir(): the initial condition is not read from a "
              "file.");
  return *this->input_file_dir_;
}

/**
 * @brief True if the initial condition is given by a SLOTH ::Coefficient that is neither a
 *        constant nor file-based, i.e. that may depend on other variables
 *
 * @tparam DIM Spatial dimension.
 */
template <int DIM>
bool InitialCondition<DIM>::requires_dependencies() const {
  if (const auto* coefficient = std::get_if<::Coefficient>(&this->data_)) {
    return !coefficient->is_scalar() && !coefficient->is_file();
  }
  return false;
}
