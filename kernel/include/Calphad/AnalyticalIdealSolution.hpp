/**
 * @file AnalyticalIdealSolution.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Analytical thermodynamic description for an ideal solution
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
#pragma once
#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

template <typename T>
class AnalyticalIdealSolution : public CalphadBase<T> {
 private:
  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  constexpr explicit AnalyticalIdealSolution(const Parameters& params);
  constexpr AnalyticalIdealSolution(const Parameters& params, bool is_KKS);

  void initialize([[maybe_unused]] const std::vector<std::tuple<std::string, std::string>>&
                      sorted_chemical_system) override {}

  void execute(const int dt, const std::set<int>& list_nodes, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::optional<std::vector<std::tuple<std::string, std::string, double>>>
                   status_phase = std::nullopt) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~AnalyticalIdealSolution();
};

#include "Calphad/AnalyticalIdealSolution.tpp"
