/**
 * @file CalphadBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for Calphad objets
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * @anchor CalphadBase
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
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Calphad/CalphadUtils.hpp"
#include "Coefficients/CommonCoefficients.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

// Previous declaration of the KKS class
template <typename T>
class KKS;

template <typename T>
class CalphadBase {
 private:
  bool is_KKS_;

 protected:
  // Smart pointer used for KKS studies
  std::shared_ptr<KKS<T>> KKS_;

  // Common methods for CALPHAD studies
  std::shared_ptr<CalphadUtils<T>> CU_;

  //  Heat capacity. Nodal values. key: [node]
  std::map<int, double> heat_capacity_;

  // Mobility for each element in a given phase.
  // Nodal values. key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> mobilities_;

  // Site fraction of a given constituant, in a given sublattice for each phase.
  // Nodal values. key: [node, phase, cons, sub]
  std::map<std::tuple<int, std::string, std::string, int>, double> site_fraction_;

  void clear_containers();

  void update_outputs(
      const int dt, const size_t nb_nodes,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system,
      const std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>&
          previous_output_system);

 public:
  // Parameters for CALPHAD problems
  const Parameters& params_;

  // Symbol of the chemical element removed from the system when initializing equilibrium
  // calculations (performed with molar fractions).
  std::string element_removed_from_ic_;

  // Containers used to store the results of the equilibrium calculations
  std::map<int, int> error_equilibrium_;

  // Chemical potential. Nodal values. key: [node, elem]
  std::map<std::tuple<int, std::string>, double> chemical_potentials_;
  std::map<std::tuple<int, std::string>, double> diffusion_chemical_potentials_;

  // Mole fraction of phase. Nodal values. key: [node, phase]
  std::map<std::tuple<int, std::string>, double> mole_fraction_of_phase_;

  // Element mole fraction by phase. Nodal values. key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> elem_mole_fraction_by_phase_;

  // Energy for each phase. Nodal values. key: [node, phase, symbol]
  std::map<std::tuple<int, std::string, std::string>, double> energies_of_phases_;

  // Driving force for each phase. Nodal values. key: [node, phase]
  std::map<std::tuple<int, std::string>, double> driving_forces_;
  std::map<std::tuple<int, std::string>, double> nucleus_;

  /// Time integral results
  std::multimap<IterationKey, SpecializedValue> time_specialized_;
  const std::multimap<IterationKey, SpecializedValue> get_time_specialized() const;
  void clear_time_specialized();

  explicit CalphadBase(const Parameters& params);
  CalphadBase(const Parameters& params, bool is_KKS);

  virtual void get_parameters();

  virtual void initialize(
      const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) = 0;

  void global_execute(
      const int dt, const double time_step, const std::vector<T>& tp_gf,
      const std::vector<std::tuple<std::string, std::string>>& chemicalsystem,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system,
      const std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>& previous_output_system,
      std::optional<const std::tuple<std::string, T, T>> phase_fields = std::nullopt,
      std::optional<const std::vector<T>> tp_gf_old = std::nullopt,
      std::optional<const std::vector<std::tuple<std::string, std::string, T, T>>> x_gf =
          std::nullopt,
      std::optional<const std::vector<std::tuple<std::string, T>>> coordinates = std::nullopt);

  virtual void execute(const int dt, const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
                       const std::vector<std::tuple<std::string, std::string>>& chemicalsystem,
                       std::optional<std::vector<std::tuple<std::string, std::string, double>>>
                           status_phase = std::nullopt) = 0;

  virtual void finalize() = 0;

  virtual ~CalphadBase();
};

#include "Calphad/CalphadBase.tpp"
