/**
 * @file SteadyOperator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Steady  Operator
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

#include <memory>
#include <string>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Operators/OperatorBase.hpp"
#include "Operators/SteadyReducedOperator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief SteadyOperator class
 *
 */
template <class T, int DIM>
class SteadyOperator : public OperatorBase<T, DIM> {
 private:
  SteadyPhaseFieldReducedOperator* steady_reduced_oper;
  void free_memory();

 public:
  SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                 const std::vector<std::string>& integrators);

  SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                 const std::vector<std::string>& integrators,
                 std::vector<AnalyticalFunctions<DIM>> source_term_name);

  SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                 const std::vector<std::string>& integrators, const Parameters& params);

  SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                 const std::vector<std::string>& integrators, const Parameters& params,
                 std::vector<AnalyticalFunctions<DIM>> source_term_name);

  void Mult([[maybe_unused]] const mfem::Vector& k, [[maybe_unused]] mfem::Vector& y)
      const override {  // Nothing to be done because of manage by steadyreducedoperator
  }
  // Virtual methods

  void initialize(const double& initial_time, Variables<T, DIM>& vars,
                  std::vector<Variables<T, DIM>*> auxvars) override;
  void SetTransientParameters(const std::vector<mfem::Vector>& u_vet) override;
  void solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
             const double& current_time, double dt, const int iter) override;
};

#include "Operators/SteadyOperator.tpp"
