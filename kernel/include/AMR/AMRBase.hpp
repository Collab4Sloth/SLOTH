/**
 * @file AMRBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class orchestrating adaptive mesh refinement (AMR): defines
 *        the common criteria (error estimator integrator, thresholds,
 *        NC limit) and the shared refine/derefine/update orchestration
 *        (InitialRefine, StepRefine, StepDerefine), delegating the actual
 *        per-strategy error estimation and mesh mutation to derived
 *        classes via the Refine()/Derefine() pure virtual methods.
 * @version 0.1
 * @date 2026-08-04
 *
 * @copyright CEA (C) 2026
 *
 * @anchor AMR
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
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "AMR/SlothErrorEstimators.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class VAR>
class AMRBase {
 public:
  using Variable_t = std::remove_reference_t<decltype(std::declval<VAR&>()[0])>;

 protected:
  mfem::ParMesh& par_mesh_;
  bool is_nc_simplices_;
  double amr_max_elem_error_{0.0};
  double amr_scale_down_factor_{0.25};
  int amr_nc_limit_{0};
  int amr_max_level_{0};
  int amr_max_preref_cycles_{1};
  bool amr_rebalance_{false};
  // Auxiliary variables of the Problem owning this AMR object. Those which are pilot
  // variables of another AMR object on the same mesh take part in the derefinement decision
  std::vector<VAR*> coupled_vars_;
  bool pending_rebalance_{false};

  SlothErrorEstimators* amr_estimator_{nullptr};

  virtual bool Refine(VAR& vars) = 0;
  virtual bool Derefine(VAR& var, std::vector<VAR*> auxvars) = 0;
  // Variables of the owning Problem associated with the refinement criterion of this object
  virtual std::vector<Variable_t*> PilotVariables(VAR& vars) = 0;

 public:
  AMRBase(mfem::ParMesh& mesh, bool is_nc_simplices);
  virtual ~AMRBase();
  AMRBase(const AMRBase&) = delete;
  AMRBase& operator=(const AMRBase&) = delete;

  void SetCriteria(SlothErrorEstimators* estimator, const Parameters& params);

  void EnsureCriteriaSet() const;

  void InitialRefine(VAR& vars, std::vector<VAR*> auxvars);
  void StepRefine(VAR& vars, std::vector<VAR*> auxvars);
  void StepDerefine(VAR& vars, std::vector<VAR*> auxvars);
  void StepRebalance(VAR& vars, std::vector<VAR*> auxvars);
  void SetCoupledVariables(const std::vector<VAR*>& coupled_vars);

  mfem::Array<mfem::Refinement> FilterRefinedCandidates(
      const mfem::Array<mfem::Refinement>& candidates);

 protected:
  void EnsureNCMeshIfNeeded();
  bool IsOnMesh(const mfem::ParFiniteElementSpace* fes) const;
  void UpdateVariables(VAR& vars, const std::vector<VAR*>& auxvars);
  bool Rebalance(VAR& vars, const std::vector<VAR*>& auxvars);
  void PrintLoadBalance(const std::string& stage) const;

  // Pilot variables of every AMR object, per mesh: variable -> AMR object it drives
  using PilotRegistry_t = std::map<const mfem::ParMesh*, std::map<const void*, AMRBase<VAR>*>>;
  static PilotRegistry_t& PilotRegistry();
  void RegisterPilots(VAR& vars);
  void AccumulateNormalizedError(Variable_t& var, mfem::Vector& combined_error);
  bool DerefineByConsensus(VAR& vars, const std::vector<VAR*>& coupled_vars);
  static Variable_t* FindVariable(VAR& vars, const std::string& name);
};

#include "AMR/AMRBase.tpp"
