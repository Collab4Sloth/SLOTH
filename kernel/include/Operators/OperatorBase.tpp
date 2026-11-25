/**
 * @file OperatorBase.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for building Steady and TimeDependent PhaseFieldOperators
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * @anchor OperatorBase
 *
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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/Coefficients.hpp"
#include "Integrators/ListIntegrators.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Build and initialize an integrator of the right-hand side
 *
 * @details A built-in integrator is built from its name; a user integrator by its own builder,
 *          with the same arguments.
 *
 * @tparam T
 * @tparam DIM
 * @param dt Time step
 * @param nlfi Name of a built-in integrator, or user integrator
 * @param u Unknowns
 * @return NLFI*
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* OperatorBase<T, DIM>::set_nlfi_ptr(
    const double dt, const IntegratorEntry<Variables<T, DIM>>& nlfi,
    const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("OperatorBase::set_nlfi_ptr");
  std::vector<mfem::ParGridFunction> vun;
  std::vector<mfem::ParGridFunction> vauxn;
  for (unsigned int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(std::move(un));
  }
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (auto auxvars : auxvar_vec->getVariables()) {
      auto fes = auxvars.get_fespace();
      mfem::ParGridFunction auxn(fes);
      auto auxvar_n = auxvars.get_second_to_last();
      auxn.SetFromTrueDofs(auxvar_n);
      vauxn.emplace_back(std::move(auxn));
    }
  }

  const Parameters& all_params = this->params_ - this->default_p_;
  SlothNLFormIntegrator<Variables<T, DIM>>* rhs_nlfi = nullptr;
  if (const auto* user = std::get_if<UserIntegrator<Variables<T, DIM>>>(&nlfi)) {
    const UserIntegratorArgs<Variables<T, DIM>> args{
        this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_};
    rhs_nlfi = user->build(args);
  } else {
    rhs_nlfi = this->get_rhs_integrator(dt, std::get<std::string>(nlfi), vun, vauxn, all_params);
  }
  const std::string& verify_msg = "OperatorBase: integrator '" +
                                  integrator_name<Variables<T, DIM>>(nlfi) +
                                  "' could not be built.";
  MFEM_VERIFY(rhs_nlfi != nullptr, verify_msg.c_str());
  rhs_nlfi->set_physical_time_step(this->current_dt_ > 0. ? this->current_dt_ : dt);
  rhs_nlfi->init();
  return rhs_nlfi;
}

/**
 * @brief Set the NonLinearFormIntegrator dedicated to boundary terms
 *
 * @tparam T
 * @tparam DIM
 * @tparam OPEBASE
 * @param nlfi Name of integrator
 * @param u Unknown
 * @param block Block to which the boundary term is applied
 * @param bdr_id Id of the boundary to which the boundary term is applied
 * @return NLFI*
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* OperatorBase<T, DIM>::set_bdr_nlfi_ptr(
    const double dt, const std::string nlfi, const std::vector<mfem::Vector>& u,
    const unsigned int block, const unsigned int bdr_id) {
  Catch_Time_Section("OperatorBase::set_bdr_nlfi_ptr");

  std::vector<mfem::ParGridFunction> vun;
  std::vector<mfem::ParGridFunction> vauxn;
  for (unsigned int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(std::move(un));
  }
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (auto auxvars : auxvar_vec->getVariables()) {
      auto fes = auxvars.get_fespace();
      mfem::ParGridFunction auxn(fes);
      auto auxvar_n = auxvars.get_second_to_last();
      auxn.SetFromTrueDofs(auxvar_n);
      vauxn.emplace_back(std::move(auxn));
    }
  }

  const Parameters& all_params = this->params_ - this->default_p_;
  auto bdr_nlfi = this->get_bdr_integrator(dt, nlfi, vun, vauxn, all_params, block, bdr_id);
  bdr_nlfi->set_physical_time_step(this->current_dt_ > 0. ? this->current_dt_ : dt);
  bdr_nlfi->init();
  return bdr_nlfi;
}

/**
 * @brief Return the total height (output=rows of Operator) of the PDE system
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatials
 * @return int
 */
template <class T, int DIM>
int OperatorBase<T, DIM>::compute_total_height(
    const std::vector<SpatialDiscretization<T, DIM>*>& spatials) {
  int total_size = 0;
  for (const auto* s : spatials) {
    total_size += s->getSize();
  }
  return total_size;
}

/**
 * @brief Return the total width (input=column of Operator) of the PDE system
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatials
 * @return int
 */
template <class T, int DIM>
int OperatorBase<T, DIM>::compute_total_width(
    const std::vector<SpatialDiscretization<T, DIM>*>& spatials) {
  int total_size = 0;
  for (const auto* s : spatials) {
    total_size += s->getSize();
  }
  return total_size;
}

/**
 * @brief Set the Coefficients associated with the current Operator
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param coefficients
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::set_coefficients(const std::vector<Coefficients>& coefficients) {
  this->coefficients_ = coefficients;
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatial
 */
template <class T, int DIM>
OperatorBase<T, DIM>::OperatorBase(
    const std::vector<IntegratorEntry<Variables<T, DIM>>>& integrators,
    std::vector<SpatialDiscretization<T, DIM>*> spatials)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(default_params_),
      RHS(NULL),
      rhs_solver_(nullptr),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      rhs_integrators_(integrators) {
  this->fes_.SetSize(spatials.size());
  this->bcs_.reserve(spatials.size());
  this->ess_tdof_list_.reserve(spatials.size());
  this->block_trueOffsets_.SetSize(spatials.size() + 1);

  this->block_trueOffsets_[0] = 0;
  int i = 0;
  for (const auto* s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->set_default_solver();
  // Sizes above are those of the current mesh: a later mesh change (e.g. pre-refinement by the
  // AMR of another Problem, before initialize()) is detected by mesh_has_changed()
  this->record_fes_sequence();
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatial
 * @param params
 */
template <class T, int DIM>
OperatorBase<T, DIM>::OperatorBase(
    const std::vector<IntegratorEntry<Variables<T, DIM>>>& integrators,
    std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(params),
      RHS(NULL),
      rhs_solver_(nullptr),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      rhs_integrators_(integrators) {
  this->fes_.SetSize(spatials.size());
  this->bcs_.reserve(spatials.size());
  this->ess_tdof_list_.reserve(spatials.size());

  this->block_trueOffsets_.SetSize(spatials.size() + 1);

  this->block_trueOffsets_[0] = 0;
  int i = 0;
  for (const auto* s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->set_default_solver();
  // Sizes above are those of the current mesh: a later mesh change (e.g. pre-refinement by the
  // AMR of another Problem, before initialize()) is detected by mesh_has_changed()
  this->record_fes_sequence();
}

/**
 * @brief  Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param initial_time
 * @param vars
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::initialize([[maybe_unused]] const double& initial_time, const double dt,
                                      Variables<T, DIM>& vars,
                                      std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("OperatorBase::initialize");

  this->auxvariables_ = auxvars;
  const auto nvars = vars.get_variables_number();
  std::vector<mfem::Vector> u_vect;
  u_vect.reserve(nvars);
  const auto auxvars_unk = this->get_auxvars_unknowns();
  for (unsigned int iv = 0; iv < nvars; iv++) {
    auto& vv = vars[iv];
    auto u = vv.get_unknown();

    this->bcs_.emplace_back(vv.get_boundary_conditions());
    // Boundaries named in the coefficients: indices given by the boundary conditions
    this->bcs_[iv]->resolve_coefficient_boundaries(this->coefficients_[iv]);
    // The mesh may have changed since the boundary conditions were built (e.g. pre-refinement
    // by the AMR of a Problem initialized before this one): essential dofs on the current mesh
    this->bcs_[iv]->UpdateEssentialDofs();
    this->ess_tdof_list_.emplace_back(this->bcs_[iv]->GetEssentialDofs());

    this->bcs_[iv]->SetDirichletBoundaryConditions(u, this->coefficients_[iv], auxvars_unk);
    vv.update(u);
    u_vect.emplace_back(u);
  }
  // Sizes of the operator on the current mesh (built on a coarser mesh if the mesh was refined
  // between the construction and the initialization, e.g. by the AMR of another Problem)
  this->sync_with_mesh();
  this->SetTransientParameters(dt, u_vect);
}

/**
 * @brief Store the sequence number of each finite element space.
 *
 * @details mfem increments FiniteElementSpace::GetSequence() at each Update()
 *          (refinement/derefinement). Comparing it with the stored value tells whether
 *          the operator must be resynchronized (see mesh_has_changed()).
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::record_fes_sequence() {
  this->fes_sequence_.resize(static_cast<size_t>(this->fes_.Size()));
  for (int i = 0; i < this->fes_.Size(); ++i) {
    this->fes_sequence_[static_cast<size_t>(i)] = this->fes_[i]->GetSequence();
  }
}

/**
 * @brief True if at least one finite element space has been updated since the last call to
 *        record_fes_sequence().
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
bool OperatorBase<T, DIM>::mesh_has_changed() const {
  if (this->fes_sequence_.size() != static_cast<size_t>(this->fes_.Size())) return true;
  for (int i = 0; i < this->fes_.Size(); ++i) {
    if (this->fes_[i]->GetSequence() != this->fes_sequence_[static_cast<size_t>(i)]) return true;
  }
  return false;
}

/**
 * @brief Resynchronize the operator with the current mesh if it has changed since the last
 *        synchronization (AMR of this Problem or of any Problem sharing the mesh).
 *
 * @details Updates the sizes, the block offsets and the essential dofs
 *          (UpdateAfterMeshChange()), then records the new state of the finite element spaces.
 *          Derived operators override this method to reset their own mesh-dependent data, after
 *          calling it.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @return true The mesh had changed: the operator has been resynchronized.
 * @return false The operator was already up to date.
 */
template <class T, int DIM>
bool OperatorBase<T, DIM>::sync_with_mesh() {
  if (!this->mesh_has_changed()) return false;
  this->UpdateAfterMeshChange();
  this->record_fes_sequence();
  return true;
}

/**
 * @brief Current unknown vectors of all the auxiliary variables, in the order of
 *        `auxvariables_` (input of the Dirichlet coefficients).
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @return std::vector<mfem::Vector> One vector of true dofs per auxiliary variable.
 */
template <class T, int DIM>
std::vector<mfem::Vector> OperatorBase<T, DIM>::get_auxvars_unknowns() const {
  std::vector<mfem::Vector> auxvars_unk;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      auxvars_unk.emplace_back(auxvar.get_unknown());
    }
  }
  return auxvars_unk;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Build the NonLinear Form Integrator associated with the RHS of the PDEs
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param dt
 * @param u
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::build_rhs_nonlinear_form(const double dt,
                                                    const std::vector<mfem::Vector>& u_vect) {
  // Catch_Time_Section("OperatorBase::build_rhs_nonlinear_form");
  if (this->RHS != nullptr) {
    delete this->RHS;
  }
  this->RHS = new SlothParBlockNonlinearForm<Variables<T, DIM>>(this->fes_);
  for (const auto& integrator : this->rhs_integrators_) {
    auto integrator_ptr = this->set_nlfi_ptr(dt, integrator, u_vect);
    this->RHS->AddDomainIntegrator(integrator_ptr);
  }

  // Add boundary integrators
  const int fes_size = this->block_trueOffsets_.Size() - 1;
  mfem::Array<int> Robin_bdr;
  mfem::Array<int> Neumann_bdr;

  auto make_boundary_marker = [](int nb_bdr, int j) {
    mfem::Array<int> marker(nb_bdr);
    marker = 0;
    marker[j] = 1;
    return marker;
  };

  // Advection: an Inflow coefficient attached to a boundary
  // activates the inflow/outflow term
  const auto has_rhs_integrator = [this](const std::string& name) {
    return std::any_of(
        this->rhs_integrators_.begin(), this->rhs_integrators_.end(),
        [&name](const auto& entry) { return integrator_name<Variables<T, DIM>>(entry) == name; });
  };
  const bool advection = has_rhs_integrator("Advection");
  const bool split_advection = has_rhs_integrator("SplitAdvection");

  // Loop over variables
  for (int i = 0; i < fes_size; ++i) {
    Robin_bdr = this->bcs_[i]->get_marker_array("Robin");
    const int nb_bdr = Robin_bdr.Size();
    Neumann_bdr = this->bcs_[i]->get_marker_array("Neumann");

    // All finite element spaces are assumed to live on the same mesh.
    // Therefore all boundary marker arrays have the same size.
    if (this->array_bdr_.empty()) {
      this->array_bdr_.resize(nb_bdr);
    }
    MFEM_VERIFY(this->array_bdr_.size() == nb_bdr, "Inconsistent number of boundary attributes.");

    const Coefficients& coefficients = this->coefficients_[i];
    const int coef_size = coefficients.size();

    // Loop over boundaries
    for (int j = 0; j < nb_bdr; j++) {
      // Add Neumann integrator
      if (Neumann_bdr[j] > 0) {
        // Check if a coefficient is given for this bc
        // (else Homogeneous Neumann)
        bool has_neumann_coeff = false;
        for (int l = 0; l < coef_size; l++) {
          const auto& coef = coefficients[l];
          if (coef.get_type() == GlossaryType::Neumann) {
            auto bdr_ids = coef.get_bdr_index_coef();
            if (std::find(bdr_ids.begin(), bdr_ids.end(), j) != bdr_ids.end()) {
              has_neumann_coeff = true;
              break;
            }
          }
        }
        // If a Neumann coefficient is specified,
        // add the Neumann integrator even if the coefficient is zero
        if (has_neumann_coeff) {
          this->array_bdr_[j] = make_boundary_marker(nb_bdr, j);
          auto integrator_ptr = this->set_bdr_nlfi_ptr(dt, "Neumann", u_vect, i, j);
          this->RHS->AddBoundaryIntegrator(integrator_ptr, this->array_bdr_[j]);
        }
      }
      // Add Robin integrator
      if (Robin_bdr[j] > 0) {
        this->array_bdr_[j] = make_boundary_marker(nb_bdr, j);
        auto integrator_ptr = this->set_bdr_nlfi_ptr(dt, "Robin", u_vect, i, j);
        this->RHS->AddBoundaryIntegrator(integrator_ptr, this->array_bdr_[j]);
      }

      // Add advection inflow/outflow integrator
      // (independent of the boundary condition type)
      bool has_inflow_coeff = false;
      for (int l = 0; l < coef_size; l++) {
        const auto& coef = coefficients[l];
        if (coef.get_type() == GlossaryType::Inflow) {
          auto bdr_ids = coef.get_bdr_index_coef();
          if (std::find(bdr_ids.begin(), bdr_ids.end(), j) != bdr_ids.end()) {
            has_inflow_coeff = true;
            break;
          }
        }
      }
      if (has_inflow_coeff) {
        MFEM_VERIFY(advection || split_advection,
                    "An Inflow coefficient is given but neither Advection nor SplitAdvection is "
                    "used. Please check your data.");
        this->array_bdr_[j] = make_boundary_marker(nb_bdr, j);
        const std::string name = split_advection ? "AdvectionBoundaryCH" : "AdvectionBoundary";
        auto integrator_ptr = this->set_bdr_nlfi_ptr(dt, name, u_vect, i, j);
        this->RHS->AddBoundaryIntegrator(integrator_ptr, this->array_bdr_[j]);
      }
      //
    }
  }
}

/**
 * @brief Configure the Newton solver
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param oper
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::SetNewtonAlgorithm(mfem::Operator* oper) {
  this->rhs_solver_ =
      new NLSolver(this->nl_solver_, this->nl_solver_params_, this->solver_, this->solver_params_,
                   this->precond_, this->precond_params_, *oper);
  this->newton_solver_.reset();
  this->newton_solver_ = this->rhs_solver_->get_nl_solver();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compute L2 error
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param it
 * @param t
 * @param dt
 * @param id_var
 * @param name
 * @param u
 * @param solution_func
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::ComputeError(
    const int& it, const double& t, const double& dt, const int id_var, const std::string& name,
    const mfem::Vector& u, std::function<double(const mfem::Vector&, double)> solution_func) {
  Catch_Time_Section("OperatorBase::ComputeError");

  mfem::ParGridFunction gf(this->fes_[id_var]);
  mfem::ParGridFunction zero(this->fes_[id_var]);
  zero = 0.0;

  gf.SetFromTrueDofs(u);
  mfem::FunctionCoefficient solution_coef(solution_func);
  solution_coef.SetTime(t);

  const auto errorL2 = gf.ComputeLpError(2., solution_coef);
  const auto errorLinf = gf.ComputeLpError(mfem::infinity(), solution_coef);
  const auto norm_solution = zero.ComputeLpError(2, solution_coef);
  const auto normalized_error = errorL2 / norm_solution;

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_L2-error[-]", errorL2));
  this->time_specialized_.emplace(
      IterationKey(it, dt, t),
      SpecializedValue(name + "_L2-error normalized[-]", normalized_error));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_Linf-error [-]", errorLinf));
}

/**
 * @brief Compute integrals of a variable over the domain
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param it
 * @param t
 * @param dt
 * @param id_var
 * @param name
 * @param u
 * @param integral_threshold
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::ComputeIntegral(const int& it, const double& t, const double& dt,
                                           const int id_var, const std::string& name,
                                           const mfem::Vector& u, const double lower_bound,
                                           const double upper_bound) {
  Catch_Time_Section("OperatorBase::ComputeIntegral");

  mfem::ParGridFunction gf(this->fes_[id_var]);
  // Values outside [lower_bound, upper_bound] set to zero, on the host
  mfem::Vector u_cut(u.Size());
  const mfem::real_t* h_u = u.HostRead();
  mfem::real_t* h_cut = u_cut.HostWrite();
  std::transform(h_u, h_u + u.Size(), h_cut, [&](const mfem::real_t value) {
    return (value > upper_bound || value < lower_bound) ? mfem::real_t(0) : value;
  });
  gf.SetFromTrueDofs(u_cut);

  double integral = 0.0;
  double domain_volume = 0.0;
  mfem::Vector vals;
  const mfem::FiniteElement* fe;
  mfem::ElementTransformation* Tr;

  for (int i = 0; i < this->fes_[id_var]->GetNE(); i++) {
    fe = this->fes_[id_var]->GetFE(i);
    const mfem::IntegrationRule* ir;

    int intorder = 2 * fe->GetOrder() + 3;
    ir = &(mfem::IntRules.Get(fe->GetGeomType(), intorder));

    mfem::real_t int_elem = 0.0;
    mfem::real_t domain_volume_elem = 0.0;

    gf.GetValues(i, *ir, vals);
    Tr = this->fes_[id_var]->GetElementTransformation(i);
    for (int j = 0; j < ir->GetNPoints(); j++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(j);
      Tr->SetIntPoint(&ip);

      int_elem += vals(j) * ip.weight * Tr->Weight();
      domain_volume_elem += ip.weight * Tr->Weight();
    }

    integral += int_elem;
    domain_volume += domain_volume_elem;
  }

  mfem::real_t global_integral = 0.0;
  mfem::real_t global_domain_volume = 0.0;
  MPI_Allreduce(&integral, &global_integral, 1, mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&domain_volume, &global_domain_volume, 1, mfem::MPITypeMap<mfem::real_t>::mpi_type,
                MPI_SUM, MPI_COMM_WORLD);

  integral = global_integral;
  domain_volume = global_domain_volume;

  const double average = integral / domain_volume;
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_integral[-]", integral));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_average[-]", average));
}

/**
 * @brief Energies of the system with the current values of the variables u, from the
 * coefficients of the operator (FreeEnergy and, if given, GradientEnergy of block 0):
 *   Density = int (f + g),  Sigma = 2 int g             (FreeEnergy f(u), GradientEnergy g)
 *   Density = int f,        Sigma = 2 int (f - f(u, 0))  (variational FreeEnergy f(u, grad u))
 * Computed element by element by SlothNLFormIntegrator::system_energies (evaluation by batch at
 * the integration points). Sigma is saved if the energies
 * have an interface part (GradientEnergy given or variational FreeEnergy).
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param it Current iteration
 * @param t Current time
 * @param dt Current time step
 * @param u Current values of the variables (true dofs)
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::ComputeEnergies(const int& it, const double& t, const double& dt,
                                           const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("OperatorBase::ComputeEnergies");
  const std::size_t nb = u.size();
  std::vector<mfem::ParGridFunction> vun;
  vun.reserve(nb);
  for (std::size_t b = 0; b < nb; ++b) {
    mfem::ParGridFunction un(this->fes_[b]);
    un.SetFromTrueDofs(u[b]);
    vun.emplace_back(std::move(un));
  }

  SystemEnergyNLFormIntegrator<Variables<T, DIM>> energies(
      this->geometry_, vun, this->params_ - this->default_p_, this->auxvariables_,
      this->coefficients_);
  MFEM_VERIFY(
      energies.prepare_system_energies(),
      "Error: FreeEnergy coefficient is required for computing energy density of the system.");

  mfem::Array<const mfem::FiniteElement*> el(static_cast<int>(nb));
  mfem::Array<const mfem::Vector*> elfun(static_cast<int>(nb));
  std::vector<mfem::Vector> el_dofs(nb);
  mfem::Array<int> vdofs;
  mfem::real_t local[2] = {0.0, 0.0};
  for (int e = 0; e < this->fes_[0]->GetNE(); ++e) {
    for (std::size_t b = 0; b < nb; ++b) {
      const auto ib = static_cast<int>(b);
      el[ib] = this->fes_[b]->GetFE(e);
      this->fes_[b]->GetElementVDofs(e, vdofs);
      vun[b].GetSubVector(vdofs, el_dofs[b]);
      elfun[ib] = &el_dofs[b];
    }
    mfem::ElementTransformation* Tr = this->fes_[0]->GetElementTransformation(e);
    const auto element = energies.system_energies(el, *Tr, elfun);
    local[0] += element.density;
    local[1] += element.interface;
  }

  mfem::real_t global[2] = {0.0, 0.0};
  MPI_Allreduce(local, global, 2, mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM,
                MPI_COMM_WORLD);
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", global[0]));
  if (energies.has_interface_energy()) {
    this->time_specialized_.emplace(IterationKey(it, dt, t),
                                    SpecializedValue("Sigma[J.m-3]", global[1]));
  }
}

/**
 * @brief Compute the position of an isovalue
 * @param it current iteration
 * @param t current time
 * @param dt current timestep
 * @param u unknown vector
 * @param iso_value value of the solution
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::ComputeIsoVal(const int& it, const double& t, const double& dt,
                                         const int var_id, const std::string& var_name,
                                         const mfem::Vector& u, const double& iso_value) {
  Catch_Time_Section("OperatorBase::ComputeIsoVal");
  std::vector<std::string> vstr = {"x", "y", "z"};
  std::multimap<IterationKey, SpecializedValue> variable_time_iso_specialized;

  if (!this->time_iso_specialized_[var_name].empty()) {
    variable_time_iso_specialized.merge(this->time_iso_specialized_[var_name]);
  }

  // Values at the local dofs (element dofs are local indices), read on the host
  mfem::ParGridFunction gf(this->fes_[var_id]);
  gf.SetFromTrueDofs(u);
  const mfem::real_t* h_gf = gf.HostRead();
  std::vector<mfem::Vector> iso_points;

  for (int i = 0; i < this->fes_[var_id]->GetNE(); i++) {
    const mfem::FiniteElement* el = this->fes_[var_id]->GetFE(i);
    mfem::Array<int> dofs;
    this->fes_[var_id]->GetElementDofs(i, dofs);

    mfem::DenseMatrix dof_coords;
    mfem::ElementTransformation* Tr = this->fes_[var_id]->GetElementTransformation(i);
    Tr->Transform(el->GetNodes(), dof_coords);
    std::vector<double> dof_val;
    dof_val.reserve(dofs.Size());
    for (int j = 0; j < dofs.Size(); j++) {
      // Local dofs of the element: values of the grid function (u holds the true dofs, whose
      // numbering differs from the local one in parallel and on nonconforming meshes)
      dof_val.emplace_back(h_gf[dofs[j]] - iso_value);
    }

    // At least one positive and one negative value
    bool has_positive_value =
        std::any_of(dof_val.begin(), dof_val.end(), [](double x) { return x > 0; });
    bool has_negative_value =
        std::any_of(dof_val.begin(), dof_val.end(), [](double x) { return x < 0; });

    if (has_positive_value && has_negative_value) {
      for (int j = 0; j < dofs.Size(); j++) {
        mfem::Vector coord1(DIM);
        double val1 = dof_val[j];
        dof_coords.GetColumn(j, coord1);
        for (int k = j + 1; k < dofs.Size(); k++) {
          double val2 = dof_val[k];

          if (val1 * val2 < 0) {
            const double abs_val1 = std::abs(val1);
            const double abs_val2 = std::abs(val2);
            double t = abs_val1 / (abs_val1 + abs_val2);

            mfem::Vector coord2(DIM), iso_coord(DIM);
            iso_coord = mfem::infinity();
            dof_coords.GetColumn(k, coord2);

            for (int d = 0; d < coord1.Size(); d++) {
              iso_coord[d] = coord1[d] + t * (coord2[d] - coord1[d]);
            }
            iso_points.emplace_back(std::move(iso_coord));
          }
        }
      }
    }
  }

  // Points of all the processes, gathered on the process 0 which writes the file
  std::vector<double> local;
  local.reserve(iso_points.size() * DIM);
  for (const auto& point : iso_points) {
    for (int d = 0; d < DIM; d++) local.push_back(point(d));
  }
  const MPI_Comm comm = this->fes_[var_id]->GetComm();
  int rank = 0, nb_ranks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nb_ranks);
  const int local_size = static_cast<int>(local.size());
  std::vector<int> sizes(rank == 0 ? nb_ranks : 0), offsets(rank == 0 ? nb_ranks : 0);
  MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, comm);
  std::vector<double> all;
  if (rank == 0) {
    int total = 0;
    for (int r = 0; r < nb_ranks; r++) {
      offsets[r] = total;
      total += sizes[r];
    }
    all.resize(static_cast<std::size_t>(total));
  }
  MPI_Gatherv(local.data(), local_size, MPI_DOUBLE, all.data(), sizes.data(), offsets.data(),
              MPI_DOUBLE, 0, comm);
  for (std::size_t i = 0; i + DIM <= all.size(); i += DIM) {
    for (int d = 0; d < DIM; d++) {
      variable_time_iso_specialized.emplace(IterationKey(it, dt, t),
                                            SpecializedValue(vstr[d], all[i + d]));
    }
  }
  this->time_iso_specialized_[var_name] = variable_time_iso_specialized;
}

/**
 * @brief Get a of a multimap of integral values calculated at a given iteration (see
 * computeEnergies and computeError)
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM>
const std::map<std::string, std::multimap<IterationKey, SpecializedValue>>
OperatorBase<T, DIM>::get_time_iso_specialized() const {
  return this->time_iso_specialized_;
}

/**
 * @brief Clear time_iso_specialized_ container
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::clear_iso_time_specialized() {
  this->time_iso_specialized_.clear();
}

/**
 * @brief Get a of a multimap of integral values calculated at a given iteration (see
 * computeEnergies and computeError)
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM>
const std::multimap<IterationKey, SpecializedValue>& OperatorBase<T, DIM>::get_time_specialized()
    const {
  return this->time_specialized_;
}

/**
 * @brief Clear time_specialized_ container
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::clear_time_specialized() {
  this->time_specialized_.clear();
}

/**
 * @brief Set the default options for the nonlinear algorithm and associated solvers
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::set_default_solver() {
  auto nl_params =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", false));
  auto s_params = Parameters(Parameter("description", "Default Solver for Newton Algorithm"));
  auto p_params =
      Parameters(Parameter("description", "Default Preconditioner for Newton Algorithm"));

  this->nl_solver_ = NLSolverType::NEWTON;
  this->nl_solver_params_ = nl_params;
  this->solver_ = HypreSolverType::HYPRE_GMRES;
  this->solver_params_ = s_params;
  this->precond_ = HyprePreconditionerType::HYPRE_ILU;
  this->precond_params_ = p_params;
}

/**
 * @brief Return a RHS NLFormIntegrator
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param integrator
 * @param vun
 * @param all_params
 * @return SlothNLFormIntegrator<Variables<T, DIM>>*
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* OperatorBase<T, DIM>::get_rhs_integrator(
    const double dt, const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
    const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params) {
  switch (Integrators::from(integrator)) {
    case Integrators::MassFlux: {
      return new MassDiffusionFluxNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::Fick: {
      return new FickNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::Fourier: {
      return new FourierNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::CahnHilliard: {
      return new CahnHilliardNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::AllenCahn: {
      return new AllenCahnNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::SplitAllenCahn: {
      return new BlockAllenCahnNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::MeltingTemperature: {
      return new MeltingTemperatureNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::MeltingCalphad: {
      return new MeltingCalphadNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::MeltingConstant: {
      return new MeltingConstantNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::LatentHeat: {
      return new LatentHeatNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::SourceTerm: {
      return new SourceTermNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::Advection: {
      return new AdvectionNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::SplitAdvection: {
      return new SplitAdvectionNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    default:
      mfem::mfem_error("RHS Integrators not found. Please check your data.");
  }
}

/**
 * @brief Return a boundary NLFormIntegrator
 *
 * @tparam T
 * @tparam DIM
 * @param integrator
 * @param vun
 * @param vauxn
 * @param all_params
 * @param block
 * @param bdr_id
 * @return SlothNLFormIntegrator<Variables<T, DIM>>*
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* OperatorBase<T, DIM>::get_bdr_integrator(
    const double dt, const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
    const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params,
    const unsigned int block, const unsigned int bdr_id) {
  switch (Integrators::from(integrator)) {
    case Integrators::Neumann: {
      return new NeumannNLFormIntegrator<Variables<T, DIM>>(this->geometry_, dt, vun, vauxn,
                                                            all_params, this->auxvariables_,
                                                            this->coefficients_, block, bdr_id);
      break;
    }
    case Integrators::Robin: {
      return new RobinNLFormIntegrator<Variables<T, DIM>>(this->geometry_, dt, vun, vauxn,
                                                          all_params, this->auxvariables_,
                                                          this->coefficients_, block, bdr_id);
      break;
    }
    case Integrators::AdvectionBoundary: {
      // AdvectionBoundary added on boundaries holding an Inflow coefficient.
      return new AdvectionBoundaryNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_,
          block, bdr_id);
    }
    case Integrators::AdvectionBoundaryCH: {
      // AdvectionBoundaryCH added on boundaries holding an Inflow coefficient.
      return new AdvectionBoundaryCHNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_,
          block, bdr_id);
    }
    default:
      mfem::mfem_error("Boundary Integrators not found. Please check your data.");
  }
}

/**
 * @brief Overload the nonlinear algorithm
 * @remark NewtonSolver is the only one implemented
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param NLSOLVER
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::overload_nl_solver(NLSolverType NLSOLVER) {
  this->nl_solver_ = NLSOLVER;
}

/**
 * @brief Overload Parameters associated with nonlinear algorithm
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param NLSOLVER
 * @param nl_params
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::overload_nl_solver(NLSolverType NLSOLVER, const Parameters& nl_params) {
  this->nl_solver_ = NLSOLVER;
  this->nl_solver_params_ = nl_params;
}

/**
 * @brief  Overload the default linear solver used for the LHS
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::overload_solver(VSolverType SOLVER) {
  this->solver_ = SOLVER;
}

/**
 * @brief  Overload the parameter for the linear solver used for the LHS
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::overload_solver(VSolverType SOLVER, const Parameters& s_params) {
  this->solver_ = SOLVER;
  this->solver_params_ = s_params;
}

/**
 * @brief Overload the preconditioner used by the solver in the NL algorithm.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param PRECOND
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::overload_preconditioner(VSolverType PRECOND) {
  this->precond_ = PRECOND;
}

/**
 * @brief  Overload the parameters of the preconditioner used by the solver in the NL algorithm.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::overload_preconditioner(VSolverType PRECOND,
                                                   const Parameters& p_params) {
  this->precond_ = PRECOND;
  this->precond_params_ = p_params;
}

/**
 * @brief Set the time for all coefficients
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param time current time
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::set_time_coefficients(double time) {
  for (auto& coefs : this->coefficients_) {
    coefs.set_time(time);
  }
}

/**
 * @brief Define the geometry of the problem (Cartesian, Axisymmetric)
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param geometry
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::setGeometry(Geometry geometry) {
  this->geometry_ = geometry;
}

/**
 * @brief Resynchronize the operator's structural state (height, width,
 *        internal work vector, and block offsets) after the underlying
 *        mesh/finite element spaces have changed size (e.g. following AMR).
 *
 * @details `this->fes_` holds raw `mfem::ParFiniteElementSpace*` pointers
 *          that are updated in place elsewhere (via `fespace->Update()`),
 *          so reading their current `GetTrueVSize()` here always reflects
 *          the latest mesh state. However, `mfem::Operator::height`/`width`,
 *          the internal work vector `z`, and `block_trueOffsets_` are only
 *          computed once in the constructor and are NOT automatically kept
 *          in sync when the mesh changes size later. This method must be
 *          called after every mesh refinement/derefinement cycle, before
 *          the operator is used again (e.g. before the next Newton solve),
 *          otherwise stale sizes silently propagate into downstream solves
 *          (mismatched block/vector sizes, out-of-bounds reads).
 *
 *          The lists of essential true dofs (Dirichlet boundaries) are
 *          recomputed as well: the true dof numbering changes with the mesh,
 *          and these lists are used to eliminate the constrained dofs in the
 *          reduced operators (held by reference) and in the mass matrices.
 *
 * @note This only updates the `height`/`width` copy inherited through
 *       `OperatorBase`. Classes that also inherit `mfem::TimeDependentOperator`
 *       (a second, independent `mfem::Operator` base due to non-virtual
 *       multiple inheritance) must separately resynchronize that copy too.
 *
 * @tparam T Finite element collection type (mfem object).
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::UpdateAfterMeshChange() {
  int total_size = 0;
  for (int i = 0; i < this->fes_.Size(); i++) {
    total_size += this->fes_[i]->GetTrueVSize();
  }

  // MFEM (mfem::TimeDependentOperator)
  this->height = total_size;
  this->width = total_size;
  // SLOTH (OperatorBase)
  this->height_ = total_size;

  this->z.SetSize(total_size);

  std::ostringstream msg;
  msg << "OperatorBase::UpdateAfterMeshChange: block_trueOffsets_ size ("
      << this->block_trueOffsets_.Size() << ") does not match fes_.Size()+1 ("
      << this->fes_.Size() + 1 << ").";
  MFEM_VERIFY(this->block_trueOffsets_.Size() == this->fes_.Size() + 1, msg.str());

  this->block_trueOffsets_[0] = 0;
  for (int i = 0; i < this->fes_.Size(); i++) {
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
  }
  this->block_trueOffsets_.PartialSum();

  // Essential true dofs on the new numbering (updated in place: the reduced operators
  // keep a reference to ess_tdof_list_)
  MFEM_VERIFY(this->bcs_.size() == this->ess_tdof_list_.size(),
              "OperatorBase::UpdateAfterMeshChange: bcs_ and ess_tdof_list_ sizes differ.");
  for (std::size_t i = 0; i < this->bcs_.size(); i++) {
    this->bcs_[i]->UpdateEssentialDofs();
    this->ess_tdof_list_[i] = this->bcs_[i]->GetEssentialDofs();
  }
}
