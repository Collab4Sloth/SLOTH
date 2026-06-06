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
 * @brief Set the NonLinearFormIntegrator dedicated to AllenCahn
 *
 * @tparam T
 * @tparam DIM
 * @tparam OPEBASE
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* OperatorBase<T, DIM>::set_nlfi_ptr(
    const std::string nlfi, const std::vector<mfem::Vector>& u) {
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
  auto rhs_nlfi = this->get_rhs_integrator(nlfi, vun, vauxn, all_params);
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
    const std::string nlfi, const std::vector<mfem::Vector>& u, const unsigned int block,
    const unsigned int bdr_id) {
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
  auto bdr_nlfi = this->get_bdr_integrator(nlfi, vun, vauxn, all_params, block, bdr_id);
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
void OperatorBase<T, DIM>::set_coefficients(const std::vector<Coefficients>& coefficients,
                                            bool enable_compute_energies) {
  this->coefficients_ = coefficients;

  if (enable_compute_energies) {
    Coefficients coefs = this->coefficients_[0];
    for (unsigned int i = 0; i < coefs.size(); i++) {
      if (coefs[i].get_type() == GlossaryType::FreeEnergy) {
        this->energy_coefficient_ = coefs[i];
      }
      if (coefs[i].get_type() == GlossaryType::GradientEnergy) {
        this->grad_energy_coefficient_ = coefs[i];
      }
    }

    MFEM_VERIFY(
        this->energy_coefficient_.has_value(),
        "Error: FreeEnergy coefficient is required for computing energy density of the system.");
  }
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatial
 */
template <class T, int DIM>
OperatorBase<T, DIM>::OperatorBase(const std::vector<std::string>& integrators,
                                   std::vector<SpatialDiscretization<T, DIM>*> spatials)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(default_params_),
      RHS(NULL),
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
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatial
 * @param source_term_name
 */
template <class T, int DIM>
OperatorBase<T, DIM>::OperatorBase(const std::vector<std::string>& integrators,
                                   std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                   const std::vector<AnalyticalFunctions<DIM>>& source_term_name)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(default_params_),
      RHS(NULL),
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
  for (const auto& src : source_term_name) {
    this->src_func_.emplace_back(std::move(src.getFunction()));
  }
  this->set_default_solver();
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
OperatorBase<T, DIM>::OperatorBase(const std::vector<std::string>& integrators,
                                   std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                   const Parameters& params)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(params),
      RHS(NULL),
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
}

/**
 * @brief  Construct a new Phase Field Operator Base< T,  DIM, NLFI>:: Phase Field Operator
 * Base object
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param spatial
 * @param params
 * @param vars
 * @param auxvars
 * @param source_term_name
 */
template <class T, int DIM>
OperatorBase<T, DIM>::OperatorBase(const std::vector<std::string>& integrators,
                                   std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                   const Parameters& params,
                                   const std::vector<AnalyticalFunctions<DIM>>& source_term_name)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(params),
      RHS(NULL),
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
  for (const auto& src : source_term_name) {
    auto s = src.getFunction();
    this->src_func_.emplace_back(s);
  }
  this->set_default_solver();
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
void OperatorBase<T, DIM>::initialize([[maybe_unused]] const double& initial_time,
                                      Variables<T, DIM>& vars,
                                      std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("OperatorBase::initialize");

  this->auxvariables_ = auxvars;
  const auto nvars = vars.get_variables_number();
  std::vector<mfem::Vector> u_vect;
  u_vect.reserve(nvars);
  for (unsigned int iv = 0; iv < nvars; iv++) {
    auto& vv = vars[iv];
    auto u = vv.get_unknown();

    this->bcs_.emplace_back(vv.get_boundary_conditions());
    this->ess_tdof_list_.emplace_back(this->bcs_[iv]->GetEssentialDofs());
    this->bcs_[iv]->SetBoundaryConditions(u);
    vv.update(u);
    u_vect.emplace_back(u);
  }
  this->SetTransientParameters(u_vect);
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
void OperatorBase<T, DIM>::build_rhs_nonlinear_form(const std::vector<mfem::Vector>& u_vect) {
  if (this->RHS != nullptr) {
    delete this->RHS;
  }
  this->RHS = new mfem::ParBlockNonlinearForm(this->fes_);
  for (const std::string& s_integrator : this->rhs_integrators_) {
    auto integrator_ptr = this->set_nlfi_ptr(s_integrator, u_vect);
    this->RHS->AddDomainIntegrator(integrator_ptr);
  }

  // Add boundary integrators
  const int fes_size = this->block_trueOffsets_.Size() - 1;
  mfem::Array<int> Robin_bdr;
  mfem::Array<int> Neumann_bdr;
  // Loop over variables
  for (int i = 0; i < fes_size; ++i) {
    Robin_bdr = this->bcs_[i]->get_marker_array("Robin");
    Neumann_bdr = this->bcs_[i]->get_marker_array("Neumann");

    const int nb_bdr = Robin_bdr.Size();
    this->array_bdr_.SetSize(nb_bdr);
    Coefficients coefficients = this->coefficients_[i];
    const int coef_size = coefficients.size();

    // Loop over boundaries
    for (auto j = 0; j < nb_bdr; j++) {
      // Add Neumann integrator
      if (Neumann_bdr[j] > 0) {
        // Check if a coefficient is given for this bc
        // (else Homogeneous Neumann)
        bool has_neumann_coeff = false;
        for (int l = 0; l < coef_size; l++) {
          auto coef = coefficients[l];
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
          this->array_bdr_ = 0;
          this->array_bdr_[j] = 1;
          auto integrator_ptr = this->set_bdr_nlfi_ptr("Neumann", u_vect, i, j);
          this->RHS->AddBoundaryIntegrator(integrator_ptr, this->array_bdr_);
        }
      }
      // Add Robin integrator
      if (Robin_bdr[j] > 0) {
        this->array_bdr_ = 0;
        this->array_bdr_[j] = 1;
        auto integrator_ptr = this->set_bdr_nlfi_ptr("Robin", u_vect, i, j);
        this->RHS->AddBoundaryIntegrator(integrator_ptr, this->array_bdr_);
      }
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
  mfem::Vector u_cut(u.Size());
  std::transform(u.begin(), u.end(), u_cut.begin(), [&](auto value) {
    if (value > upper_bound || value < lower_bound) {
      return 0.0;
    }
    return value;
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

    int intorder = 2 * fe->GetOrder() + 3;  // <----------
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
 * @brief Compute energy density
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param it
 * @param t
 * @param dt
 * @param u
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::ComputeEnergies(const int& it, const double& t, const double& dt,
                                           const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("OperatorBase::ComputeEnergies");
  std::vector<mfem::ParGridFunction> vun;
  vun.reserve(u.size());
  for (size_t j = 0; j < u.size(); ++j) {
    mfem::ParGridFunction un(this->fes_[j]);
    un.SetFromTrueDofs(u[j]);
    vun.emplace_back(std::move(un));
  }

  std::vector<mfem::ParGridFunction> vaux;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      vaux.emplace_back(auxvar.get_gf());
    }
  }
  const int m = vun.size();
  const int l = vaux.size();
  std::vector<mfem::real_t> v1(m);
  std::vector<mfem::real_t> v2(l);

  mfem::ParGridFunction gf(this->fes_[0]);
  const bool semi_implicit_coef = (*this->energy_coefficient_).is_semi_implicit();

  for (int i = 0; i < gf.Size(); ++i) {
    for (int j = 0; j < m; ++j) v1[j] = vun[j][i];
    for (int j = 0; j < l; ++j) v2[j] = vaux[j][i];
    gf[i] = semi_implicit_coef ? (*this->energy_coefficient_).compute(v1, v1, v2, -1)
                               : (*this->energy_coefficient_).compute(v1, v2, -1);
  }

  std::vector<mfem::ParGridFunction> vgf;
  vgf.reserve(m);
  std::vector<mfem::ParGridFunction> vgf_aux;
  vgf_aux.reserve(l);
  for (int k = 0; k < m; ++k) {
    mfem::ParGridFunction gf(this->fes_[0]);
    for (int i = 0; i < gf.Size(); ++i) {
      gf[i] = vun[k][i];
    }
    vgf.emplace_back(std::move(gf));
  }
  for (int k = 0; k < l; ++k) {
    mfem::ParGridFunction gf(this->fes_[0]);
    for (int i = 0; i < gf.Size(); ++i) {
      gf[i] = vaux[k][i];
    }
    vgf_aux.emplace_back(std::move(gf));
  }

  mfem::real_t integral = 0.0;
  mfem::real_t sigma = 0.0;
  const mfem::FiniteElement* fe;
  mfem::ElementTransformation* Tr;
  mfem::Vector vals;

  for (int i = 0; i < this->fes_[0]->GetNE(); i++) {
    fe = this->fes_[0]->GetFE(i);

    unsigned int dim = fe->GetDim();
    std::vector<mfem::real_t> grad_vun_ip(m * dim);
    std::vector<mfem::real_t> grad_vaux_ip(l * dim);
    mfem::Vector grad_tmp;
    grad_tmp.SetSize(dim);
    const mfem::IntegrationRule* ir;

    int intorder = 2 * fe->GetOrder() + 3;  // <----------
    ir = &(mfem::IntRules.Get(fe->GetGeomType(), intorder));

    mfem::real_t int_fo = 0.0;
    mfem::real_t int_grad = 0.0;

    gf.GetValues(i, *ir, vals);

    Tr = this->fes_[0]->GetElementTransformation(i);

    for (int j = 0; j < ir->GetNPoints(); j++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(j);
      Tr->SetIntPoint(&ip);

      // CCI Grad contrib
      if (this->grad_energy_coefficient_.has_value()) {
        for (unsigned int k = 0; k < vgf.size(); ++k) {
          vgf[k].GetGradient(*Tr, grad_tmp);
          for (unsigned int ii = 0; ii < dim; ii++) grad_vun_ip[k * dim + ii] = grad_tmp[ii];
        }
        for (unsigned int k = 0; k < vgf_aux.size(); ++k) {
          vgf_aux[k].GetGradient(*Tr, grad_tmp);
          for (unsigned int ii = 0; ii < dim; ii++) grad_vaux_ip[k * dim + ii] = grad_tmp[ii];
        }
        int_grad += (*this->grad_energy_coefficient_).compute(grad_vun_ip, grad_vaux_ip, dim) *
                    ip.weight * Tr->Weight();
      }
      int_fo += vals(j) * ip.weight * Tr->Weight();
    }
    integral += int_fo + int_grad;
    sigma += 2.0 * int_grad;
  }

  mfem::real_t global_integral = 0.0;

  MPI_Allreduce(&integral, &global_integral, 1, mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM,
                MPI_COMM_WORLD);

  integral = global_integral;

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", integral));
  if (this->grad_energy_coefficient_.has_value()) {
    mfem::real_t global_sigma = 0.0;
    MPI_Allreduce(&sigma, &global_sigma, 1, mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM,
                  MPI_COMM_WORLD);

    sigma = global_sigma;
    this->time_specialized_.emplace(IterationKey(it, dt, t),
                                    SpecializedValue("Sigma[J.m-3]", sigma));
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

  mfem::ParGridFunction gf(this->fes_[var_id]);
  gf.SetFromTrueDofs(u);
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
      dof_val.emplace_back(u(dofs[j]) - iso_value);
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
      for (size_t i = 0; i < iso_points.size(); i++) {
        for (size_t d = 0; d < DIM; d++) {
          variable_time_iso_specialized.emplace(IterationKey(it, dt, t),
                                                SpecializedValue(vstr[d], iso_points[i](d)));
        }
      }
    }
  }
  this->time_iso_specialized_[var_name] = variable_time_iso_specialized;
}

/**
 * @brief Get the source term by equation of the PDEs
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void OperatorBase<T, DIM>::get_source_term(
    const int id_block, const std::function<double(const mfem::Vector&, double)>& src_func,
    mfem::Vector& source_term, mfem::ParLinearForm* RHSS) const {
  mfem::FunctionCoefficient src(src_func);
  src.SetTime(this->current_time_ + this->current_dt_);

  RHSS->AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
  RHSS->Assemble();

  source_term.SetSize(this->fes_[id_block]->GetTrueVSize());
  RHSS->ParallelAssemble(source_term);

  source_term.SetSubVector(this->ess_tdof_list_[id_block], 0.);
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
    const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
    const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params) {
  switch (Integrators::from(integrator)) {
    case Integrators::MassFlux: {
      return new MassDiffusionFluxNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::Fick: {
      return new FickNLFormIntegrator<Variables<T, DIM>>(this->geometry_, vun, vauxn, all_params,
                                                         this->auxvariables_, this->coefficients_);
    }
    case Integrators::Fourier: {
      return new FourierNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::CahnHilliard: {
      return new CahnHilliardNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::AllenCahn: {
      return new AllenCahnNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::SplitAllenCahn: {
      return new BlockAllenCahnNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::MeltingTemperature: {
      return new MeltingTemperatureNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::MeltingCalphad: {
      return new MeltingCalphadNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::MeltingConstant: {
      return new MeltingConstantNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
    }
    case Integrators::LatentHeat: {
      return new LatentHeatNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
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
    const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
    const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params,
    const unsigned int block, const unsigned int bdr_id) {
  switch (Integrators::from(integrator)) {
    case Integrators::Neumann: {
      return new NeumannNLFormIntegrator<Variables<T, DIM>>(this->geometry_, vun, vauxn, all_params,
                                                            this->auxvariables_,
                                                            this->coefficients_, block, bdr_id);
      break;
    }
    case Integrators::Robin: {
      return new RobinNLFormIntegrator<Variables<T, DIM>>(this->geometry_, vun, vauxn, all_params,
                                                          this->auxvariables_, this->coefficients_,
                                                          block, bdr_id);
      break;
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
 * @remark only used for explicit time-scheme
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
 * @remark only used for explicit time-scheme
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
  if (this->energy_coefficient_.has_value()) {
    (*this->energy_coefficient_).set_time(time);
  }
  if (this->grad_energy_coefficient_.has_value()) {
    (*this->grad_energy_coefficient_).set_time(time);
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
