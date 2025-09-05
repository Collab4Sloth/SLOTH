/**
 * @file PhaseFieldOperator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief AllenCahn operator (Base, Steady and TimeDependent)
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
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
#include <memory>
#include <utility>
#include <vector>

#include "Operators/SteadyOperatorBase.hpp"
#include "Operators/TransientOperatorBase.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Base class for AllenCahn
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
class PhaseFieldOperatorBase : public OPEBASE<T, DIM, NLFI, LHS_NLFI> {
 protected:
  Parameters mobility_params_;
  double omega_, lambda_;

 public:
  template <typename... Args>
  explicit PhaseFieldOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                  Args&&... args)
      : OPEBASE<T, DIM, NLFI, LHS_NLFI>(spatials, std::forward<Args>(args)...) {
    this->get_parameters();
  }
  template <typename... Args>
  PhaseFieldOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                         const Parameters& params, Args&&... args)
      : OPEBASE<T, DIM, NLFI, LHS_NLFI>(spatials, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;
  NLFI* set_nlfi_ptr(const double dt, const std::vector<mfem::Vector>& u) override;
  void get_parameters() override;
  void ComputeEnergies(const int& it, const double& t, const double& dt,
                       const std::vector<mfem::Vector>& u) override;

  ~PhaseFieldOperatorBase();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Destroy the PhaseFieldOperatorBase< T, DIM, NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::~PhaseFieldOperatorBase() {}

/**
 * @brief Set the NonLinearFormIntegrator dedicated to AllenCahn
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
NLFI* PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::set_nlfi_ptr(
    const double dt, const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("PhaseFieldOperatorBase::set_nlfi_ptr");
  std::vector<mfem::ParGridFunction> vun;
  for (int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }

  const Parameters& all_params = this->mobility_params_ + this->params_ - this->default_p_;

  NLFI* nlfi_ptr = new NLFI(vun, all_params, this->auxvariables_);

  return nlfi_ptr;
}

/**
 * @brief Get parameters for use in the current operator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
void PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::get_parameters() {
  this->omega_ = this->params_.template get_param_value<double>("omega");
  this->lambda_ = this->params_.template get_param_value<double>("lambda");
}

/**
 * @brief  Compute Energy contribution at the end of time step
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param it
 * @param dt
 * @param t
 * @param u
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
void PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::ComputeEnergies(
    const int& it, const double& t, const double& dt, const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("PhaseFieldOperatorBase::ComputeEnergies");

  std::vector<mfem::ParGridFunction> vun;
  vun.reserve(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(std::move(un));
  }

  std::vector<mfem::ParGridFunction*> vun_ptr;
  for (auto& un : vun) {
    vun_ptr.push_back(&un);
  }

  auto gphi = this->nlfi_ptr_->get_energy(vun_ptr, this->omega_);

  std::unique_ptr<GradientEnergyCoefficient> gradGphi =
      this->nlfi_ptr_->get_grad_energy(vun_ptr, 0.5 * this->lambda_);
  auto energy_coeff = gphi.get();

  mfem::ParGridFunction gf(this->fes_[0]);
  gf.ProjectCoefficient(*energy_coeff);

  mfem::ParGridFunction sigf(this->fes_[0]);
  auto grad_energy_coeff = gradGphi.get();
  sigf.ProjectCoefficient(*grad_energy_coeff);

  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeLpError(1, zero);
  const auto interfacial_energy = sigf.ComputeLpError(1, zero);

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", energy + interfacial_energy));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Sigma[J.m-3]", 2. * interfacial_energy));
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class SteadyPhaseFieldOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI = mfem::BlockNonlinearFormIntegrator>
class SteadyPhaseFieldOperator final
    : public PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyPhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials, Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyOperatorBase>(
            spatials, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyPhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                           const Parameters& params, Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyOperatorBase>(
            spatials, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  void overload_mobility(const Parameters& p_params);

  ~SteadyPhaseFieldOperator() {}
};

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyPhaseFieldOperator<T, DIM, NLFI, LHS_NLFI>::set_default_properties() {
  this->mobility_params_ = Parameters(Parameter("mob", 1.0));
}
/**
 * @brief Overload the default options for mobility
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyPhaseFieldOperator<T, DIM, NLFI, LHS_NLFI>::overload_mobility(
    const Parameters& p_params) {
  this->mobility_params_ = p_params;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class PhaseFieldOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
class PhaseFieldOperator final
    : public PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, TransientOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  PhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials, Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, TransientOperatorBase>(
            spatials, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  PhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params,
                     Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI, TransientOperatorBase>(
            spatials, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }

  void overload_mobility(const Parameters& p_params);
  void get_mass_coefficient(const mfem::Vector& u) override;

  ~PhaseFieldOperator() {}
};
/**
 * @brief Overload the MassMatrix coefficient definition
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param u
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void PhaseFieldOperator<T, DIM, NLFI, LHS_NLFI>::get_mass_coefficient(const mfem::Vector& u) {
  TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::get_mass_coefficient(u);
}

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void PhaseFieldOperator<T, DIM, NLFI, LHS_NLFI>::set_default_properties() {
  this->mobility_params_ = Parameters(Parameter("mob", 1.0));
}
/**
 * @brief Overload the default options for mobility
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void PhaseFieldOperator<T, DIM, NLFI, LHS_NLFI>::overload_mobility(const Parameters& p_params) {
  this->mobility_params_ = p_params;
}
