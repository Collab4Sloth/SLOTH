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

#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Operators/SteadyOperatorBase.hpp"
#include "Operators/TransientOperatorBase.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
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
template <class T, int DIM, template <class, int> class OPEBASE>
class PhaseFieldOperatorBase : public OPEBASE<T, DIM> {
 protected:
  Parameters mobility_params_;
  double omega_, lambda_;

 public:
  template <typename... Args>
  explicit PhaseFieldOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                  const std::vector<std::string>& integrators, Args&&... args)
      : OPEBASE<T, DIM>(spatials, integrators, std::forward<Args>(args)...) {
    this->get_parameters();
  }
  template <typename... Args>
  PhaseFieldOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                         const std::vector<std::string>& integrators, const Parameters& params,
                         Args&&... args)
      : OPEBASE<T, DIM>(spatials, integrators, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;
  SlothNLFormIntegrator<Variables<T, DIM>>* set_nlfi_ptr(
      const std::string nlfi, const double dt, const std::vector<mfem::Vector>& u) override;
  void get_parameters() override;

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
 * @tparam OPEBASE
 */
template <class T, int DIM, template <class, int> class OPEBASE>
PhaseFieldOperatorBase<T, DIM, OPEBASE>::~PhaseFieldOperatorBase() {}

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
template <class T, int DIM, template <class, int> class OPEBASE>
SlothNLFormIntegrator<Variables<T, DIM>>* PhaseFieldOperatorBase<T, DIM, OPEBASE>::set_nlfi_ptr(
    const std::string nlfi, const double dt, const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("PhaseFieldOperatorBase::set_nlfi_ptr");
  std::vector<mfem::ParGridFunction> vun;
  for (unsigned int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }

  const Parameters& all_params = this->mobility_params_ + this->params_ - this->default_p_;
  auto rhs_nlfi = this->get_rhs_integrator(nlfi, vun, all_params);
  rhs_nlfi->init();
  return rhs_nlfi;
}

/**
 * @brief Get parameters for use in the current operator
 *
 * @tparam T
 * @tparam DIM
 * @tparam OPEBASE
 * @param params
 */
template <class T, int DIM, template <class, int> class OPEBASE>
void PhaseFieldOperatorBase<T, DIM, OPEBASE>::get_parameters() {
  // this->omega_ = this->params_.template get_param_value<double>("omega");
  // this->lambda_ = this->params_.template get_param_value<double>("lambda");
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class SteadyPhaseFieldOperator
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
class SteadyPhaseFieldOperator final : public PhaseFieldOperatorBase<T, DIM, SteadyOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyPhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                           const std::vector<std::string>& integrators, Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, SteadyOperatorBase>(spatials, integrators,
                                                           std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyPhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                           const std::vector<std::string>& integrators, const Parameters& params,
                           Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, SteadyOperatorBase>(spatials, integrators, params,
                                                           std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  void overload_mobility(const Parameters& p_params);

  ~SteadyPhaseFieldOperator() {}
};

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void SteadyPhaseFieldOperator<T, DIM>::set_default_properties() {
  this->mobility_params_ = Parameters(Parameter("mob", 1.0));
}
/**
 * @brief Overload the default options for mobility
 *
 * @tparam T
 * @tparam DIM
 * @param p_params
 */
template <class T, int DIM>
void SteadyPhaseFieldOperator<T, DIM>::overload_mobility(const Parameters& p_params) {
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
 */
template <class T, int DIM>
class PhaseFieldOperator final : public PhaseFieldOperatorBase<T, DIM, TransientOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  PhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const std::vector<std::string>& integrators, Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, TransientOperatorBase>(spatials, integrators,
                                                              std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  PhaseFieldOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const std::vector<std::string>& integrators, const Parameters& params,
                     Args&&... args)
      : PhaseFieldOperatorBase<T, DIM, TransientOperatorBase>(spatials, integrators, params,
                                                              std::forward<Args>(args)...) {
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
 * @tparam OPEBASE
 * @param u
 */
template <class T, int DIM>
void PhaseFieldOperator<T, DIM>::get_mass_coefficient(const mfem::Vector& u) {
  TransientOperatorBase<T, DIM>::get_mass_coefficient(u);
}

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void PhaseFieldOperator<T, DIM>::set_default_properties() {
  this->mobility_params_ = Parameters(Parameter("mob", 1.0));
}
/**
 * @brief Overload the default options for mobility
 *
 * @tparam T
 * @tparam DIM
 * @param p_params
 */
template <class T, int DIM>
void PhaseFieldOperator<T, DIM>::overload_mobility(const Parameters& p_params) {
  this->mobility_params_ = p_params;
}
