/**
 * @file DiffusionOperator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Diffusion operator (Base, Steady and TimeDependent)
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

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Operators/SteadyOperatorBase.hpp"
#include "Operators/TransientOperatorBase.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

/**
 * @brief Base class for Mass Diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, template <class, int> class OPEBASE>
class DiffusionOperatorBase : public OPEBASE<T, DIM> {
 protected:
  Parameters diffusion_params_;

 public:
  template <typename... Args>
  DiffusionOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                        const std::vector<std::string>& integrators, Args&&... args)
      : OPEBASE<T, DIM>(spatials, integrators, std::forward<Args>(args)...) {}

  template <typename... Args>
  DiffusionOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                        const std::vector<std::string>& integrators, const Parameters& params,
                        Args&&... args)
      : OPEBASE<T, DIM>(spatials, integrators, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;

  SlothNLFormIntegrator<Variables<T, DIM>>* set_nlfi_ptr(
      const std::string nlfi, const double dt, const std::vector<mfem::Vector>& u) override;
  void get_parameters() override;

  ~DiffusionOperatorBase();
};

/**
 * @brief Destroy the DiffusionOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, template <class, int> class OPEBASE>
DiffusionOperatorBase<T, DIM, OPEBASE>::~DiffusionOperatorBase() {}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief  Set the NonLinearFormIntegrator dedicated to diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM, template <class, int> class OPEBASE>
SlothNLFormIntegrator<Variables<T, DIM>>* DiffusionOperatorBase<T, DIM, OPEBASE>::set_nlfi_ptr(
    const std::string nlfi, const double dt, const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("DiffusionOperatorBase::set_nlfi_ptr");

  std::vector<mfem::ParGridFunction> vun;
  for (unsigned int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }

  const Parameters& all_params = this->diffusion_params_ + this->params_ - this->default_p_;

  auto rhs_nlfi = this->get_rhs_integrator(nlfi, vun, all_params);
  rhs_nlfi->init();
  return rhs_nlfi;
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
template <class T, int DIM, template <class, int> class OPEBASE>
void DiffusionOperatorBase<T, DIM, OPEBASE>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description", "Diffusion operator");
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class SteadyDiffusionOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM>
class SteadyDiffusionOperator final : public DiffusionOperatorBase<T, DIM, SteadyOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyDiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                          const std::vector<std::string>& integrators, const Parameters& params,
                          Args&&... args)
      : DiffusionOperatorBase<T, DIM, SteadyOperatorBase>(spatials, integrators, params,
                                                          std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyDiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                          const std::vector<std::string>& integrators, Args&&... args)
      : DiffusionOperatorBase<T, DIM, SteadyOperatorBase>(spatials, integrators,
                                                          std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  void overload_diffusion(const Parameters& p_params);

  ~SteadyDiffusionOperator() {}
};
/**
 * @brief Set the default options for properties
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM>
void SteadyDiffusionOperator<T, DIM>::set_default_properties() {
  this->diffusion_params_ = Parameters(Parameter("D", 1.0), Parameter("D_stab", 0.0));
}

/**
 * @brief Overload the default options for diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param DIFFU
 * @param p_params
 */
template <class T, int DIM>
void SteadyDiffusionOperator<T, DIM>::overload_diffusion(const Parameters& p_params) {
  this->diffusion_params_ = p_params;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class DiffusionOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM>
class DiffusionOperator final : public DiffusionOperatorBase<T, DIM, TransientOperatorBase> {
 protected:
  Parameters density_params_;
  void set_default_properties() override;
  void get_mass_coefficient(const mfem::Vector& u) override;

 public:
  template <typename... Args>
  DiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& integrators, const Parameters& params,
                    Args&&... args)
      : DiffusionOperatorBase<T, DIM, TransientOperatorBase>(spatials, integrators, params,
                                                             std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  DiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& integrators, Args&&... args)
      : DiffusionOperatorBase<T, DIM, TransientOperatorBase>(spatials, integrators,
                                                             std::forward<Args>(args)...) {
    this->set_default_properties();
  }

  void overload_density(const Parameters& p_params);
  void overload_diffusion(const Parameters& p_params);

  ~DiffusionOperator() {}
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
template <class T, int DIM>
void DiffusionOperator<T, DIM>::get_mass_coefficient(const mfem::Vector& u) {
  if (this->MassCoeff_ != nullptr) {
    delete this->MassCoeff_;
  }
  if (this->mass_gf_ != nullptr) {
    delete this->mass_gf_;
  }
  this->mass_gf_ = new mfem::ParGridFunction(this->fes_[0]);
  this->mass_gf_->SetFromTrueDofs(u);
  // auto density_coefficient = new DensityCoefficient<0, DENS>(this->mass_gf_,
  // this->density_params_);
  auto one = new mfem::ConstantCoefficient(1.0);
  // this->MassCoeff_ = new mfem::ProductCoefficient(*density_coefficient, *one);
  this->MassCoeff_ = new mfem::ProductCoefficient(*one, *one);
}

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM>
void DiffusionOperator<T, DIM>::set_default_properties() {
  this->density_params_ = Parameters(Parameter("rho", 1.0));

  this->diffusion_params_ = Parameters(Parameter("D", 1.0));
}
/**
 * @brief Overload the default options for density
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM>
void DiffusionOperator<T, DIM>::overload_density(const Parameters& p_params) {
  this->density_params_ = p_params;
}

/**
 * @brief Overload the default options for diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM>
void DiffusionOperator<T, DIM>::overload_diffusion(const Parameters& p_params) {
  this->diffusion_params_ = p_params;
}
