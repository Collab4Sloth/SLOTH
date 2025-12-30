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

  SlothNLFormIntegrator<Variables<T, DIM>>* set_nlfi_ptr(
      const std::string nlfi, const std::vector<mfem::Vector>& u) override;
  void get_parameters() override;
};

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
    const std::string nlfi, const std::vector<mfem::Vector>& u) {
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
 public:
  template <typename... Args>
  SteadyDiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                          const std::vector<std::string>& integrators, const Parameters& params,
                          Args&&... args)
      : DiffusionOperatorBase<T, DIM, SteadyOperatorBase>(spatials, integrators, params,
                                                          std::forward<Args>(args)...) {}
  template <typename... Args>
  SteadyDiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                          const std::vector<std::string>& integrators, Args&&... args)
      : DiffusionOperatorBase<T, DIM, SteadyOperatorBase>(spatials, integrators,
                                                          std::forward<Args>(args)...) {}
};

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
 public:
  template <typename... Args>
  DiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& integrators, const Parameters& params,
                    Args&&... args)
      : DiffusionOperatorBase<T, DIM, TransientOperatorBase>(spatials, integrators, params,
                                                             std::forward<Args>(args)...) {}
  template <typename... Args>
  DiffusionOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& integrators, Args&&... args)
      : DiffusionOperatorBase<T, DIM, TransientOperatorBase>(spatials, integrators,
                                                             std::forward<Args>(args)...) {}
};
