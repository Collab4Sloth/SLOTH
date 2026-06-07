/**
 * @file LatentHeatNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Latent heat contribution for heat transfer equation
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
#include <list>
#include <memory>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new LatentHeatNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base SLOTH
 * nonlinear form integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 *
 */
template <class VARS>
LatentHeatNLFormIntegrator<VARS>::LatentHeatNLFormIntegrator(
    Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "LatentHeat";
  this->latent_time_step_ = this->params_.template get_param_value<double>("latent_time_step");
  this->check_variables_consistency();
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void LatentHeatNLFormIntegrator<VARS>::check_variables_consistency() {
  // Phase-field variable (AC)
  this->phase_field_index_ = -1;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "LatentHeatNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "PHI") {
      this->phase_field_index_ = i;
      break;
    }
  }
}

/**
 * @brief Initialize the integrator.
 *
 * This method performs all necessary setup steps for the integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
 * 3. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' is populated with the
 *      required coefficient types for this integrator.
 *
 */
template <class VARS>
void LatentHeatNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
    this->vaux_old_gf_.emplace_back(std::move(this->aux_old_gf_[i]));
  }
}

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector  by element
 *
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elvect  Array of vectors where the computed element residual contributions
 *                will be stored.
 *
 * @note Users typically do not call this function directly; it is invoked
 *       internally during the assembly of the global nonlinear form.
 */
template <class VARS>
void LatentHeatNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    [[maybe_unused]] const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);
  std::vector<double> vaux_gf_at_ip(this->vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Catch_Time_Section("LatentHeatNLFormIntegrator:AssembleElementVector");
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    // elvect.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      for (size_t k = 0; k < this->vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = this->vaux_gf_[k].GetValue(Tr, ip);
        vaux_gf_at_ip[k + num_blocks] = this->vaux_old_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= mfem::CylindricalRadialCoefficient().Eval(Tr, ip);
      }

      const double latent_heat =
          this->get_latent_heat_at_ip(blk, std::span<const double>(u_values),
                                      std::span<const double>(vaux_gf_at_ip)) *
          weight_coef;
      add(*elvect[blk], latent_heat, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix by element.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elmat   Array of dense matrices where the computed element Jacobian contributions
 *                will be stored.
 *
 * @note Users typically do not call this function directly; it is invoked
 *       internally during the assembly of the global nonlinear form.
 */
template <class VARS>
void LatentHeatNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el,
    [[maybe_unused]] mfem::ElementTransformation& Tr,
    [[maybe_unused]] const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;
  }
}

/**
 * @brief Retrieve and store the coefficients required by the integrator.
 *
 * This method collects the coefficients of type  `Mobility`, and `PhaseFieldPotential`
 * from each coefficients and adds them to the corresponding internal storage:
 * - 'mobility' stores the mobility coefficients,
 * - 'interpolation_potential' stores the PhaseFieldPotential coefficients.
 *
 * Only coefficients with ID 0 are considered for each block.
 *
 * @remark This method can be overridden in derived classes to provide
 *         custom behavior for retrieving coefficients.
 *
 * @tparam VARS Template parameter defining the variabls used in the integrator.
 *
 */
template <class VARS>
void LatentHeatNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
  }
}

/**
 * @brief Compute the value of latent heat at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed mobility at integration point
 */
template <class VARS>
double LatentHeatNLFormIntegrator<VARS>::get_latent_heat_at_ip(
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  std::span<const double> local_auxvalues(aux_values.begin(), aux_values.begin() + this->nb_blk_);
  std::span<const double> local_auxvalues_n(aux_values.begin() + this->nb_blk_, aux_values.end());
  const double mobility_value = this->get_mob_at_ip(blk, values, local_auxvalues);

  const double phi = local_auxvalues[this->phase_field_index_];
  const double phin = local_auxvalues_n[this->phase_field_index_];

  double square_time_derivative = (phi - phin) / this->latent_time_step_;
  square_time_derivative *= square_time_derivative;
  const double epsilon = 1.e-10;
  const double latent_heat = -square_time_derivative / std::max(epsilon, mobility_value);

  return latent_heat;
}

/**
 * @brief Compute the value of the mobility at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed mobility at integration point
 */
template <class VARS>
double LatentHeatNLFormIntegrator<VARS>::get_mob_at_ip(
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  const double mobi = this->compute_coefficient(mobility[blk], values, aux_values);

  return mobi;
}
