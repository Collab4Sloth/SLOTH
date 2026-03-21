/**
 * @file AllenCahnNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the Allen Cahn equation
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

#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AllenCahnNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for AllenCahn-type
 * problems. It forwards the provided previous solution fields, simulation
 * parameters, auxiliary variables, and coefficients to the base
 * SLOTH nonlinear form integrator.
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
AllenCahnNLFormIntegrator<VARS>::AllenCahnNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "AllenCahn";

  this->check_variables_consistency();
}

/**
 * @brief Initialize the AllenCahn integrator.
 *
 * This method performs all necessary setup steps for a AllenCahn-type
 * nonlinear form integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
 * 3. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' is populated with the
 *      required coefficient types for this AllenCahn integrator.
 *
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void AllenCahnNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "AllenCahnNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      temperature_found = true;
      break;
    }
  }
  if (this->params_.has_parameter("ScaleMobilityByTemperature")) {
    this->scale_mobility_by_temperature_ =
        this->params_.template get_param_value<bool>("ScaleMobilityByTemperature");
    if (this->scale_mobility_by_temperature_) {
      MFEM_VERIFY(
          temperature_found,
          "AllenCahnNLFormIntegrator: "
          "Temperature variable required to scale mobility, but not found in auxiliary variables");
    }
  }
}

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector of the AllenCahn-type nonlinear problem by element
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
void AllenCahnNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);
  std::vector<double> vaux_gf_at_ip(this->vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Catch_Time_Section("AllenCahnNLFormIntegrator:AssembleElementVector");
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

      // Get aux values at ip TODO(cci) (move in method)
      for (size_t k = 0; k < vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
      // given u (elfun), compute grad(u)
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      const double coef_mobi =
          this->compute_coefficient(mobility[blk], std::span<const double>(u_values),
                                    std::span<const double>(vaux_gf_at_ip)) *
          ip.weight * Tr.Weight();
      const double lamb = this->compute_coefficient(lambda[blk], std::span<const double>(u_values),
                                                    std::span<const double>(vaux_gf_at_ip));
      gradU *= coef_mobi * lamb;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // Given u, compute (w'(u), psi), psi is shape function
      const double ww =
          coef_mobi * this->compute_gradient_coefficient(double_well_energy[blk], blk,
                                                         std::span<const double>(u_values),
                                                         std::span<const double>(vaux_gf_at_ip));
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix of a AllenCahn-type nonlinear problem by element.
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
void AllenCahnNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("AllenCahnNLFormIntegrator::AssembleElementGrad");
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    // elmat.SetSize(nd);
    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    // Diag
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      Tr.SetIntPoint(&ip);

      // Get aux values at ip TODO(cci) (move in method)
      for (size_t k = 0; k < vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      // Laplacian : compute (grad(u), grad(psi)), psi is shape function.
      const double coef_mobi =
          this->compute_coefficient(mobility[blk], std::span<const double>(u_values),
                                    std::span<const double>(vaux_gf_at_ip)) *
          ip.weight * Tr.Weight();
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      const double lamb = this->compute_coefficient(lambda[blk], std::span<const double>(u_values),
                                                    std::span<const double>(vaux_gf_at_ip));
      AddMult_a_AAt(coef_mobi * lamb, gradPsi, *elmats(blk, blk));

      // Compute w'(u)*(du,psi), psi is shape function ( // w''(u))
      double fun_val =
          coef_mobi * this->compute_hessian_coefficient(double_well_energy[blk], blk, blk,
                                                        std::span<const double>(u_values),
                                                        std::span<const double>(vaux_gf_at_ip));

      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));  // w'(u)*(du, psi)
    }

    // Loop splitted in two blocks, more efficient than a if condition?
    // Off-diag before blk
    for (int jblk = 0; jblk < blk; ++jblk) {
      elmats(blk, jblk)->SetSize(nd);
      *elmats(blk, jblk) = 0.0;
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const mfem::IntegrationPoint& ip = ir->IntPoint(i);
        el[blk]->CalcShape(ip, Psi);
        Tr.SetIntPoint(&ip);

        // Get aux values at ip TODO(cci) (move in method)
        for (size_t k = 0; k < vaux_gf_.size(); ++k) {
          vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
        }
        // Get values
        for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
          u_values[off_blk] = (*elfun[off_blk]) * Psi;
          u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
        }

        const double coef_mobi =
            this->compute_coefficient(mobility[jblk], std::span<const double>(u_values),
                                      std::span<const double>(vaux_gf_at_ip)) *
            ip.weight * Tr.Weight();

        double fun_val =
            coef_mobi * this->compute_hessian_coefficient(double_well_energy[jblk], blk, jblk,
                                                          std::span<const double>(u_values),
                                                          std::span<const double>(vaux_gf_at_ip));
        AddMult_a_VVt(fun_val, Psi, *elmats(blk, jblk));
      }
    }
    // Off-diag after blk
    for (int jblk = blk + 1; jblk < num_blocks; ++jblk) {
      elmats(blk, jblk)->SetSize(nd);
      *elmats(blk, jblk) = 0.0;
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const mfem::IntegrationPoint& ip = ir->IntPoint(i);
        el[blk]->CalcShape(ip, Psi);
        Tr.SetIntPoint(&ip);

        // Get aux values at ip TODO(cci) (move in method)
        for (size_t k = 0; k < vaux_gf_.size(); ++k) {
          vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
        }
        // Get values
        for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
          u_values[off_blk] = (*elfun[off_blk]) * Psi;
          u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
        }

        const double coef_mobi =
            this->compute_coefficient(mobility[jblk], std::span<const double>(u_values),
                                      std::span<const double>(vaux_gf_at_ip)) *
            ip.weight * Tr.Weight();

        double fun_val =
            coef_mobi * this->compute_hessian_coefficient(double_well_energy[jblk], blk, jblk,
                                                          std::span<const double>(u_values),
                                                          std::span<const double>(vaux_gf_at_ip));
        AddMult_a_VVt(fun_val, Psi, *elmats(blk, jblk));
      }
    }
  }
}

/**
 * @brief Retrieve and store the coefficients required by the Allen-Cahn integrator.
 *
 * This method collects the coefficients of type `Capillary`, `Mobility`, and `FreeEnergy`
 * from each coefficients and adds them to the corresponding internal storage:
 * - 'lambda' stores the capillary coefficients,
 * - 'mobility' stores the mobility coefficients,
 * - 'double_well_energy' stores the free energy coefficients.
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
void AllenCahnNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Capillary, 0).has_value()) {
      lambda.add(*(this->get_coefficient(i, GlossaryType::Capillary, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::FreeEnergy, 0).has_value()) {
      double_well_energy.add(*(this->get_coefficient(i, GlossaryType::FreeEnergy, 0)));
    }
  }
}

/**
 * @brief Return the value of the coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_

 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).

 * @return The computed scalar value of the coefficient.
 */
template <class VARS>
double AllenCahnNLFormIntegrator<VARS>::compute_coefficient(
    Coefficient coef, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  if (coef.is_scalar()) {
    return coef.compute();
  } else {
    std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
    std::span<const double> un(values.begin() + this->nb_blk_, values.end());

    const std::span<const double>& input = coef.is_implicit() ? u : un;
    return coef.compute(input, aux_values);
  }
}

/**
 * @brief Compute the value of a specific component of the gradient of a coefficient.
 *
 * This method evaluates the gradient of the given coefficient with respect to the
 * variable corresponding to a specified index 'blk'. By default, the 'values' vector
 * contains {u, u_old}, and auxiliary variables remain accessible via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose gradient is to be computed.
 * @param blk    Index of the gradient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the gradient component of the coefficient.
 *
 * @note This function can be overridden in derived classes.
 */
template <class VARS>
double AllenCahnNLFormIntegrator<VARS>::compute_gradient_coefficient(
    Coefficient coef, const int blk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
  std::span<const double> un(values.begin() + this->nb_blk_, values.end());

  const std::span<const double>& input = coef.is_implicit() ? u : un;
  return coef.compute_gradient(blk, input, aux_values);
}

/**
 * @brief Compute the value of a specific component of the Hessian of a  coefficient.
 *
 * This method evaluates the Hessian of the given coefficient with respect to the
 * variables corresponding to the specified indexes 'iblk' and 'jblk'. By default,
 * the 'values' vector contains {u, u_old}, and auxiliary variables remain accessible
 * via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose Hessian is to be computed.
 * @param iblk   Index of the row  of the Hessian component.
 * @param jblk   Index of the column  of the Hessian component.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the Hessian component for the given coefficient.
 *
 * @note This function can be overridden in derived classes to implement
 *       more complex or nonlinear Hessian evaluations.
 */
template <class VARS>
double AllenCahnNLFormIntegrator<VARS>::compute_hessian_coefficient(
    Coefficient coef, const int iblk, const int jblk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_hessian(iblk, jblk, u, aux_values);
  }
  return coef_value;
}
