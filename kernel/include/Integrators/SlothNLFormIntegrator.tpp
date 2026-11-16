/**
 * @file SlothNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class prodiving common methods to all NLFormIntegrators
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * @anchor SlothNLFormIntegrator
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
#include <cstddef>
#include <list>
#include <optional>
#include <span>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Coefficients/Coefficient.hpp"
#include "Coefficients/Coefficients.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Options/ProblemsOptions.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new SlothNLFormIntegrator::SlothNLFormIntegrator object
 *
 * @tparam VARS
 * @param params
 * @param auxvars
 */
template <class VARS>
SlothNLFormIntegrator<VARS>::SlothNLFormIntegrator(Geometry geometry, const double time_step,
                                                   const std::vector<mfem::ParGridFunction> u_old,
                                                   const std::vector<mfem::ParGridFunction> aux_old,
                                                   const Parameters& params,
                                                   std::vector<VARS*> auxvars,
                                                   const std::vector<Coefficients>& coefficients)
    : geometry_(geometry),
      time_step_(time_step),
      u_old_(u_old),
      aux_old_gf_(aux_old),
      params_(params),
      coefficients_(coefficients) {
  this->nb_blk_ = this->u_old_.size();

  this->manage_auxiliary_variables(auxvars);

  this->aux_gf_ = this->get_aux_gf();
  this->aux_infos_ = this->get_aux_infos();
}

/**
 * @brief Return the vector of grid functions associated with the auxiliary variables
 *
 * @tparam VARS
 * @return std::vector<mfem::ParGridFunction>
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::manage_auxiliary_variables(std::vector<VARS*> auxvars) {
  this->auxvariables_ = auxvars;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      this->vect_aux_gf_.emplace_back(std::move(auxvar.get_gf()));

      std::vector<std::string> var_info = auxvar.get_additional_variable_info();

      this->vect_aux_infos_.emplace_back(std::move(var_info));
    }
  }

  this->nb_vaux_ = this->vect_aux_gf_.size();
}

/**
 * @brief Return a vector of GridFunction associated with auxiliary variables
 * @remark Order of the vector is implicitly the same as the order of auxiliary variables
 *
 * @tparam VARS
 * @return std::vector<mfem::ParGridFunction>
 */
template <class VARS>
std::vector<mfem::ParGridFunction> SlothNLFormIntegrator<VARS>::get_aux_gf() {
  return vect_aux_gf_;
}

/**
 * @brief Return a vector of the additional information (vector of string) associated with auxiliary
 * variables.
 * @remark Order of the vector is implicitly the same as the order of auxiliary variables
 * @remark The last element of the vector is the name of the variable
 *
 * @tparam VARS
 * @return std::vector<std::vector<std::string>>
 */
template <class VARS>
std::vector<std::vector<std::string>> SlothNLFormIntegrator<VARS>::get_aux_infos() {
  return this->vect_aux_infos_;
}

/**
 * @brief Verify that the coefficients associated with the integrator match the expected glossary
 *        types.
 *
 * The check is order-independent and ignores duplicates. If any required glossary type is missing,
 * the function aborts execution via MFEM_VERIFY.
 *
 * @tparam VARS Template parameter defining the variable set handled by the integrator.
 *
 * @param expected_types List of glossary types required by the current SLOTH-based integrator.
 *
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::check_coefficient_types(std::list<GlossaryType> expected_types) {
  for (const auto& coefficients : this->coefficients_) {
    auto vect_types = coefficients.get_types();
    std::list<GlossaryType> TestedGlossaryType;
    TestedGlossaryType.assign(vect_types.begin(), vect_types.end());
    expected_types.sort();
    TestedGlossaryType.sort();

    bool expected_types_found = std::ranges::includes(TestedGlossaryType, expected_types);

    MFEM_VERIFY(expected_types_found,
                "Error in " + this->integrator_name_ +
                    " at least one coefficient does not match with the expected list of "
                    "GlossaryType for the current SLOTH-based Integrator. Please check your data.");
  }
}

/**
 * @brief Verify that the coefficients associated with the integrator contain the expected custom
 *        quantities (see CustomQuantity()).
 *
 * @details To be called by custom integrators, as check_coefficient_types() for the built-in
 *          quantities. Aborts via MFEM_VERIFY if a name is missing in at least one block.
 *
 * @tparam VARS Template parameter defining the variable set handled by the integrator.
 * @param expected_names Names of the custom quantities required by the integrator.
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::check_custom_coefficient_names(
    const std::list<std::string>& expected_names) {
  for (const auto& coefficients : this->coefficients_) {
    for (const auto& name : expected_names) {
      bool found = false;
      for (unsigned int i = 0; i < coefficients.size() && !found; i++) {
        found = coefficients[i].get_type() == GlossaryType::Custom &&
                coefficients[i].get_custom_name() == name;
      }
      MFEM_VERIFY(found, "Error in " + this->integrator_name_ + ": custom quantity '" + name +
                             "' not found in the coefficients. Please check your data.");
    }
  }
}

/**
 * @brief Retrieve a coefficient of a custom quantity by name and identifier.
 *
 * @tparam VARS Template parameter defining the variable set handled by the integrator.
 * @param blk          Index of the block.
 * @param custom_name  Name of the custom quantity (see CustomQuantity()).
 * @param id           Identifier of the coefficient.
 * @param bdr_id       Optional boundary id for boundary coefficients.
 * @return An std::optional containing the matching Coefficient if found; std::nullopt otherwise.
 */
template <class VARS>
std::optional<Coefficient> SlothNLFormIntegrator<VARS>::get_coefficient(
    const int blk, const std::string& custom_name, unsigned int id, std::optional<int> bdr_id) {
  const Coefficients& coefficients = this->coefficients_[blk];
  for (unsigned int i = 0; i < coefficients.size(); i++) {
    const auto& coef = coefficients[i];
    if (coef.get_type() != GlossaryType::Custom || coef.get_custom_name() != custom_name ||
        coef.get_id() != id) {
      continue;
    }
    if (!bdr_id.has_value()) {
      return coef;
    }
    const auto bdr_index = coef.get_bdr_index_coef();
    if (std::find(bdr_index.begin(), bdr_index.end(), bdr_id) != bdr_index.end()) {
      return coef;
    }
  }
  return std::nullopt;
}

/**
 * @brief Retrieve a coefficient by type and identifier.
 *
 * This function searches the list of coefficients associated with the
 * specified block (blk) and returns the coefficient that matches
 * the given glossary type, identifier and if given, the boundary id.
 *
 * @tparam VARS Template parameter defining the variables.
 *
 * @param blk     Index of the block.
 * @param type    Type of coefficient.
 * @param id      Identifier of the coefficient.
 * @param bdr_id  Optional boundary id for boundary coefficients.
 *
 * @return An std::optional containing the matching Coefficient if found;
 *         std::nullopt otherwise.
 */
template <class VARS>
std::optional<Coefficient> SlothNLFormIntegrator<VARS>::get_coefficient(const int blk,
                                                                        GlossaryType type,
                                                                        unsigned int id,
                                                                        std::optional<int> bdr_id) {
  const Coefficients& coefficients = this->coefficients_[blk];

  for (unsigned int i = 0; i < coefficients.size(); i++) {
    const auto& coef = coefficients[i];
    if (coef.get_type() == type && coef.get_id() == id) {
      if (bdr_id.has_value()) {
        auto bdr_index = coef.get_bdr_index_coef();
        bool bdr_id_found =
            std::find(bdr_index.begin(), bdr_index.end(), bdr_id) != bdr_index.end();
        if (bdr_id_found) return coef;
      } else {
        return coef;
      }
    }
  }
  return std::nullopt;
}

/**
 * @brief Set the current integration point.
 *
 * Wraps Tr.SetIntPoint(&ip) and keeps a pointer to the transformation so that
 * file-based coefficients can retrieve the physical coordinates of the point.
 * The physical coordinates are computed lazily (only if requested).
 *
 * @tparam VARS Template parameter defining the variables.
 * @param Tr Element transformation.
 * @param ip Integration point (reference coordinates).
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::set_int_point(mfem::ElementTransformation& Tr,
                                                const mfem::IntegrationPoint& ip) {
  Tr.SetIntPoint(&ip);
  this->current_tr_ = &Tr;
  this->x_ip_valid_ = false;
}

/**
 * @brief Return the physical coordinates of the current integration point.
 *
 * Computed once per integration point, then cached until the next call to set_int_point().
 *
 * @tparam VARS Template parameter defining the variables.
 * @return Physical coordinates of the current integration point.
 */
template <class VARS>
std::span<const double> SlothNLFormIntegrator<VARS>::get_physical_point() {
  MFEM_VERIFY(this->current_tr_ != nullptr,
              "Error in " + this->integrator_name_ +
                  ": no integration point set. Use set_int_point(Tr, ip) instead of "
                  "Tr.SetIntPoint(&ip).");
  if (!this->x_ip_valid_) {
    this->current_tr_->Transform(this->current_tr_->GetIntPoint(), this->x_ip_);
    this->x_ip_valid_ = true;
  }
  return std::span<const double>(this->x_ip_.GetData(), this->x_ip_.Size());
}

//////////////////////////////////////////////////////////////////////////////
// Coefficient evaluation (by reference: no copy of the Coefficient object)
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Return the value of the coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient (may be modified: coordinates are set before evaluation).
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the coefficient.
 */
template <class VARS>
double SlothNLFormIntegrator<VARS>::eval_coefficient(Coefficient& coef,
                                                     const std::span<const double>& values,
                                                     const std::span<const double>& aux_values,
                                                     const int var_blk) {
  if (coef.uses_coordinates()) {
    coef.set_coordinates(this->get_physical_point());
  }
  if (coef.is_file()) {
    return coef.compute(this->get_physical_point(), std::span<const double>());
  }
  if (coef.is_scalar()) {
    return coef.compute();
  } else {
    const auto [u, un] = this->coefficient_inputs(coef, values, var_blk);
    if (coef.is_semi_implicit()) {
      return coef.compute(u, un, aux_values);
    } else {
      return coef.compute(coef.is_implicit() ? u : un, aux_values);
    }
  }
}

/**
 * @brief Compute the value of a specific component of the gradient of a coefficient.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose gradient is to be computed.
 * @param blk    Index of the gradient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the gradient component of the coefficient.
 */
template <class VARS>
double SlothNLFormIntegrator<VARS>::eval_gradient_coefficient(
    Coefficient& coef, const int blk, const std::span<const double>& values,
    const std::span<const double>& aux_values, const int var_blk) {
  if (coef.uses_coordinates()) {
    coef.set_coordinates(this->get_physical_point());
  }
  if (coef.is_scalar() || coef.is_explicit() || coef.is_file()) {
    return 0.0;
  }
  // Coefficient of one variable: derivative with respect to its own block only
  const bool single = coef.is_single_variable();
  if (single && blk != var_blk) return 0.0;
  const int i = single ? 0 : blk;
  const auto [u, un] = this->coefficient_inputs(coef, values, var_blk);
  if (coef.is_semi_implicit()) {
    return coef.compute_gradient(i, u, un, aux_values);
  } else {
    return coef.compute_gradient(i, coef.is_implicit() ? u : un, aux_values);
  }
}

/**
 * @brief Compute the value of a specific component of the gradient of an energy coefficient.
 *
 * This method must be applied when the time discretization concerns the first derivative of the
 * coefficient (as the free energy and the interpolation function). In that case, when time scheme
 * is explicit, first derivative are not equal to zero.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose gradient is to be computed.
 * @param blk    Index of the gradient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the gradient component of the coefficient.
 */
template <class VARS>
double SlothNLFormIntegrator<VARS>::eval_gradient_energy_coefficient(
    Coefficient& coef, const int blk, const std::span<const double>& values,
    const std::span<const double>& aux_values, const int var_blk) {
  if (coef.uses_coordinates()) {
    coef.set_coordinates(this->get_physical_point());
  }
  if (coef.is_scalar() || coef.is_file()) {
    return 0.0;
  }
  // Coefficient of one variable: derivative with respect to its own block only
  const bool single = coef.is_single_variable();
  if (single && blk != var_blk) return 0.0;
  const int i = single ? 0 : blk;
  const auto [u, un] = this->coefficient_inputs(coef, values, var_blk);
  if (coef.is_semi_implicit()) {
    return coef.compute_gradient(i, u, un, aux_values);
  } else {
    return coef.compute_gradient(i, coef.is_implicit() ? u : un, aux_values);
  }
}

/**
 * @brief Compute the value of a specific component of the Hessian of a  coefficient.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose Hessian is to be computed.
 * @param iblk   Index of the row  of the Hessian component.
 * @param jblk   Index of the column of the Hessian component.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the Hessian component for the given coefficient.
 */
template <class VARS>
double SlothNLFormIntegrator<VARS>::eval_hessian_coefficient(
    Coefficient& coef, const int iblk, const int jblk, const std::span<const double>& values,
    const std::span<const double>& aux_values, const int var_blk) {
  if (coef.uses_coordinates()) {
    coef.set_coordinates(this->get_physical_point());
  }
  if (coef.is_scalar() || coef.is_explicit() || coef.is_file()) {
    return 0.0;
  }
  // Coefficient of one variable: only the second derivative with respect to its own block
  const bool single = coef.is_single_variable();
  if (single && (iblk != var_blk || jblk != var_blk)) return 0.0;
  const int i = single ? 0 : iblk;
  const int j = single ? 0 : jblk;
  const auto [u, un] = this->coefficient_inputs(coef, values, var_blk);
  if (coef.is_semi_implicit()) {
    return coef.compute_hessian(i, j, u, un, aux_values);
  } else {
    return coef.compute_hessian(i, j, u, aux_values);
  }
}

/**
 * @brief {u_old + dt k, u_old} from {k, u_old}, for the integrators of the time derivative
 * (see the declaration).
 */
template <class VARS>
std::span<const double> SlothNLFormIntegrator<VARS>::values_at_end_of_step(
    const std::span<const double>& values) {
  const std::size_t nb = static_cast<std::size_t>(this->nb_blk_);
  MFEM_VERIFY(values.size() == 2 * nb, "Error in " + this->integrator_name_ +
                                           ": {k, u_old} expected at the integration point.");
  this->u_end_.resize(2 * nb);
  for (std::size_t b = 0; b < nb; ++b) {
    this->u_end_[b] = values[nb + b] + this->time_step_ * values[b];
    this->u_end_[nb + b] = values[nb + b];
  }
  return std::span<const double>(this->u_end_);
}

/**
 * @brief Values given to a coefficient at the current point: all the blocks {u, u_old}, or
 * only the values of block var_blk for a coefficient of one variable.
 *
 * @param values {u, u_old} at the current point (u_old may be absent)
 */
template <class VARS>
std::pair<std::span<const double>, std::span<const double>>
SlothNLFormIntegrator<VARS>::coefficient_inputs(const Coefficient& coef,
                                                const std::span<const double>& values,
                                                const int var_blk) const {
  const std::size_t nb = std::min(static_cast<std::size_t>(this->nb_blk_), values.size());
  std::span<const double> u = values.first(nb);
  std::span<const double> un = values.subspan(nb);
  if (coef.is_single_variable()) {
    const auto k = static_cast<std::size_t>(var_blk);
    MFEM_VERIFY(k < u.size(), "Error in " + this->integrator_name_ +
                                  ": coefficient of one variable evaluated for block " +
                                  std::to_string(var_blk) + ", which does not exist.");
    u = u.subspan(k, 1);
    un = k < un.size() ? un.subspan(k, 1) : std::span<const double>();
  }
  return {u, un};
}

//////////////////////////////////////////////////////////////////////////////
// Evaluation of {u, u_old} and of the auxiliary fields at integration points
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief True if two finite elements are interchangeable for the evaluation of a field
 *        (same type, geometry, order, number of dofs and map/range types).
 *
 * @tparam VARS Template parameter defining the variables.
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::same_finite_element(const mfem::FiniteElement& a,
                                                      const mfem::FiniteElement& b) {
  if (&a == &b) return true;
  return typeid(a) == typeid(b) && a.GetGeomType() == b.GetGeomType() &&
         a.GetOrder() == b.GetOrder() && a.GetDof() == b.GetDof() &&
         a.GetMapType() == b.GetMapType() && a.GetRangeType() == b.GetRangeType();
}

/**
 * @brief True if all the blocks share the same finite element.
 *
 * @tparam VARS Template parameter defining the variables.
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::same_finite_elements(
    const mfem::Array<const mfem::FiniteElement*>& el) {
  for (int b = 1; b < el.Size(); ++b) {
    if (!same_finite_element(*el[0], *el[b])) return false;
  }
  return true;
}

/**
 * @brief Prepare the evaluation of {u, u_old} and of the auxiliary fields on the current element.
 *
 * Must be called at the beginning of each AssembleElementVector/AssembleElementGrad that uses
 * values_at_ip().
 *
 * When all the blocks share the same finite element, the values at every point of `ir` are
 * computed here once for all the blocks: the element dofs of u_old and of the auxiliary fields
 * are extracted once per element, and the value at a point is the dot product with the shape
 * functions, i.e. exactly what GridFunction::GetValue(Tr, ip) computes. An auxiliary field whose
 * finite element differs from the one of the blocks is evaluated with GetValue.
 *
 * Otherwise nothing is cached and values_at_ip() evaluates the fields at each call.
 *
 * @tparam VARS Template parameter defining the variables.
 * @param el     Finite elements of the blocks.
 * @param Tr     Element transformation (its integration point is modified).
 * @param elfun  Current element dofs of the blocks.
 * @param aux_gf Auxiliary fields, in the order expected by the coefficients.
 * @param ir     Integration rule used by the caller.
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::prepare_element_values(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const std::vector<mfem::ParGridFunction>& aux_gf,
    const mfem::IntegrationRule& ir) {
  this->ipv_ir_ = nullptr;
  this->ipv_element_ = Tr.ElementNo;

  const int nb = el.Size();
  const std::size_t unb = static_cast<std::size_t>(nb);
  const std::size_t naux = aux_gf.size();
  if (!same_finite_elements(el) || this->u_old_.size() < unb) return;

  const mfem::FiniteElement& fe = *el[0];
  // Shape functions in reference space (GetValue uses CalcPhysShape otherwise)
  if (fe.GetMapType() != mfem::FiniteElement::VALUE) return;
  const int elem = Tr.ElementNo;

  // Element dofs of u_old (u_old lives on the finite element space of its block)
  this->ipv_uold_dofs_.resize(unb);
  for (std::size_t b = 0; b < unb; ++b) {
    const mfem::FiniteElementSpace* fes = this->u_old_[b].FESpace();
    if (elem >= fes->GetNE() || !same_finite_element(*fes->GetFE(elem), fe)) return;
    this->u_old_[b].GetElementDofValues(elem, this->ipv_uold_dofs_[b]);
  }

  // Element dofs of the auxiliary fields sharing the finite element of the blocks
  this->ipv_aux_dofs_.resize(naux);
  this->ipv_aux_dot_.assign(naux, 0);
  bool any_getvalue = false;
  for (std::size_t k = 0; k < naux; ++k) {
    const mfem::FiniteElementSpace* fes = aux_gf[k].FESpace();
    if (elem < fes->GetNE() && same_finite_element(*fes->GetFE(elem), fe)) {
      aux_gf[k].GetElementDofValues(elem, this->ipv_aux_dofs_[k]);
      this->ipv_aux_dot_[k] = 1;
    } else {
      any_getvalue = true;
    }
  }

  // Values at every integration point
  const int nip = ir.GetNPoints();
  this->ipv_shape_.SetSize(fe.GetDof());
  this->ipv_u_table_.resize(static_cast<std::size_t>(nip) * 2 * unb);
  this->ipv_aux_table_.resize(static_cast<std::size_t>(nip) * naux);
  for (int i = 0; i < nip; ++i) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(i);
    fe.CalcShape(ip, this->ipv_shape_);

    double* u = this->ipv_u_table_.data() + static_cast<std::size_t>(i) * 2 * unb;
    for (std::size_t b = 0; b < unb; ++b) {
      u[b] = (*elfun[static_cast<int>(b)]) * this->ipv_shape_;
      u[b + unb] = this->ipv_uold_dofs_[b] * this->ipv_shape_;
    }

    if (any_getvalue) this->set_int_point(Tr, ip);
    double* aux = this->ipv_aux_table_.data() + static_cast<std::size_t>(i) * naux;
    for (std::size_t k = 0; k < naux; ++k) {
      aux[k] = this->ipv_aux_dot_[k] ? this->ipv_aux_dofs_[k] * this->ipv_shape_
                                     : aux_gf[k].GetValue(Tr, ip);
    }
  }

  this->ipv_ir_ = &ir;
  this->ipv_nb_ = nb;
  this->ipv_naux_ = naux;
}

/**
 * @brief {u, u_old} and auxiliary values at the integration point `i` of `ir`.
 *
 * Returns the values computed by prepare_element_values() when available. Otherwise the values
 * are evaluated here: current values with the shape functions `Psi` of the
 * current block, u_old and auxiliary fields with GridFunction::GetValue (the integration point of
 * Tr must then be set by the caller).
 *
 * The returned spans remain valid until the next call to values_at_ip() or
 * prepare_element_values().
 *
 * @tparam VARS Template parameter defining the variables.
 */
template <class VARS>
typename SlothNLFormIntegrator<VARS>::IpValues SlothNLFormIntegrator<VARS>::values_at_ip(
    const mfem::IntegrationRule& ir, const int i, const mfem::Vector& Psi,
    const mfem::Array<const mfem::Vector*>& elfun, const std::vector<mfem::ParGridFunction>& aux_gf,
    mfem::ElementTransformation& Tr) {
  const int nb = elfun.Size();
  const std::size_t unb = static_cast<std::size_t>(nb);
  const std::size_t naux = aux_gf.size();

  if (this->ipv_ir_ == &ir && this->ipv_element_ == Tr.ElementNo && this->ipv_nb_ == nb &&
      this->ipv_naux_ == naux) {
    const std::size_t ui = static_cast<std::size_t>(i);
    return {std::span<const double>(this->ipv_u_table_.data() + ui * 2 * unb, 2 * unb),
            std::span<const double>(this->ipv_aux_table_.data() + ui * naux, naux)};
  }

  // No values cached for this element: evaluation at the point
  const mfem::IntegrationPoint& ip = ir.IntPoint(i);
  this->ipv_u_.resize(2 * unb);
  this->ipv_aux_.resize(naux);
  for (std::size_t k = 0; k < naux; ++k) {
    this->ipv_aux_[k] = aux_gf[k].GetValue(Tr, ip);
  }
  for (std::size_t b = 0; b < unb; ++b) {
    this->ipv_u_[b] = (*elfun[static_cast<int>(b)]) * Psi;
    this->ipv_u_[b + unb] = this->u_old_[b].GetValue(Tr, ip);
  }
  return {std::span<const double>(this->ipv_u_), std::span<const double>(this->ipv_aux_)};
}

/**
 * @brief Gradient of a field at the current integration point of the current element.
 *
 * Same computation as SlothGridFunction::GetGradient (element dofs times the physical gradients
 * of the shape functions), without copying the field into a SlothGridFunction.
 *
 * @tparam VARS Template parameter defining the variables.
 * @param gf      Field (defined on the mesh of the current element).
 * @param Tr      Element transformation (only the element number is used).
 * @param gradPsi Physical gradients of the shape functions at the current integration point.
 * @param grad    Gradient of the field (output).
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::element_gradient(const mfem::ParGridFunction& gf,
                                                   const mfem::ElementTransformation& Tr,
                                                   const mfem::DenseMatrix& gradPsi,
                                                   mfem::Vector& grad) {
  gf.GetElementDofValues(Tr.ElementNo, this->grad_dofs_);
  gradPsi.MultTranspose(this->grad_dofs_, grad);
}

/**
 * @brief Indicates whether the Problem uses an axisymmetric geometry.
 *
 * @tparam VAR Type representing the problem variables.
 *
 * @retval true The problem geometry is axisymmetric.
 * @retval false The problem geometry is not axisymmetric.
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::isAxisymmetric() const noexcept {
  return this->geometry_ == Geometry::Axisymmetric;
}

/**
 * @brief Set the physical time step of the current step
 *
 * @tparam VAR Type representing the problem variables.
 * @param dt the physical time
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::set_physical_time_step(const double dt) {
  this->physical_time_step_ = dt;
}

/**
 * @brief Return the physical time step of the current step
 *
 * @tparam VAR Type representing the problem variables.
 */
template <class VARS>
double SlothNLFormIntegrator<VARS>::get_physical_time_step() const {
  return this->physical_time_step_ > 0. ? this->physical_time_step_ : this->time_step_;
}
/**
 * @brief Evaluates a coefficient at all the integration points of an element.
 *
 * The coefficient is evaluated by batch with its evaluation structure when it has one
 * (FunctionCoefficient::eval_batch: coefficients generated from a JSON file), and point by point
 * otherwise (hand-written FunctionCoefficient...). A constant coefficient is evaluated once, a
 * coefficient from a file at the physical points of the element (eval_file_at_points).
 * In all cases, out holds the same values as the point by point evaluation (eval_coefficient,
 * eval_gradient_coefficient, eval_gradient_energy_coefficient, eval_hessian_coefficient):
 *  - the values of the variables are selected according to the scheme of the coefficient: u
 *    (implicit), u_old (explicit), u and u_old (semi-implicit);
 *  - a coefficient of one variable is evaluated on the block var_blk only, and its derivatives
 *    with respect to the other blocks are zero;
 *  - the second derivatives of an explicit coefficient are zero.
 *
 * @param coef Coefficient to evaluate
 * @param var_blk Block of the variable of a coefficient of one variable (see eval_coefficient)
 * @param outputs Requested outputs: CoefValue, CoefFirstDerivative, CoefSecondDerivative (with |)
 * @param fe Finite element of the blocks
 * @param Tr Transformation of the element
 * @param ir Integration rule
 * @param elfun Values of the dofs of each block on the element
 * @param aux_gf Auxiliary fields
 * @param out Values at the points
 * @param values_at Values of the variables: unknowns {u, u_old}, or end of the step
 *        {u_old + dt k, u_old} for the integrators of the time derivative
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::eval_coefficient_at_points(
    Coefficient& coef, const int var_blk, const unsigned int outputs, const mfem::FiniteElement& fe,
    mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const std::vector<mfem::ParGridFunction>& aux_gf,
    CoefficientAtPoints& out, const ValuesAt values_at) {
  const auto np = static_cast<std::size_t>(ir.GetNPoints());
  const auto nb = static_cast<std::size_t>(elfun.Size());
  const bool want_value = (outputs & CoefValue) != 0u;
  const bool want_first_derivative = (outputs & CoefFirstDerivative) != 0u;
  const bool want_second_derivative = (outputs & CoefSecondDerivative) != 0u;

  // Values prepared at all the points of the local mesh (see prepare_quadrature_points)
  if (this->read_prepared_coefficient(coef, var_blk, outputs, aux_gf, values_at, Tr.ElementNo, ir,
                                      out)) {
    return;
  }

  out.nb = elfun.Size();
  out.is_explicit = coef.is_explicit();
  out.by_batch = false;
  out.F.assign(want_value ? np : 0, 0.);
  out.dF.assign(want_first_derivative ? np * nb : 0, 0.);
  out.d2F.assign(want_second_derivative ? np * nb * nb : 0, 0.);

  // Constant coefficient: same value at all the points, zero derivatives
  if (coef.is_scalar()) {
    if (want_value) std::ranges::fill(out.F, coef.compute());
    return;
  }
  // Coefficient from a file: tabulated values at the physical points, zero derivatives
  if (coef.is_file()) {
    if (want_value) this->eval_file_at_points(coef, Tr, ir, out);
    return;
  }

  const ElementValues values = this->values_at_points(fe, Tr, ir, elfun, aux_gf, values_at);
  const bool coordinates = coef.uses_coordinates() && coef.batch_nb_variables() > 0;
  out.by_batch = this->eval_by_batch(
      coef, var_blk, outputs, ir.GetNPoints(), values,
      coordinates ? this->element_coordinates(Tr, ir) : std::span<const double>(), out);
  if (!out.by_batch) {
    this->eval_point_by_point(coef, var_blk, outputs, Tr, ir, values, out);
  }
}

/**
 * @brief Values {u, u_old} and auxiliary fields at the integration points of an element.
 *
 * Tables of prepare_element_values if they are prepared for this element and this rule (same
 * layout, no copy), computed point by point with values_at_ip otherwise. For the integrators of
 * the time derivative (ValuesAt::EndOfStep), {k, u_old} is replaced by {u_old + dt k, u_old}, as
 * in values_at_end_of_step.
 *
 * @return Views on the values: valid until the next evaluation
 */
template <class VARS>
typename SlothNLFormIntegrator<VARS>::ElementValues SlothNLFormIntegrator<VARS>::values_at_points(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const std::vector<mfem::ParGridFunction>& aux_gf,
    const ValuesAt values_at) {
  ElementValues values;
  values.nb = static_cast<std::size_t>(elfun.Size());
  values.naux = aux_gf.size();
  const int nip = ir.GetNPoints();
  const auto np = static_cast<std::size_t>(nip);
  const std::size_t width = 2 * values.nb;  // [u, u_old] at one point

  const bool table = this->ipv_ir_ == &ir && this->ipv_element_ == Tr.ElementNo &&
                     this->ipv_nb_ == elfun.Size() && this->ipv_naux_ == values.naux;
  if (table) {
    values.u_uold = std::span<const double>(this->ipv_u_table_).first(np * width);
    values.aux = std::span<const double>(this->ipv_aux_table_).first(np * values.naux);
  } else {
    this->ipv_shape_.SetSize(fe.GetDof());
    this->cap_u_.resize(np * width);
    this->cap_aux_.resize(np * values.naux);
    for (int q = 0; q < nip; ++q) {
      const mfem::IntegrationPoint& ip = ir.IntPoint(q);
      this->set_int_point(Tr, ip);
      fe.CalcShape(ip, this->ipv_shape_);
      const auto ipv = this->values_at_ip(ir, q, this->ipv_shape_, elfun, aux_gf, Tr);
      const auto uq = static_cast<std::size_t>(q);
      std::ranges::copy(ipv.u, this->cap_u_.begin() + static_cast<std::ptrdiff_t>(uq * width));
      std::ranges::copy(ipv.aux,
                        this->cap_aux_.begin() + static_cast<std::ptrdiff_t>(uq * values.naux));
    }
    values.u_uold = this->cap_u_;
    values.aux = this->cap_aux_;
  }

  if (values_at == ValuesAt::EndOfStep) {
    // [k, u_old] -> [u_old + dt k, u_old], as values_at_end_of_step
    this->cap_u_end_of_step_.resize(np * width);
    const std::size_t nb = values.nb;
    for (std::size_t q = 0; q < np; ++q) {
      const std::span<const double> k_uold = values.u_uold_at(q);
      const std::span<double> u_end_of_step =
          std::span<double>(this->cap_u_end_of_step_).subspan(q * width, width);
      for (std::size_t b = 0; b < nb; ++b) {
        u_end_of_step[b] = k_uold[nb + b] + this->time_step_ * k_uold[b];
        u_end_of_step[nb + b] = k_uold[nb + b];
      }
    }
    values.u_uold = this->cap_u_end_of_step_;
  }
  return values;
}

/**
 * @brief Evaluation by batch of a coefficient with an evaluation structure (see
 * FunctionCoefficient::eval_batch), at npts points: the integration points of an element, or all
 * the points of the local mesh (see prepare_quadrature_points).
 *
 * The values given to the coefficient are selected according to its scheme and, for a
 * coefficient of one variable, its block; the derivatives returned with respect to its variables
 * are then stored with respect to the blocks of the integrator. out.F, out.dF and out.d2F are
 * sized by the caller.
 *
 * @param npts Number of points
 * @param values {u, u_old} and auxiliary fields at the points
 * @param coordinates Physical coordinates of the points (dim per point), empty if the coefficient
 *        does not use them
 * @return false if the coefficient has no evaluation structure (nothing is computed)
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::eval_by_batch(const Coefficient& coef, const int var_blk,
                                                const unsigned int outputs, const int npts,
                                                const ElementValues& values,
                                                const std::span<const double>& coordinates,
                                                CoefficientAtPoints& out) {
  if (coef.batch_nb_variables() <= 0) {
    return false;
  }
  const bool want_value = (outputs & CoefValue) != 0u;
  const bool want_first_derivative = (outputs & CoefFirstDerivative) != 0u;
  const bool want_second_derivative = (outputs & CoefSecondDerivative) != 0u;
  const int nip = npts;
  const auto np = static_cast<std::size_t>(nip);
  const std::size_t nb = values.nb;

  const bool single = coef.is_single_variable();
  MFEM_VERIFY(!single || (var_blk >= 0 && static_cast<std::size_t>(var_blk) < nb),
              "Error in " + this->integrator_name_ +
                  ": coefficient of one variable evaluated for block " + std::to_string(var_blk) +
                  ", which does not exist.");
  // Values given to the coefficient: n per point, as in the point by point evaluation; its
  // derivatives are taken with respect to its nc first variables
  const std::size_t n = single ? 1 : nb;
  const auto nc = static_cast<std::size_t>(coef.batch_nb_variables());
  MFEM_VERIFY(nc <= n, "Error in " + this->integrator_name_ + ": coefficient of " +
                           std::to_string(nc) + " variables evaluated with " + std::to_string(n) +
                           " values.");
  const bool semi = coef.is_semi_implicit();
  const bool previous = coef.is_explicit();
  // Block of the k-th variable of the coefficient
  const auto block = [single, var_blk](const std::size_t k) {
    return single ? static_cast<std::size_t>(var_blk) : k;
  };

  // Inputs written on the host (HostWrite: the device copy, if any, is no longer valid)
  out.values.SetSize(static_cast<int>(np * n));
  out.explicit_values.SetSize(semi ? static_cast<int>(np * n) : 0);
  const std::span<double> h_values(out.values.HostWrite(), np * n);
  const std::span<double> h_explicit(semi ? out.explicit_values.HostWrite() : nullptr,
                                     semi ? np * n : 0);
  for (std::size_t q = 0; q < np; ++q) {
    const std::span<const double> u_uold = values.u_uold_at(q);
    for (std::size_t k = 0; k < n; ++k) {
      const std::size_t b = block(k);
      h_values[q * n + k] = previous ? u_uold[nb + b] : u_uold[b];
      if (semi) h_explicit[q * n + k] = u_uold[nb + b];
    }
  }
  out.auxiliary_values.SetSize(static_cast<int>(values.aux.size()));
  if (!values.aux.empty()) {
    std::ranges::copy(values.aux, out.auxiliary_values.HostWrite());
  }
  out.coordinates.SetSize(static_cast<int>(coordinates.size()));
  if (!coordinates.empty()) {
    std::ranges::copy(coordinates, out.coordinates.HostWrite());
  }

  // Second derivatives of an explicit coefficient: zero, not computed
  const bool batch_second_derivative = want_second_derivative && !previous;
  out.batch_F.SetSize(want_value ? nip : 0);
  out.batch_gradient.SetSize(want_first_derivative ? static_cast<int>(np * nc) : 0);
  out.batch_hessian.SetSize(batch_second_derivative ? static_cast<int>(np * nc * nc) : 0);
  if (!coef.eval_batch(nip, out.values, out.explicit_values, out.auxiliary_values, out.coordinates,
                       want_value ? &out.batch_F : nullptr,
                       want_first_derivative ? &out.batch_gradient : nullptr,
                       batch_second_derivative ? &out.batch_hessian : nullptr)) {
    return false;
  }

  // Outputs written where the loop ran: read on the host, stored with respect to the blocks
  if (want_value) {
    std::ranges::copy(std::span<const double>(out.batch_F.HostRead(), np), out.F.begin());
  }
  if (want_first_derivative) {
    const std::span<const double> gradient(out.batch_gradient.HostRead(), np * nc);
    for (std::size_t q = 0; q < np; ++q) {
      for (std::size_t k = 0; k < nc; ++k) {
        out.dF[q * nb + block(k)] = gradient[q * nc + k];
      }
    }
  }
  if (batch_second_derivative) {
    const std::span<const double> hessian(out.batch_hessian.HostRead(), np * nc * nc);
    for (std::size_t q = 0; q < np; ++q) {
      for (std::size_t k = 0; k < nc; ++k) {
        for (std::size_t l = 0; l < nc; ++l) {
          out.d2F[(q * nb + block(k)) * nb + block(l)] = hessian[(q * nc + k) * nc + l];
        }
      }
    }
  }
  return true;
}

/**
 * @brief Evaluation of a coefficient from a file at all the integration points of an element: one
 * call of FileCoefficient::eval_batch at the physical points, computed as in get_physical_point
 * (same values as eval_coefficient point by point).
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::eval_file_at_points(const Coefficient& coef,
                                                      mfem::ElementTransformation& Tr,
                                                      const mfem::IntegrationRule& ir,
                                                      CoefficientAtPoints& out) {
  const int nip = ir.GetNPoints();
  const int dim = Tr.GetSpaceDim();
  out.coordinates.SetSize(nip * dim);
  double* x = out.coordinates.HostWrite();
  for (int q = 0; q < nip; ++q) {
    this->set_int_point(Tr, ir.IntPoint(q));
    std::ranges::copy(this->get_physical_point(), x + static_cast<std::ptrdiff_t>(q) * dim);
  }
  out.batch_F.SetSize(nip);
  const mfem::Vector none;
  MFEM_VERIFY(
      coef.eval_batch(nip, none, none, none, out.coordinates, &out.batch_F, nullptr, nullptr),
      "Error in " + this->integrator_name_ + ": evaluation of a coefficient from a file.");
  std::ranges::copy(std::span<const double>(out.batch_F.HostRead(), static_cast<std::size_t>(nip)),
                    out.F.begin());
  out.by_batch = true;
}

/**
 * @brief Evaluation point by point of a coefficient without evaluation structure (hand-written
 * FunctionCoefficient...), at all the integration points of an element, with
 * eval_coefficient, eval_gradient_energy_coefficient and eval_hessian_coefficient.
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::eval_point_by_point(Coefficient& coef, const int var_blk,
                                                      const unsigned int outputs,
                                                      mfem::ElementTransformation& Tr,
                                                      const mfem::IntegrationRule& ir,
                                                      const ElementValues& values,
                                                      CoefficientAtPoints& out) {
  const bool want_value = (outputs & CoefValue) != 0u;
  const bool want_first_derivative = (outputs & CoefFirstDerivative) != 0u;
  const bool want_second_derivative = (outputs & CoefSecondDerivative) != 0u;
  const int nb = static_cast<int>(values.nb);

  for (int q = 0; q < ir.GetNPoints(); ++q) {
    this->set_int_point(Tr, ir.IntPoint(q));
    const auto uq = static_cast<std::size_t>(q);
    const std::span<const double> u_uold = values.u_uold_at(uq);
    const std::span<const double> aux = values.aux_at(uq);
    if (want_value) {
      out.F[uq] = this->eval_coefficient(coef, u_uold, aux, var_blk);
    }
    for (int i = 0; i < nb && want_first_derivative; ++i) {
      out.dF[uq * values.nb + static_cast<std::size_t>(i)] =
          this->eval_gradient_energy_coefficient(coef, i, u_uold, aux, var_blk);
    }
    for (int i = 0; i < nb && want_second_derivative; ++i) {
      for (int j = 0; j < nb; ++j) {
        out.d2F[(uq * values.nb + static_cast<std::size_t>(i)) * values.nb +
                static_cast<std::size_t>(j)] =
            this->eval_hessian_coefficient(coef, i, j, u_uold, aux, var_blk);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Coefficients evaluated at all the quadrature points of the local mesh
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Physical coordinates of the integration points of an element (dim per point), as
 * ElementTransformation::Transform.
 *
 * @return View valid until the next call
 */
template <class VARS>
std::span<const double> SlothNLFormIntegrator<VARS>::element_coordinates(
    mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir) {
  Tr.Transform(ir, this->cap_points_);  // dim x nip, column q = point q
  return std::span<const double>(this->cap_points_.GetData(),
                                 static_cast<std::size_t>(this->cap_points_.Height()) *
                                     static_cast<std::size_t>(ir.GetNPoints()));
}

/**
 * @brief Integration rule of the element (default: order 2 p + order of the transformation).
 *
 * Must be the rule of the assembly of the integrator: the coefficients prepared on the local mesh
 * are read during the assembly only for this rule.
 */
template <class VARS>
const mfem::IntegrationRule& SlothNLFormIntegrator<VARS>::integration_rule(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr) const {
  return mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + Tr.OrderW());
}

/**
 * @brief Declares a coefficient to evaluate at all the quadrature points of the local mesh before
 * each assembly (see prepare_quadrature_points). Called once, in init().
 *
 * During the assembly, eval_coefficient_at_points returns the prepared values of the element for
 * the same arguments (coef, var_blk, aux_gf, values_at) and outputs included in those prepared.
 * Only the coefficients with an evaluation structure (generated from a JSON file) are prepared;
 * the others are evaluated element by element during the assembly, as without declaration.
 *
 * An integrator declaring all the coefficients it evaluates can skip the values at the points of
 * the element during the assembly when they are not needed (see values_prepared_on_mesh).
 *
 * @param coef Coefficient (must remain at the same address: member of the integrator)
 * @param var_blk Block of the variable of a coefficient of one variable (see eval_coefficient)
 * @param residual_outputs Outputs used by the residual (CoefValue, CoefFirstDerivative...)
 * @param jacobian_outputs Outputs used by the Jacobian
 * @param aux_gf Auxiliary fields (must remain at the same address: member of the integrator)
 * @param values_at Values of the variables (see eval_coefficient_at_points)
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::prepare_on_mesh(Coefficient& coef, const int var_blk,
                                                  const unsigned int residual_outputs,
                                                  const unsigned int jacobian_outputs,
                                                  const std::vector<mfem::ParGridFunction>& aux_gf,
                                                  const ValuesAt values_at) {
  // Constant coefficient or coefficient from a file: no values at the points needed
  if (coef.is_scalar() || coef.is_file()) {
    return;
  }
  // No evaluation structure: evaluated with the values at the points of the element
  if (coef.batch_nb_variables() <= 0) {
    this->mesh_needs_element_values_ = true;
    return;
  }
  MeshCoefficient mesh_coefficient;
  mesh_coefficient.coef = &coef;
  mesh_coefficient.var_blk = var_blk;
  mesh_coefficient.residual_outputs = residual_outputs;
  mesh_coefficient.jacobian_outputs = jacobian_outputs;
  mesh_coefficient.aux_gf = &aux_gf;
  mesh_coefficient.values_at = values_at;
  this->mesh_coefficients_.push_back(std::move(mesh_coefficient));

  // One table of values per set of auxiliary fields and per kind of values
  const bool known = std::ranges::any_of(this->mesh_values_, [&](const MeshValues& v) {
    return v.aux_gf == &aux_gf && v.values_at == values_at;
  });
  if (!known) {
    MeshValues mesh_values;
    mesh_values.aux_gf = &aux_gf;
    mesh_values.values_at = values_at;
    this->mesh_values_.push_back(std::move(mesh_values));
  }
}

/**
 * @brief Evaluates the declared coefficients (see prepare_on_mesh) at all the quadrature points of
 * the local mesh. Called by SlothParBlockNonlinearForm before each assembly.
 *
 * {u, u_old} and the auxiliary fields are first computed at the points of the integration rule:
 *  - in the domain, by interpolation on the whole local mesh (values_by_interpolation:
 *    mfem::QuadratureInterpolator), when the mesh has one type of element; same values as during
 *    the assembly up to the rounding errors;
 *  - otherwise element by element as during the assembly (values_by_element:
 *    prepare_element_values, values_at_points), from the element dofs of x extracted as in
 *    mfem::BlockNonlinearForm.
 * Each coefficient is then evaluated at all the points in one call (eval_by_batch).
 *
 * An element whose blocks do not share the same finite element is not prepared (evaluated during
 * the assembly).
 *
 * Boundary integrator (bdr_marker not null): the boundary elements of the marked attributes are
 * prepared, as in mfem::BlockNonlinearForm::MultBlocked; {u, u_old} and the auxiliary fields are
 * evaluated point by point (values_at_ip), as during the assembly of the boundary integrators.
 *
 * @param fes Finite element spaces of the blocks
 * @param x Values at the local dofs of each block
 * @param assembly Residual or Jacobian: outputs to prepare
 * @param bdr_marker Boundary attributes of a boundary integrator, nullptr for a domain integrator
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::prepare_quadrature_points(
    const mfem::Array<mfem::FiniteElementSpace*>& fes, const mfem::BlockVector& x,
    const Assembly assembly, const mfem::Array<int>* bdr_marker) {
  Catch_Time_Section("SlothNLFormIntegrator::prepare_quadrature_points");
  const int nb = fes.Size();
  const bool boundary = bdr_marker != nullptr;

  // Outputs of a coefficient for this assembly: the values at the points (and the coordinates) are
  // computed only if used by a coefficient
  const auto outputs_of = [assembly](const MeshCoefficient& m) {
    return assembly == Assembly::Residual ? m.residual_outputs : m.jacobian_outputs;
  };
  const bool coordinates =
      std::ranges::any_of(this->mesh_coefficients_, [&outputs_of](const MeshCoefficient& m) {
        return outputs_of(m) != 0u && m.coef->uses_coordinates();
      });
  this->mesh_prepared_ = false;
  for (auto& m : this->mesh_coefficients_) {
    m.prepared_outputs = 0u;
  }
  for (auto& mesh_values : this->mesh_values_) {
    mesh_values.used = std::ranges::any_of(this->mesh_coefficients_, [&](const MeshCoefficient& m) {
      return outputs_of(m) != 0u && m.aux_gf == mesh_values.aux_gf &&
             m.values_at == mesh_values.values_at;
    });
    mesh_values.nb = static_cast<std::size_t>(nb);
    mesh_values.u_uold.clear();
    mesh_values.aux.clear();
  }
  if (std::ranges::none_of(this->mesh_values_, [](const MeshValues& v) { return v.used; })) {
    return;  // nothing to prepare for this assembly
  }

  // {u, u_old}, auxiliary fields (and coordinates) at all the points of the local mesh
  if (boundary || !this->values_by_interpolation(fes, x, coordinates)) {
    this->values_by_element(fes, x, bdr_marker, coordinates);
  }

  // Each coefficient at all the points of the local mesh, in one call
  const std::size_t np = this->mesh_offsets_.back();
  const auto unb = static_cast<std::size_t>(nb);
  for (auto& m : this->mesh_coefficients_) {
    const unsigned int outputs = outputs_of(m);
    if (outputs == 0u) continue;  // not used by this assembly
    const auto it = std::ranges::find_if(this->mesh_values_, [&m](const MeshValues& v) {
      return v.aux_gf == m.aux_gf && v.values_at == m.values_at;
    });
    ElementValues values;
    values.nb = unb;
    values.naux = m.aux_gf->size();
    values.u_uold = it->u_uold;
    values.aux = it->aux;

    CoefficientAtPoints& out = m.at_points;
    out.nb = nb;
    out.is_explicit = m.coef->is_explicit();
    out.F.assign((outputs & CoefValue) != 0u ? np : 0, 0.);
    out.dF.assign((outputs & CoefFirstDerivative) != 0u ? np * unb : 0, 0.);
    out.d2F.assign((outputs & CoefSecondDerivative) != 0u ? np * unb * unb : 0, 0.);
    const std::span<const double> points = m.coef->uses_coordinates()
                                               ? std::span<const double>(this->mesh_coordinates_)
                                               : std::span<const double>();
    out.by_batch =
        this->eval_by_batch(*m.coef, m.var_blk, outputs, static_cast<int>(np), values, points, out);
    m.prepared_outputs = out.by_batch ? outputs : 0u;
  }
  this->mesh_prepared_ = true;
}

/**
 * @brief {u, u_old}, auxiliary fields (and coordinates) at the points of the local mesh, element by
 * element, as during the assembly: element tables (prepare_element_values) in the domain, point
 * by point on the boundary. See prepare_quadrature_points.
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::values_by_element(
    const mfem::Array<mfem::FiniteElementSpace*>& fes, const mfem::BlockVector& x,
    const mfem::Array<int>* bdr_marker, const bool coordinates) {
  const int nb = fes.Size();
  const bool boundary = bdr_marker != nullptr;
  const mfem::Mesh& mesh = *fes[0]->GetMesh();
  const int ne = boundary ? fes[0]->GetNBE() : fes[0]->GetNE();

  this->mesh_rules_.assign(static_cast<std::size_t>(ne), nullptr);
  this->mesh_offsets_.assign(static_cast<std::size_t>(ne) + 1, 0);
  this->mesh_coordinates_.clear();
  this->mesh_elfun_.resize(static_cast<std::size_t>(nb));
  mfem::Array<const mfem::FiniteElement*> el(nb);
  mfem::Array<const mfem::Vector*> elfun(nb);
  mfem::DofTransformation doftrans;

  for (int e = 0; e < ne; ++e) {
    const auto ue = static_cast<std::size_t>(e);
    this->mesh_offsets_[ue + 1] = this->mesh_offsets_[ue];
    if (boundary && (*bdr_marker)[mesh.GetBdrAttribute(e) - 1] == 0) continue;
    for (int s = 0; s < nb; ++s) el[s] = boundary ? fes[s]->GetBE(e) : fes[s]->GetFE(e);
    if (!same_finite_elements(el)) continue;

    mfem::ElementTransformation& Tr =
        boundary ? *fes[0]->GetBdrElementTransformation(e) : *fes[0]->GetElementTransformation(e);
    const mfem::IntegrationRule& ir = this->integration_rule(*el[0], Tr);
    // Element dofs of each block, as in mfem::BlockNonlinearForm::MultBlocked
    for (int s = 0; s < nb; ++s) {
      mfem::Vector& elfun_s = this->mesh_elfun_[static_cast<std::size_t>(s)];
      if (boundary) {
        fes[s]->GetBdrElementVDofs(e, this->mesh_vdofs_, doftrans);
      } else {
        fes[s]->GetElementVDofs(e, this->mesh_vdofs_, doftrans);
      }
      x.GetBlock(s).GetSubVector(this->mesh_vdofs_, elfun_s);
      doftrans.InvTransformPrimal(elfun_s);
      elfun[s] = &elfun_s;
    }

    // {u, u_old} and auxiliary fields at the points, as during the assembly: element tables in
    // the domain, point by point on the boundary
    for (auto& mesh_values : this->mesh_values_) {
      if (!mesh_values.used) continue;
      if (boundary) {
        this->ipv_ir_ = nullptr;
      } else {
        this->prepare_element_values(el, Tr, elfun, *mesh_values.aux_gf, ir);
      }
      const ElementValues values =
          this->values_at_points(*el[0], Tr, ir, elfun, *mesh_values.aux_gf, mesh_values.values_at);
      mesh_values.u_uold.insert(mesh_values.u_uold.end(), values.u_uold.begin(),
                                values.u_uold.end());
      mesh_values.aux.insert(mesh_values.aux.end(), values.aux.begin(), values.aux.end());
    }
    if (coordinates) {
      const std::span<const double> points = this->element_coordinates(Tr, ir);
      this->mesh_coordinates_.insert(this->mesh_coordinates_.end(), points.begin(), points.end());
    }
    this->mesh_rules_[ue] = &ir;
    this->mesh_offsets_[ue + 1] += static_cast<std::size_t>(ir.GetNPoints());
  }
}

/**
 * @brief {u, u_old} and auxiliary fields at all the points of the local mesh (domain), by
 * interpolation of the local dofs: element restriction (local dofs element by element), then
 * mfem::QuadratureInterpolator, on the host or on the device. u_old is interpolated once (it does
 * not change during the life of the integrator). Same values as element by element (see
 * values_by_element) up to the rounding errors.
 *
 * Requires a mesh with one type of element, the same finite element for all the blocks, and
 * scalar spaces supported by mfem::QuadratureInterpolator: the same integration rule is then used
 * for all the elements.
 *
 * @return false if these conditions are not met (nothing is computed)
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::values_by_interpolation(
    const mfem::Array<mfem::FiniteElementSpace*>& fes, const mfem::BlockVector& x,
    const bool coordinates) {
  const int nb = fes.Size();
  const auto unb = static_cast<std::size_t>(nb);
  const mfem::Mesh* mesh = fes[0]->GetMesh();
  const int ne = fes[0]->GetNE();
  if (ne == 0 || mesh->GetNumGeometries(mesh->Dimension()) != 1 || this->u_old_.size() < unb) {
    return false;
  }
  const auto supported = [mesh](const mfem::FiniteElementSpace& space) {
    return space.GetMesh() == mesh && space.GetVDim() == 1 && !space.IsVariableOrder() &&
           mfem::QuadratureInterpolator::SupportsFESpace(space);
  };
  mfem::Array<const mfem::FiniteElement*> el(nb);
  for (int s = 0; s < nb; ++s) {
    if (!supported(*fes[s]) || !supported(*this->u_old_[static_cast<std::size_t>(s)].FESpace())) {
      return false;
    }
    el[s] = fes[s]->GetFE(0);
  }
  if (!same_finite_elements(el)) {
    return false;
  }
  for (const auto& mesh_values : this->mesh_values_) {
    if (!mesh_values.used) continue;
    for (const auto& gf : *mesh_values.aux_gf) {
      if (!supported(*gf.FESpace())) return false;
    }
  }

  // Same rule for all the elements (one type of element, same finite element)
  const mfem::IntegrationRule& ir =
      this->integration_rule(*el[0], *fes[0]->GetElementTransformation(0));
  const int nq = ir.GetNPoints();
  const std::size_t np = static_cast<std::size_t>(ne) * static_cast<std::size_t>(nq);

  // Local dofs -> dofs element by element -> values at the points (layout: element by element)
  const auto interpolate = [this, &ir](const mfem::FiniteElementSpace& space,
                                       const mfem::Vector& local, mfem::Vector& at_points) {
    const mfem::ElementDofOrdering ordering = mfem::UsesTensorBasis(space)
                                                  ? mfem::ElementDofOrdering::LEXICOGRAPHIC
                                                  : mfem::ElementDofOrdering::NATIVE;
    const mfem::Operator& restriction = *space.GetElementRestriction(ordering);
    this->mesh_e_vec_.SetSize(restriction.Height());
    restriction.Mult(local, this->mesh_e_vec_);
    at_points.SetSize(space.GetNE() * ir.GetNPoints());
    space.GetQuadratureInterpolator(ir)->Values(this->mesh_e_vec_, at_points);
  };

  this->mesh_q_u_.resize(unb);
  for (int s = 0; s < nb; ++s) {
    interpolate(*fes[s], x.GetBlock(s), this->mesh_q_u_[static_cast<std::size_t>(s)]);
  }
  if (this->mesh_q_uold_ir_ != &ir || this->mesh_q_uold_.size() != unb) {
    this->mesh_q_uold_.resize(unb);
    for (std::size_t b = 0; b < unb; ++b) {
      interpolate(*this->u_old_[b].FESpace(), this->u_old_[b], this->mesh_q_uold_[b]);
    }
    this->mesh_q_uold_ir_ = &ir;
  }

  // Tables of the coefficients, point by point (layout of ElementValues)
  for (auto& mesh_values : this->mesh_values_) {
    if (!mesh_values.used) continue;
    const std::size_t width = 2 * unb;
    mesh_values.u_uold.resize(np * width);
    for (std::size_t b = 0; b < unb; ++b) {
      const double* u = this->mesh_q_u_[b].HostRead();
      const double* u_old = this->mesh_q_uold_[b].HostRead();
      for (std::size_t q = 0; q < np; ++q) {
        double* u_uold = mesh_values.u_uold.data() + q * width;
        // Integrators of the time derivative: {u_old + dt k, u_old} (see values_at_end_of_step)
        u_uold[b] = mesh_values.values_at == ValuesAt::EndOfStep
                        ? u_old[q] + this->time_step_ * u[q]
                        : u[q];
        u_uold[unb + b] = u_old[q];
      }
    }
    const std::vector<mfem::ParGridFunction>& aux_gf = *mesh_values.aux_gf;
    const std::size_t naux = aux_gf.size();
    mesh_values.aux.resize(np * naux);
    this->mesh_q_aux_.resize(naux);
    for (std::size_t k = 0; k < naux; ++k) {
      interpolate(*aux_gf[k].FESpace(), aux_gf[k], this->mesh_q_aux_[k]);
      const double* aux = this->mesh_q_aux_[k].HostRead();
      for (std::size_t q = 0; q < np; ++q) {
        mesh_values.aux[q * naux + k] = aux[q];
      }
    }
  }

  this->mesh_rules_.assign(static_cast<std::size_t>(ne), &ir);
  this->mesh_offsets_.resize(static_cast<std::size_t>(ne) + 1);
  for (std::size_t e = 0; e <= static_cast<std::size_t>(ne); ++e) {
    this->mesh_offsets_[e] = e * static_cast<std::size_t>(nq);
  }
  this->mesh_coordinates_.clear();
  if (coordinates) {
    for (int e = 0; e < ne; ++e) {
      const std::span<const double> points =
          this->element_coordinates(*fes[0]->GetElementTransformation(e), ir);
      this->mesh_coordinates_.insert(this->mesh_coordinates_.end(), points.begin(), points.end());
    }
  }
  return true;
}

/**
 * @brief The values prepared by prepare_quadrature_points are no longer used: the coefficients
 * are evaluated element by element until the next preparation. Called by
 * SlothParBlockNonlinearForm at the end of each assembly.
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::clear_quadrature_points() {
  this->mesh_prepared_ = false;
  for (auto& m : this->mesh_coefficients_) {
    m.prepared_outputs = 0u;
  }
}

/**
 * @brief True if the declared coefficients (see prepare_on_mesh) are prepared for this element and
 * this rule, and none of them needs the values at the points of the element: an integrator
 * declaring all the coefficients it evaluates then does not compute these values
 * (prepare_element_values) during the assembly.
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::values_prepared_on_mesh(const int elem,
                                                          const mfem::IntegrationRule& ir) const {
  const auto ue = static_cast<std::size_t>(elem);
  return this->mesh_prepared_ && !this->mesh_needs_element_values_ && elem >= 0 &&
         ue < this->mesh_rules_.size() && this->mesh_rules_[ue] == &ir;
}

/**
 * @brief {u, u_old} and auxiliary fields at the integration points of an element, read in the
 * values prepared at all the points of the local mesh (see prepare_quadrature_points): same values
 * as values_at_ip, up to the rounding errors when they are interpolated on the whole mesh.
 *
 * Prepared when a coefficient is declared with these auxiliary fields and this kind of values
 * (see prepare_on_mesh).
 *
 * @return Views on the values of the element, valid until the next preparation; std::nullopt if
 * the values are not prepared for this element and this rule
 */
template <class VARS>
std::optional<typename SlothNLFormIntegrator<VARS>::ElementValues>
SlothNLFormIntegrator<VARS>::prepared_values(const int elem, const mfem::IntegrationRule& ir,
                                             const std::vector<mfem::ParGridFunction>& aux_gf,
                                             const ValuesAt values_at) const {
  const auto ue = static_cast<std::size_t>(elem);
  if (!this->mesh_prepared_ || elem < 0 || ue >= this->mesh_rules_.size() ||
      this->mesh_rules_[ue] != &ir) {
    return std::nullopt;
  }
  const auto it = std::ranges::find_if(this->mesh_values_, [&](const MeshValues& v) {
    return v.aux_gf == &aux_gf && v.values_at == values_at;
  });
  if (it == this->mesh_values_.end() || !it->used) {
    return std::nullopt;
  }
  ElementValues values;
  values.nb = it->nb;
  values.naux = aux_gf.size();
  const std::size_t first = this->mesh_offsets_[ue];
  const auto np = static_cast<std::size_t>(ir.GetNPoints());
  values.u_uold =
      std::span<const double>(it->u_uold).subspan(first * 2 * values.nb, np * 2 * values.nb);
  values.aux = std::span<const double>(it->aux).subspan(first * values.naux, np * values.naux);
  return values;
}

/**
 * @brief Values of a coefficient at the integration points of an element, read in the values
 * prepared at all the points of the local mesh (see prepare_quadrature_points).
 *
 * @return false if the coefficient is not prepared for these arguments, these outputs, this
 * element and this rule (out is not modified)
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::read_prepared_coefficient(
    const Coefficient& coef, const int var_blk, const unsigned int outputs,
    const std::vector<mfem::ParGridFunction>& aux_gf, const ValuesAt values_at, const int elem,
    const mfem::IntegrationRule& ir, CoefficientAtPoints& out) const {
  const auto ue = static_cast<std::size_t>(elem);
  if (!this->mesh_prepared_ || elem < 0 || ue >= this->mesh_rules_.size() ||
      this->mesh_rules_[ue] != &ir) {
    return false;
  }
  const auto it = std::ranges::find_if(this->mesh_coefficients_, [&](const MeshCoefficient& m) {
    return m.coef == &coef && m.var_blk == var_blk && m.aux_gf == &aux_gf &&
           m.values_at == values_at && (m.prepared_outputs & outputs) == outputs &&
           m.prepared_outputs != 0u;
  });
  if (it == this->mesh_coefficients_.end()) {
    return false;
  }

  // Points of the element: [first, first + np) in the tables of the local mesh
  const CoefficientAtPoints& prepared = it->at_points;
  const std::size_t first = this->mesh_offsets_[ue];
  const auto np = static_cast<std::size_t>(ir.GetNPoints());
  const auto nb = static_cast<std::size_t>(prepared.nb);
  const auto copy = [](const std::vector<double>& from, const std::size_t begin,
                       const std::size_t size, std::vector<double>& to) {
    const auto b = from.begin() + static_cast<std::ptrdiff_t>(begin);
    to.assign(b, b + static_cast<std::ptrdiff_t>(size));
  };
  out.nb = prepared.nb;
  out.is_explicit = prepared.is_explicit;
  out.by_batch = true;
  if ((outputs & CoefValue) != 0u) {
    copy(prepared.F, first, np, out.F);
  } else {
    out.F.clear();
  }
  if ((outputs & CoefFirstDerivative) != 0u) {
    copy(prepared.dF, first * nb, np * nb, out.dF);
  } else {
    out.dF.clear();
  }
  if ((outputs & CoefSecondDerivative) != 0u) {
    copy(prepared.d2F, first * nb * nb, np * nb * nb, out.d2F);
  } else {
    out.d2F.clear();
  }
  return true;
}
