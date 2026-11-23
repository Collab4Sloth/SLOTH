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
#include "Coefficients/FieldsAtPoints.hpp"
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
  out.by_batch = this->eval_by_batch(coef, var_blk, outputs, Tr, ir, values, out);
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
 * in values_at_end_of_step. For the energies of the system (ValuesAt::Current), {u, u_old} is
 * replaced by {u, u}.
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
  } else if (values_at == ValuesAt::Current) {
    // [u, u_old] -> [u, u] (same work array as the end of the step)
    this->cap_u_end_of_step_.resize(np * width);
    const std::size_t nb = values.nb;
    for (std::size_t q = 0; q < np; ++q) {
      const std::span<const double> u_uold = values.u_uold_at(q);
      const std::span<double> u_u =
          std::span<double>(this->cap_u_end_of_step_).subspan(q * width, width);
      for (std::size_t b = 0; b < nb; ++b) {
        u_u[b] = u_uold[b];
        u_u[nb + b] = u_uold[b];
      }
    }
    values.u_uold = this->cap_u_end_of_step_;
  }
  return values;
}

/**
 * @brief Evaluation by batch of a coefficient with an evaluation structure (see
 * FunctionCoefficient::eval_batch), at all the integration points of an element.
 *
 * The values given to the coefficient are selected according to its scheme and, for a
 * coefficient of one variable, its block; the derivatives returned with respect to its variables
 * are then stored with respect to the blocks of the integrator.
 *
 * @return false if the coefficient has no evaluation structure (nothing is computed)
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::eval_by_batch(const Coefficient& coef, const int var_blk,
                                                const unsigned int outputs,
                                                mfem::ElementTransformation& Tr,
                                                const mfem::IntegrationRule& ir,
                                                const ElementValues& values,
                                                CoefficientAtPoints& out) {
  if (coef.batch_nb_variables() <= 0) {
    return false;
  }
  const bool want_value = (outputs & CoefValue) != 0u;
  const bool want_first_derivative = (outputs & CoefFirstDerivative) != 0u;
  const bool want_second_derivative = (outputs & CoefSecondDerivative) != 0u;
  const int nip = ir.GetNPoints();
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
  if (coef.uses_coordinates()) {
    Tr.Transform(ir, this->cap_points_);  // dim x nip, column q = point q
    const std::span<const double> points(this->cap_points_.GetData(),
                                         static_cast<std::size_t>(this->cap_points_.Height()) * np);
    out.coordinates.SetSize(static_cast<int>(points.size()));
    std::ranges::copy(points, out.coordinates.HostWrite());
  } else {
    out.coordinates.SetSize(0);
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

/**
 * @brief Values and physical gradients of the auxiliary fields at the integration points of an
 * element, with the layout of fields_at_points (Coefficients/FieldsAtPoints.hpp)
 *
 * From their element dofs if all of them share the finite element fe, point by point otherwise.
 *
 * @param fe Finite element of the element
 * @param Tr Transformation of the element
 * @param ir Integration rule
 * @param aux_gf Auxiliary fields
 * @param values Values at the points (resized to nip * naux)
 * @param gradients Gradients at the points (resized to nip * naux * dim), or nullptr if they are
 *        not needed
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::auxiliary_fields_at_points(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const std::vector<mfem::ParGridFunction>& aux_gf, mfem::Vector& values,
    mfem::Vector* gradients) {
  const int elem = Tr.ElementNo;
  const std::size_t naux = aux_gf.size();
  this->afp_dofs_.resize(naux);
  bool aux_dofs = true;
  for (std::size_t k = 0; k < naux && aux_dofs; ++k) {
    const mfem::FiniteElementSpace* fes = aux_gf[k].FESpace();
    aux_dofs = elem < fes->GetNE() && this->same_finite_element(*fes->GetFE(elem), fe);
    if (aux_dofs) aux_gf[k].GetElementDofValues(elem, this->afp_dofs_[k]);
  }
  if (aux_dofs) {
    fields_at_points(fe, Tr, ir, this->afp_dofs_, values, gradients);
    return;
  }
  const int nip = ir.GetNPoints();
  const int sdim = Tr.GetSpaceDim();
  const auto unip = static_cast<std::size_t>(nip);
  const auto usdim = static_cast<std::size_t>(sdim);
  values.SetSize(static_cast<int>(unip * naux));
  if (gradients != nullptr) gradients->SetSize(static_cast<int>(unip * naux * usdim));
  double* v = values.HostWrite();
  double* g = gradients != nullptr ? gradients->HostWrite() : nullptr;
  this->afp_grad_.SetSize(sdim);
  for (int q = 0; q < nip; ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    this->set_int_point(Tr, ip);
    const auto uq = static_cast<std::size_t>(q);
    for (std::size_t k = 0; k < naux; ++k) {
      v[uq * naux + k] = aux_gf[k].GetValue(Tr, ip);
      if (g == nullptr) continue;
      aux_gf[k].GetGradient(Tr, this->afp_grad_);
      for (std::size_t c = 0; c < usdim; ++c) {
        g[(uq * naux + k) * usdim + c] = this->afp_grad_(static_cast<int>(c));
      }
    }
  }
}

/**
 * @brief Coefficients of the energies of the system: FreeEnergy and, if given, GradientEnergy of
 * block 0. Must be called before system_energies.
 *
 * @return false if there is no FreeEnergy coefficient
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::prepare_system_energies() {
  this->system_free_energy_.reset();
  this->system_gradient_energy_.reset();
  if (this->coefficients_.empty()) return false;
  this->system_free_energy_ = this->get_coefficient(0, GlossaryType::FreeEnergy, 0);
  this->system_gradient_energy_ = this->get_coefficient(0, GlossaryType::GradientEnergy, 0);
  if (!this->system_free_energy_.has_value()) return false;
  MFEM_VERIFY(
      !(this->system_free_energy_->is_variational() && this->system_gradient_energy_.has_value()),
      "Error: a variational FreeEnergy coefficient already contains the gradient terms: a "
      "GradientEnergy coefficient must not be given.");
  return true;
}

/**
 * @brief True if the energies of the system have an interface part: variational FreeEnergy, or
 * GradientEnergy given
 */
template <class VARS>
bool SlothNLFormIntegrator<VARS>::has_interface_energy() const {
  return this->system_gradient_energy_.has_value() ||
         (this->system_free_energy_.has_value() && this->system_free_energy_->is_variational());
}

/**
 * @brief Energies of the system on an element, from the coefficients of block 0 (see
 * prepare_system_energies), with the current values of the variables (elfun) and of the
 * auxiliary fields:
 *
 * - FreeEnergy f(u) and GradientEnergy g(grad u):
 *     density = int (f + g),        interface = 2 int g
 * - variational FreeEnergy f(u, grad u):
 *     density = int f(u, grad u),   interface = 2 int (f(u, grad u) - f(u, 0))
 *
 * Explicit and semi-implicit coefficients are evaluated with the current values (explicit
 * variables equal to the variables). Integration rule of order 2 p + 3.
 */
template <class VARS>
typename SlothNLFormIntegrator<VARS>::SystemEnergies SlothNLFormIntegrator<VARS>::system_energies(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun) {
  MFEM_VERIFY(this->system_free_energy_.has_value(),
              "Error: energies of the system without FreeEnergy coefficient (see "
              "prepare_system_energies).");
  MFEM_VERIFY(same_finite_elements(el),
              "Error: energies of the system: the variables must share the same finite element.");
  const mfem::FiniteElement& fe = *el[0];
  const mfem::IntegrationRule& ir = mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + 3);
  // Tables of prepare_element_values: values of the assembly, not valid here
  this->ipv_ir_ = nullptr;
  this->ipv_element_ = -1;

  const auto nb = static_cast<std::size_t>(elfun.Size());
  this->se_dofs_.resize(nb);
  for (std::size_t b = 0; b < nb; ++b) this->se_dofs_[b] = *elfun[static_cast<int>(b)];

  SystemEnergies energies;
  if (this->system_free_energy_->is_variational()) {
    this->variational_system_energies(fe, Tr, ir, energies.density, energies.interface);
  } else {
    this->classic_system_energies(fe, Tr, ir, elfun, energies.density, energies.interface);
  }
  return energies;
}

/**
 * @brief Energies of the system with a FreeEnergy f(u) and, if given, a GradientEnergy g(grad u)
 * (see system_energies)
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::classic_system_energies(
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, double& density, double& interface) {
  this->eval_coefficient_at_points(*this->system_free_energy_, 0, CoefValue, fe, Tr, ir, elfun,
                                   this->aux_gf_, this->se_free_energy_, ValuesAt::Current);

  const bool with_gradient = this->system_gradient_energy_.has_value();
  const int dim = Tr.GetSpaceDim();
  const auto ud = static_cast<std::size_t>(dim);
  const std::size_t nb = this->se_dofs_.size();
  const std::size_t naux = this->aux_gf_.size();
  if (with_gradient) {
    fields_at_points(fe, Tr, ir, this->se_dofs_, this->se_values_, &this->se_gradients_);
    this->auxiliary_fields_at_points(fe, Tr, ir, this->aux_gf_, this->se_aux_values_,
                                     &this->se_aux_gradients_);
  }
  const double* grad = with_gradient ? this->se_gradients_.HostRead() : nullptr;
  const double* aux_grad = with_gradient && naux > 0 ? this->se_aux_gradients_.HostRead() : nullptr;

  for (int q = 0; q < ir.GetNPoints(); ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    Tr.SetIntPoint(&ip);
    const double w = ip.weight * Tr.Weight();
    double g = 0.;
    if (with_gradient) {
      const auto uq = static_cast<std::size_t>(q);
      const std::span<const double> grad_q(grad + uq * nb * ud, nb * ud);
      const std::span<const double> aux_grad_q(
          aux_grad != nullptr ? aux_grad + uq * naux * ud : nullptr, naux * ud);
      g = this->system_gradient_energy_->compute(grad_q, aux_grad_q, dim);
    }
    density += (this->se_free_energy_.value(q) + g) * w;
    interface += 2. * g * w;
  }
}

/**
 * @brief Energies of the system with a variational FreeEnergy f(u, grad u) (see
 * system_energies). f(u, 0) gathers the terms which do not depend on the gradients: the
 * interface energy is twice the energy of the gradient terms, as with a GradientEnergy
 * coefficient.
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::variational_system_energies(const mfem::FiniteElement& fe,
                                                              mfem::ElementTransformation& Tr,
                                                              const mfem::IntegrationRule& ir,
                                                              double& density, double& interface) {
  const Coefficient& coef = *this->system_free_energy_;
  const int sdim = Tr.GetSpaceDim();
  MFEM_VERIFY(coef.variational_dimension() == sdim,
              "Error: dimension of the variational FreeEnergy coefficient (" +
                  std::to_string(coef.variational_dimension()) +
                  ") different from the dimension of the mesh (" + std::to_string(sdim) + ").");
  const bool semi = coef.is_semi_implicit();
  const bool with_aux_gradients = coef.uses_auxiliary_gradients();
  const bool with_exp_gradients = semi && coef.uses_explicit_gradients();
  const int nip = ir.GetNPoints();

  fields_at_points(fe, Tr, ir, this->se_dofs_, this->se_values_, &this->se_gradients_);
  this->auxiliary_fields_at_points(fe, Tr, ir, this->aux_gf_, this->se_aux_values_,
                                   with_aux_gradients ? &this->se_aux_gradients_ : nullptr);
  if (!with_aux_gradients) this->se_aux_gradients_.SetSize(0);
  if (coef.uses_coordinates()) {
    Tr.Transform(ir, this->se_points_);  // sdim x nip, column q = point q
    this->se_coordinates_.SetSize(this->se_points_.Height() * nip);
    std::copy(this->se_points_.GetData(), this->se_points_.GetData() + this->se_coordinates_.Size(),
              this->se_coordinates_.HostWrite());
  } else {
    this->se_coordinates_.SetSize(0);
  }

  // f(u, grad u), then f(u, 0): gradients set to zero
  const mfem::Vector none;
  this->se_F_.SetSize(nip);
  this->se_F0_.SetSize(nip);
  VariationalBatch outputs;
  outputs.F = &this->se_F_;
  MFEM_VERIFY(coef.eval_variational_batch(nip, this->se_values_, semi ? this->se_values_ : none,
                                          this->se_aux_values_, this->se_gradients_,
                                          with_exp_gradients ? this->se_gradients_ : none,
                                          this->se_aux_gradients_, this->se_coordinates_, outputs),
              "Error: evaluation of the variational FreeEnergy coefficient.");
  this->se_zero_gradients_.SetSize(this->se_gradients_.Size());
  this->se_zero_gradients_ = 0.0;
  this->se_zero_aux_gradients_.SetSize(this->se_aux_gradients_.Size());
  this->se_zero_aux_gradients_ = 0.0;
  outputs.F = &this->se_F0_;
  MFEM_VERIFY(coef.eval_variational_batch(
                  nip, this->se_values_, semi ? this->se_values_ : none, this->se_aux_values_,
                  this->se_zero_gradients_, with_exp_gradients ? this->se_zero_gradients_ : none,
                  this->se_zero_aux_gradients_, this->se_coordinates_, outputs),
              "Error: evaluation of the variational FreeEnergy coefficient.");

  const double* f = this->se_F_.HostRead();
  const double* f0 = this->se_F0_.HostRead();
  for (int q = 0; q < nip; ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    Tr.SetIntPoint(&ip);
    const double w = ip.weight * Tr.Weight();
    density += f[q] * w;
    interface += 2. * (f[q] - f0[q]) * w;
  }
}

/**
 * @brief Derivatives of a variational free energy f at all the integration points of an element,
 * by batch (eval_variational_batch). The variables of f are the first nb_variables blocks:
 * first derivatives with respect to the variables and to the components of their gradients and,
 * for the Jacobian of an energy which is not explicit, second derivatives. Isotropic energy
 * (isotropic true): kappa_b in place of the derivatives with respect to the gradients
 * (kappa_b grad u_b and kappa_b I).
 *
 * Values given to the energy: values and gradients of u (implicit energy) or of u_old (explicit
 * energy); for a semi-implicit energy, variables u and explicit variables u_old. Auxiliary
 * variables: values and gradients of the auxiliary fields.
 *
 * @param f Variational free energy
 * @param nb_variables Number of variables of f (first blocks)
 * @param isotropic Isotropic energy: kappa_b computed in place of the derivatives with respect
 *        to the gradients
 * @param fe Finite element of the blocks
 * @param Tr Transformation of the element
 * @param ir Integration rule
 * @param elfun Values of the dofs of each block
 * @param aux_gf Auxiliary fields
 * @param jacobian Second derivatives required by the Jacobian
 * @param v Derivatives at the points
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::variational_at_points(
    const Coefficient& f, const int nb_variables, const bool isotropic,
    const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
    const mfem::Array<const mfem::Vector*>& elfun, const std::vector<mfem::ParGridFunction>& aux_gf,
    const bool jacobian, VariationalAtPoints& v) {
  const int nb = nb_variables;
  const auto unb = static_cast<std::size_t>(nb);
  const int nip = ir.GetNPoints();
  const int sdim = Tr.GetSpaceDim();
  const int elem = Tr.ElementNo;
  MFEM_VERIFY(f.variational_dimension() == sdim,
              "Error in " + this->integrator_name_ + ": dimension of the free energy (" +
                  std::to_string(f.variational_dimension()) +
                  ") different from the dimension of the mesh (" + std::to_string(sdim) + ").");
  MFEM_VERIFY(f.variational_nb_variables() == nb,
              "Error in " + this->integrator_name_ + ": free energy of " +
                  std::to_string(f.variational_nb_variables()) + " variables for " +
                  std::to_string(nb) + " blocks: the variables of the free energy are the blocks.");
  MFEM_VERIFY(nb <= elfun.Size(), "Error in " + this->integrator_name_ +
                                      ": more variables of the free energy than blocks.");
  const bool explicit_energy = f.is_explicit();
  const bool semi_implicit = f.is_semi_implicit();
  const bool explicit_gradients = semi_implicit && f.uses_explicit_gradients();
  const bool auxiliary_gradients = f.uses_auxiliary_gradients();
  v.nb = nb;
  v.dim = sdim;
  v.isotropic = isotropic;

  // Element dofs of u and, if needed, of u_old
  v.dofs_u.resize(unb);
  for (std::size_t b = 0; b < unb; ++b) {
    v.dofs_u[b] = *elfun[static_cast<int>(b)];
  }
  if (explicit_energy || semi_implicit) {
    v.dofs_u_old.resize(unb);
    for (std::size_t b = 0; b < unb; ++b) {
      const mfem::FiniteElementSpace* fes = this->u_old_[b].FESpace();
      MFEM_VERIFY(elem < fes->GetNE() && this->same_finite_element(*fes->GetFE(elem), fe),
                  "Error in " + this->integrator_name_ +
                      ": the previous values of the blocks must be on the finite element space of "
                      "the blocks.");
      this->u_old_[b].GetElementDofValues(elem, v.dofs_u_old[b]);
    }
  }

  // Variables (u or u_old) and explicit variables (u_old) with their gradients
  fields_at_points(fe, Tr, ir, explicit_energy ? v.dofs_u_old : v.dofs_u, v.values, &v.gradients);
  if (semi_implicit) {
    fields_at_points(fe, Tr, ir, v.dofs_u_old, v.explicit_values,
                     explicit_gradients ? &v.explicit_gradients : nullptr);
  } else {
    v.explicit_values.SetSize(0);
  }
  if (!explicit_gradients) v.explicit_gradients.SetSize(0);

  // Auxiliary variables with their gradients
  this->auxiliary_fields_at_points(fe, Tr, ir, aux_gf, v.auxiliary_values,
                                   auxiliary_gradients ? &v.auxiliary_gradients : nullptr);
  if (!auxiliary_gradients) v.auxiliary_gradients.SetSize(0);

  // Physical coordinates of the points
  if (f.uses_coordinates()) {
    Tr.Transform(ir, v.points);  // sdim x nip, column q = point q
    v.coordinates.SetSize(v.points.Height() * nip);
    std::copy(v.points.GetData(), v.points.GetData() + v.coordinates.Size(),
              v.coordinates.HostWrite());
  } else {
    v.coordinates.SetSize(0);
  }

  // Derivatives of the free energy (second ones: no contribution of an explicit energy)
  const int nd = nb * sdim;
  v.with_second_derivatives = jacobian && !explicit_energy;
  v.df_dvar.SetSize(nip * nb);
  VariationalBatch outputs;
  outputs.df_dvar = &v.df_dvar;
  if (isotropic) {
    v.kappa.SetSize(nip * nb);
    outputs.d2f_dgrad2_isotropic = &v.kappa;
  } else {
    v.df_dgrad.SetSize(nip * nd);
    outputs.df_dgrad = &v.df_dgrad;
  }
  if (v.with_second_derivatives) {
    v.d2f_dvar2.SetSize(nip * nb * nb);
    outputs.d2f_dvar2 = &v.d2f_dvar2;
    if (!isotropic) {
      v.d2f_dvar_dgrad.SetSize(nip * nb * nd);
      v.d2f_dgrad2.SetSize(nip * nd * nd);
      outputs.d2f_dvar_dgrad = &v.d2f_dvar_dgrad;
      outputs.d2f_dgrad2 = &v.d2f_dgrad2;
    }
  }
  MFEM_VERIFY(
      f.eval_variational_batch(nip, v.values, v.explicit_values, v.auxiliary_values, v.gradients,
                               v.explicit_gradients, v.auxiliary_gradients, v.coordinates, outputs),
      "Error in " + this->integrator_name_ + ": evaluation of the free energy.");
}

/**
 * @brief Adds to elvect the derivative of int f in the direction of the test functions psi of
 * the variable b of a variational free energy, at the point q, multiplied by a:
 *   a (df/du_b psi + df/d(grad u_b) . grad psi)
 * (derivatives of variational_at_points, energy which is not isotropic).
 *
 * @param v Derivatives of the free energy at the points
 * @param q Integration point
 * @param b Variable of f (block)
 * @param a Factor (mobility, integration weight, sign...)
 * @param Psi Shape functions at the point
 * @param gradPsi Physical gradients of the shape functions at the point
 * @param elvect Element vector
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::add_variational_vector(VariationalAtPoints& v, const int q,
                                                         const int b, const double a,
                                                         const mfem::Vector& Psi,
                                                         const mfem::DenseMatrix& gradPsi,
                                                         mfem::Vector& elvect) {
  const int sdim = v.dim;
  const int ng = v.nb * sdim;  // components of the gradients of the variables
  const double* df_dvar = v.df_dvar.HostRead();
  const double* df_dgrad = v.df_dgrad.HostRead();
  v.flux.SetSize(sdim);
  // (df/d(grad u_b), grad psi)
  for (int k = 0; k < sdim; ++k) v.flux(k) = df_dgrad[q * ng + b * sdim + k];
  gradPsi.AddMult_a(a, v.flux, elvect);
  // (df/du_b, psi)
  elvect.Add(a * df_dvar[q * v.nb + b], Psi);
}

/**
 * @brief Adds to mat the second derivative of int f in the directions of the test functions psi
 * of the variable b and psi_j of the variable j, at the point q, multiplied by a:
 *   a (d2f/du_b du_j psi psi_j + (d2f/du_b d(grad u_j) . grad psi_j) psi
 *      + (d2f/d(grad u_b) du_j . grad psi) psi_j
 *      + grad psi . d2f/d(grad u_b) d(grad u_j) grad psi_j)
 * (second derivatives of variational_at_points, energy which is not isotropic nor explicit).
 *
 * @param v Derivatives of the free energy at the points
 * @param q Integration point
 * @param b Variable of f of the test functions (row block)
 * @param j Variable of f of the trial functions (column block)
 * @param a Factor (mobility, integration weight, sign...)
 * @param Psi Shape functions at the point
 * @param gradPsi Physical gradients of the shape functions at the point
 * @param mat Element matrix of the blocks (b, j)
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::add_variational_grad(VariationalAtPoints& v, const int q,
                                                       const int b, const int j, const double a,
                                                       const mfem::Vector& Psi,
                                                       const mfem::DenseMatrix& gradPsi,
                                                       mfem::DenseMatrix& mat) {
  const int nb = v.nb;
  const int sdim = v.dim;
  const int nd = gradPsi.Height();
  const int ng = nb * sdim;  // components of the gradients of the variables
  const double* d2f_dvar2 = v.d2f_dvar2.HostRead();
  const double* d2f_dvar_dgrad = v.d2f_dvar_dgrad.HostRead();
  const double* d2f_dgrad2 = v.d2f_dgrad2.HostRead();
  v.flux.SetSize(sdim);
  v.row.SetSize(nd);
  v.column.SetSize(nd);
  v.block_hessian.SetSize(sdim, sdim);
  v.gradPsi_hessian.SetSize(nd, sdim);
  // d2f/du_b du_j psi psi_j
  mfem::AddMult_a_VVt(a * d2f_dvar2[(q * nb + b) * nb + j], Psi, mat);
  // (d2f/du_b d(grad u_j) . grad psi_j) psi
  for (int k = 0; k < sdim; ++k) v.flux(k) = d2f_dvar_dgrad[(q * nb + b) * ng + j * sdim + k];
  gradPsi.Mult(v.flux, v.column);
  mfem::AddMult_a_VWt(a, Psi, v.column, mat);
  // (d2f/d(grad u_b) du_j . grad psi) psi_j
  for (int k = 0; k < sdim; ++k) v.flux(k) = d2f_dvar_dgrad[(q * nb + j) * ng + b * sdim + k];
  gradPsi.Mult(v.flux, v.row);
  mfem::AddMult_a_VWt(a, v.row, Psi, mat);
  // grad psi . d2f/d(grad u_b) d(grad u_j) grad psi_j
  for (int k = 0; k < sdim; ++k) {
    for (int l = 0; l < sdim; ++l) {
      v.block_hessian(k, l) = d2f_dgrad2[(q * ng + b * sdim + k) * ng + j * sdim + l];
    }
  }
  mfem::Mult(gradPsi, v.block_hessian, v.gradPsi_hessian);
  mfem::AddMult_a_ABt(a, v.gradPsi_hessian, gradPsi, mat);
}
