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
