/**
 * @file VelocityField.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Imposed advection velocity shared by the advection integrators
 * @version 0.1
 * @date 2026-10-19
 *
 * @copyright CEA (C) 2026
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
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Imposed advection velocity of a block, used by the volume (Advection, SplitAdvection)
 * and boundary (AdvectionBoundary, AdvectionBoundaryCH) integrators.
 *
 * Two modes, selected from the coefficients of the block:
 *  - Components: v = (v_x, v_y[, v_z]) given by Glossary::Velocity_X/Y/Z (type Velocity,
 *    ID = component). All the components of the space dimension are required.
 *  - Gradient:   v = alpha grad(a), with alpha given by Glossary::VelocityPrefactor and a an
 *    auxiliary variable whose symbol is given by the parameter "velocity_gradient"
 *    (e.g. Parameter("velocity_gradient", "T")).
 *
 * In both modes, the coefficients can depend on the coordinates, the time and the auxiliary
 * variables, not on the unknowns (no dv/du in the Jacobian). The auxiliary variable a is frozen
 * in the Jacobian (staggered coupling).
 *
 * The gradient is evaluated with mfem::GridFunction::GetGradient, which also handles boundary
 * element transformations (gradient taken in the adjacent element).
 */
class VelocityField {
 public:
  enum class Mode { None, Components, Gradient };

 private:
  Mode mode_{Mode::None};
  int sdim_{0};
  std::vector<Coefficient> components_;   ///< Components mode
  std::optional<Coefficient> prefactor_;  ///< Gradient mode: alpha
  int grad_aux_{-1};                      ///< Gradient mode: index of a among the aux fields
  bool needs_aux_values_{false};
  mutable mfem::Vector grad_;

 public:
  VelocityField() = default;

  /**
   * @brief Configure the velocity of a block from its coefficients.
   *
   * @param name         Integrator name (error messages).
   * @param blk          Block index (error messages).
   * @param sdim         Space dimension.
   * @param get_coef     Callable (GlossaryType, id) -> std::optional<Coefficient> for the block.
   * @param grad_symbol  Symbol of the auxiliary variable for the Gradient mode ("" if none).
   * @param aux_infos    Additional information of the auxiliary variables (symbol = last item).
   * @return true if the block has a velocity (Components or Gradient mode).
   */
  template <class GetCoef>
  bool configure(const std::string& name, int blk, int sdim, GetCoef&& get_coef,
                 const std::string& grad_symbol,
                 const std::vector<std::vector<std::string>>& aux_infos) {
    this->mode_ = Mode::None;
    this->sdim_ = sdim;
    this->components_.clear();
    this->prefactor_.reset();
    this->grad_aux_ = -1;
    this->needs_aux_values_ = false;

    const auto depends_on_aux = [](const Coefficient& c) { return !c.is_scalar() && !c.is_file(); };
    const std::string where = name + "NLFormIntegrator: block " + std::to_string(blk);

    for (int d = 0; d < sdim; ++d) {
      auto coef = get_coef(GlossaryType::Velocity, d);
      if (coef.has_value()) {
        this->needs_aux_values_ = this->needs_aux_values_ || depends_on_aux(*coef);
        this->components_.push_back(*coef);
      }
    }
    auto prefactor = get_coef(GlossaryType::VelocityPrefactor, 0);

    if (!prefactor.has_value()) {
      if (this->components_.empty()) {
        return false;
      }
      MFEM_VERIFY(static_cast<int>(this->components_.size()) == sdim,
                  where + " holds only some of the velocity components (" + std::to_string(sdim) +
                      " expected). Please check your data.");
      this->mode_ = Mode::Components;
      return true;
    }

    MFEM_VERIFY(this->components_.empty(),
                where +
                    ": Velocity components and VelocityPrefactor are both given. Please "
                    "check your data.");
    MFEM_VERIFY(!grad_symbol.empty(),
                where +
                    ": a VelocityPrefactor coefficient is given but the parameter "
                    "\"velocity_gradient\" (symbol of the auxiliary variable) is missing.");
    const std::string symbol = toUpperCase(grad_symbol);
    for (std::size_t i = 0; i < aux_infos.size(); ++i) {
      if (!aux_infos[i].empty() && toUpperCase(aux_infos[i].back()) == symbol) {
        this->grad_aux_ = static_cast<int>(i);
        break;
      }
    }
    MFEM_VERIFY(this->grad_aux_ >= 0, where + ": auxiliary variable \"" + grad_symbol +
                                          "\" (parameter velocity_gradient) not found.");
    this->prefactor_ = *prefactor;
    this->needs_aux_values_ = depends_on_aux(*prefactor);
    this->mode_ = Mode::Gradient;
    return true;
  }

  Mode mode() const { return this->mode_; }
  bool active() const { return this->mode_ != Mode::None; }
  /// The coefficients need the values of the auxiliary fields at the integration point
  bool needs_aux_values() const { return this->needs_aux_values_; }

  /**
   * @brief Velocity at the current integration point (Tr must be set at ip).
   *
   * @param Tr       Element (or boundary element) transformation, set at ip.
   * @param aux_gf   Auxiliary grid functions (same order as aux_infos in configure()).
   * @param coef     Callable (const Coefficient&) -> double evaluating a coefficient at ip.
   * @param v        Velocity (resized to sdim).
   */
  template <class CoefEval>
  void eval(mfem::ElementTransformation& Tr, const std::vector<mfem::ParGridFunction>& aux_gf,
            CoefEval&& coef, mfem::Vector& v) const {
    v.SetSize(this->sdim_);
    switch (this->mode_) {
      case Mode::Components:
        for (int d = 0; d < this->sdim_; ++d) {
          v(d) = coef(this->components_[d]);
        }
        break;
      case Mode::Gradient: {
        aux_gf[this->grad_aux_].GetGradient(Tr, this->grad_);
        const double alpha = coef(*this->prefactor_);
        for (int d = 0; d < this->sdim_; ++d) {
          v(d) = alpha * this->grad_(d);
        }
        break;
      }
      case Mode::None:
      default:
        v = 0.0;
        break;
    }
  }
};
