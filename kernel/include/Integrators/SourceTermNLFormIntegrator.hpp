/**
 * @file SourceTermNLFormIntegrator.hpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @brief VF of a user-defined source term
 * @version 0.4
 * @date 2026-10-16
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
#include <list>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class dedicated to user-defined source terms.
 *
 * For each block i, with the source term f_i (coefficient of type SourceTerm, ID 0):
 *   R_i  = - (f_i, psi_i)
 *   J_ij = - (df_i/du_j phi_j, psi_i)
 *
 * f_i can be a constant, a function of the unknowns, of the auxiliary variables, of the time
 * ("constants": "(t:T)") and of the spatial coordinates ("coordinates"), or read from a file.
 * In axisymmetric geometry, the r weight is applied by the integrator.
 *
 * A SourceTerm coefficient is required for every block (use a zero constant for blocks without
 * source term). The Jacobian is only assembled for source terms that can depend on the current
 * solution (neither constant, explicit nor file-based).
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class SourceTermNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  std::list<GlossaryType> expected_list_{GlossaryType::SourceTerm};

  // Work arrays reused across elements and integration points (one entry per block)
  std::vector<mfem::Vector> shape_;       ///< shape functions of each block at ip
  std::vector<mfem::Vector> u_old_dofs_;  ///< local dofs of u_old for each block
  std::vector<double> u_values_;          ///< [u_0..u_{n-1}, u_old_0..u_old_{n-1}] at ip
  std::vector<double> aux_values_;        ///< auxiliary fields at ip

  void prepare_element(const mfem::Array<const mfem::FiniteElement*>& el,
                       mfem::ElementTransformation& Tr);

  double evaluate_at_ip(const mfem::Array<const mfem::FiniteElement*>& el,
                        const mfem::Array<const mfem::Vector*>& elfun,
                        mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip);

  static const mfem::IntegrationRule& get_integration_rule(
      const mfem::Array<const mfem::FiniteElement*>& el, const mfem::ElementTransformation& Tr);

 protected:
  Coefficients source_term_;          ///< one source term per block (positional)
  std::vector<bool> has_derivative_;  ///< per block: df/du may be non-zero
  bool any_derivative_{false};        ///< at least one block has a non-zero Jacobian
  bool needs_fields_{false};          ///< at least one source term is not constant (u, aux needed)

  std::vector<mfem::ParGridFunction> vaux_gf_;

  void get_coefficients() override;

 public:
  void init() override;
  SourceTermNLFormIntegrator(Geometry geometry, const double time_step,
                             const std::vector<mfem::ParGridFunction>& u_old,
                             const std::vector<mfem::ParGridFunction>& aux_old,
                             const Parameters& params, std::vector<VARS*> auxvars,
                             const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
};

#include "Integrators/SourceTermNLFormIntegrator.tpp"
