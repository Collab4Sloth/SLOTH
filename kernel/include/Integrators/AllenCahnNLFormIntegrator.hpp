/**
 * @file AllenCahnNLFormIntegrator.hpp
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
#include <optional>
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
 * @brief Right-hand side of the Allen-Cahn equations (see the TimeDerivative integrator for the
 * left-hand side), with the free energy given in one of two ways:
 *
 * - classic: Capillary coefficient lambda_b and FreeEnergy W_b(u) of each block,
 *     R_b(psi) = int M_b (lambda_b grad u_b . grad psi + dW_b/du_b psi);
 * - variational: FreeEnergy f(u, grad u) of the first block, variational coefficient
 *   ("variational": true in the JSON file) of the n blocks and of their gradients, without
 *   Capillary coefficient (the gradient terms are part of f),
 *     R_b(psi) = int M_b (df/du_b psi + df/d(grad u_b) . grad psi).
 *   An isotropic f (sum_b kappa_b/2 |grad u_b|^2 + g(u), not explicit) is assembled with the
 *   operations of the classic form (kappa_b in place of lambda_b), the other ones with all the
 *   derivatives of f with respect to the blocks and to their gradients.
 *
 * Coefficients of each block: Mobility M_b, FreeEnergy and, in the classic form, Capillary.
 * Scheme of a variational free energy: implicit, explicit (values and gradients of u_old, no
 * contribution to the Jacobian) or semi-implicit.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class AllenCahnNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  /// Base class of the nonlinear integrators of SLOTH
  using SlothNLBase = SlothNLFormIntegrator<VARS>;

  std::list<GlossaryType> expected_list_{GlossaryType::Mobility, GlossaryType::Capillary,
                                         GlossaryType::FreeEnergy};
  /// Expected coefficients with a variational free energy (no Capillary)
  std::list<GlossaryType> expected_variational_list_{GlossaryType::Mobility,
                                                     GlossaryType::FreeEnergy};
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  std::vector<double> mob_at_ip_;   ///< mobility of each block at the current point
  std::vector<double> dmob_at_ip_;  ///< dM/du_j at the current point

  // Coefficients of each block at the integration points of the current element (blocks with
  // the same finite element, see SlothNLFormIntegrator::eval_coefficient_at_points)
  using CoefficientAtPoints = typename SlothNLBase::CoefficientAtPoints;
  std::vector<CoefficientAtPoints> mobility_at_points_;
  std::vector<CoefficientAtPoints> lambda_at_points_;
  std::vector<CoefficientAtPoints> energy_at_points_;
  mfem::Vector grad_u_mob_, g_mob_;  ///< work vectors of add_mobility_derivative

  // Theta-scheme (theta != 1): dofs of theta u + (1 - theta) u_old on the current element, used
  // in the gradient term
  std::vector<mfem::Vector> theta_dofs_;
  mfem::Vector uold_dofs_;  ///< work vector: element dofs of u_old
  void prepare_theta_dofs(mfem::ElementTransformation& Tr,
                          const mfem::Array<const mfem::Vector*>& elfun);
  const mfem::Vector& gradient_term_dofs(const mfem::Array<const mfem::Vector*>& elfun,
                                         const int blk) const;

  // Variational free energy
  bool variational_ = false;  ///< FreeEnergy of the first block variational
  /// Variational free energy isotropic and not explicit: classic assembly with kappa_b
  bool isotropic_ = false;
  /// Derivatives of the variational free energy at the points of the current element
  typename SlothNLBase::VariationalAtPoints energy_;
  std::vector<mfem::Vector> residual_integrand_;  ///< integrand of R_b at one point (each block)

  void coefficients_at_points(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                              const mfem::IntegrationRule& ir,
                              const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian);
  void mobility_at_points(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                          const mfem::IntegrationRule& ir,
                          const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian);
  double integration_weight(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip);
  void variational_element_vector(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                                  const mfem::IntegrationRule& ir,
                                  const mfem::Array<mfem::Vector*>& elvect);
  void variational_element_grad(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                                const mfem::IntegrationRule& ir,
                                const mfem::Array2D<mfem::DenseMatrix*>& elmats);

 protected:
  Coefficients lambda;
  Coefficients mobility;
  Coefficients double_well_energy;
  std::optional<Coefficient> free_energy;  ///< variational free energy

  std::vector<mfem::ParGridFunction> vaux_gf_;

  void get_coefficients() override;

  virtual void check_variables_consistency();

  void add_mobility_derivative(const int blk, const double lamb, const double w_prime,
                               const mfem::DenseMatrix& gradPsi, const mfem::Vector& Psi,
                               const mfem::Vector& elfun_blk, const double weight_coef,
                               const mfem::Array2D<mfem::DenseMatrix*>& elmats);

 public:
  void init() override;
  /// Theta-scheme on the gradient term (isotropic gradient energy)
  bool supports_theta() const override { return true; }
  AllenCahnNLFormIntegrator(Geometry geometry, const double time_step,
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

#include "Integrators/AllenCahnNLFormIntegrator.tpp"
