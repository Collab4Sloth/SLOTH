/**
 * @file CahnHilliardNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the CahnHilliard equations
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
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class dedicated to the VF of the RHS of the Cahn-Hilliard equation (splitted form) (see
 * SplitTimeDerivative integrator for the LHS). Classic form: blocks phi (block 0) and mu (block 1).
 *
 * Free energy given in one of two ways:
 * - classic: Capillary coefficient lambda and FreeEnergy W(phi) of block 0,
 *     R_phi(psi) = int (mu - W'(phi)) psi - lambda grad phi . grad psi;
 * - variational: n concentrations c_0..c_{n-1} (blocks 0..n-1) and n chemical potentials
 *   mu_0..mu_{n-1} (blocks n..2n-1); FreeEnergy f(c, grad c) of block 0, variational coefficient
 *   ("variational": true in the JSON file) of the n concentrations and of their gradients,
 *   without Capillary coefficient (the gradient terms are part of f),
 *     R_{c_i}(psi)  = int mu_i psi - df/dc_i psi - df/d(grad c_i) . grad psi,
 *     R_{mu_i}(psi) = int M_i grad mu_i . grad psi   (Mobility of block n + i).
 *   An isotropic f (sum_i kappa_i/2 |grad c_i|^2 + g(c), not explicit) is assembled with the
 *   operations of the classic form (kappa_i in place of lambda).
 * Classic form: R_mu(psi) = int M grad mu . grad psi.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class CahnHilliardNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  /// Base class of the nonlinear integrators of SLOTH
  using SlothNLBase = SlothNLFormIntegrator<VARS>;

  std::list<GlossaryType> expected_list_{GlossaryType::Mobility, GlossaryType::Capillary,
                                         GlossaryType::FreeEnergy};
  /// Expected coefficients with a variational free energy (no Capillary)
  std::list<GlossaryType> expected_variational_list_{GlossaryType::Mobility,
                                                     GlossaryType::FreeEnergy};
  bool variational_ = false;  ///< FreeEnergy of block 0 variational
  /// Variational free energy isotropic and not explicit: classic assembly with kappa
  bool isotropic_ = false;
  /// Derivatives of the variational free energy at the points of the current element
  typename SlothNLBase::VariationalAtPoints energy_;
  /// Mobility of each chemical potential mu_i at the points of the current element
  std::vector<typename SlothNLBase::CoefficientAtPoints> potential_mobility_;
  mfem::Vector grad_mu_, adv_;  ///< work vectors of variational_element_grad

  int nb_concentrations(const int nb_blocks) const;
  void potential_mobility_at_points(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                                    const mfem::IntegrationRule& ir,
                                    const mfem::Array<const mfem::Vector*>& elfun, const int n,
                                    const bool jacobian);
  void variational_element_vector(const mfem::Array<const mfem::FiniteElement*>& el,
                                  mfem::ElementTransformation& Tr,
                                  const mfem::Array<const mfem::Vector*>& elfun,
                                  const mfem::Array<mfem::Vector*>& elvect);
  void variational_element_grad(const mfem::Array<const mfem::FiniteElement*>& el,
                                mfem::ElementTransformation& Tr,
                                const mfem::Array<const mfem::Vector*>& elfun,
                                const mfem::Array2D<mfem::DenseMatrix*>& elmats);
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  void check_variables_consistency();

  // Coefficients at the integration points of the current element (see
  // SlothNLFormIntegrator::eval_coefficient_at_points)
  using CoefficientAtPoints = typename SlothNLBase::CoefficientAtPoints;
  CoefficientAtPoints lambda_at_points_;    ///< capillary coefficient (block 0)
  CoefficientAtPoints energy_at_points_;    ///< free energy (block 0)
  CoefficientAtPoints mobility_at_points_;  ///< mobility (block 1)

  // Theta-scheme (theta != 1): dofs of theta c + (1 - theta) c_old on the current element for the
  // concentration blocks, used in the gradient term of the equation of mu. The chemical
  // potentials are not interpolated (unknowns at the time n + theta).
  std::vector<mfem::Vector> theta_dofs_;
  mfem::Vector uold_dofs_;  ///< work vector: element dofs of c_old
  void prepare_theta_dofs(mfem::ElementTransformation& Tr,
                          const mfem::Array<const mfem::Vector*>& elfun, const int n);
  const mfem::Vector& gradient_term_dofs(const mfem::Array<const mfem::Vector*>& elfun,
                                         const int blk) const;

 protected:
  Coefficients lambda;
  Coefficients mobility;
  Coefficients double_well_energy;
  std::optional<Coefficient> free_energy;  ///< variational free energy of phi

  std::vector<mfem::ParGridFunction> vaux_gf_;

  void get_coefficients() override;

 public:
  void init() override;
  /// Theta-scheme on the gradient term of the equation of mu (isotropic gradient energy)
  bool supports_theta() const override { return true; }
  CahnHilliardNLFormIntegrator(Geometry geometry, const double time_step,
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

#include "Integrators/CahnHilliardNLFormIntegrator.tpp"
