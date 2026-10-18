/**
 * @file AdvectionNLFormIntegrator.hpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the conservative advection term div(v u)
 * @version 0.5
 * @date 2026-10-17
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
#include <span>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Volume part of the conservative advection term div(v u_c), integrated by parts.
 *
 * The velocity v of an advected variable c is imposed through the coefficients of its block:
 * Glossary::Velocity_X, Velocity_Y (and Velocity_Z in 3D), i.e. type Velocity with ID = index of
 * the component. They can be constants or functions of the coordinates ("coordinates"), of the
 * time ("constants": "(t:T)") and of the auxiliary variables, but NOT of the unknowns (the
 * Jacobian contains no dv/du term). In axisymmetric geometry, components are (v_r, v_z).
 *
 * Each term (r, c) adds to the residual of block r (linear in u_c):
 *   R_r  = - ( u_c v, grad(psi_r) )
 *   J_rc = - ( v.grad(psi_r), phi_c )
 *
 * No div(v) is needed: valid for non-solenoidal velocities, int(u_c) is conserved exactly on a
 * closed domain (v.n = 0). On open boundaries (v.n != 0), an Inflow coefficient attached to the
 * boundary activates the inflow/outflow term (AdvectionBoundaryNLFormIntegrator); otherwise the
 * natural condition is a zero advective flux.
 *
 * Default (Advection): one term (b, b) per block having Velocity coefficients.
 * Derived classes redefine the terms through build_terms() (see SplitAdvection) and may
 * redefine the velocity through compute_velocity().
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class AdvectionNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  // Work arrays reused across elements and integration points
  std::vector<mfem::Vector> shape_;       ///< shape functions of each block at ip
  std::vector<mfem::Vector> u_old_dofs_;  ///< local dofs of u_old for each block
  std::vector<double> u_values_;          ///< [u_0..u_{n-1}, u_old_0..u_old_{n-1}] at ip
  std::vector<double> aux_values_;        ///< auxiliary fields at ip
  mfem::DenseMatrix dshape_;              ///< physical gradients of the row-block shapes
  mfem::Vector v_;                        ///< velocity at ip
  mfem::Vector adv_;                      ///< v . grad(psi_r)

  void prepare_element(const mfem::Array<const mfem::FiniteElement*>& el,
                       mfem::ElementTransformation& Tr);
  double evaluate_at_ip(const mfem::Array<const mfem::FiniteElement*>& el,
                        const mfem::Array<const mfem::Vector*>& elfun,
                        mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip);
  void compute_advection_shape(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr);
  static const mfem::IntegrationRule& get_integration_rule(
      const mfem::Array<const mfem::FiniteElement*>& el, const mfem::ElementTransformation& Tr);

 protected:
  /// Advection term: residual row 'row', advected variable 'var', velocity of block 'vel'
  struct AdvectionTerm {
    int row;
    int var;
    int vel;
  };

  int sdim_{0};                                     ///< space dimension (number of components)
  std::vector<std::vector<Coefficient>> velocity_;  ///< per block: components (empty if none)
  std::vector<AdvectionTerm> terms_;
  bool needs_aux_{false};  ///< at least one velocity component may depend on auxiliary fields

  std::vector<mfem::ParGridFunction> vaux_gf_;

  void get_coefficients() override;

  /// Define the advection terms (called in init(), after get_coefficients())
  virtual void build_terms();

  /// Velocity of block 'blk' at the current integration point (v is already sized)
  virtual void compute_velocity(int blk, const std::span<const double>& u,
                                const std::span<const double>& aux, mfem::Vector& v);

 public:
  void init() override;
  AdvectionNLFormIntegrator(Geometry geometry, const double time_step,
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

#include "Integrators/AdvectionNLFormIntegrator.tpp"
