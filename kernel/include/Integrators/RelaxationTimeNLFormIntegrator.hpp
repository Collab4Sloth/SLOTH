/**
 * @file RelaxationTimeNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative with a relaxation time tau(u, grad u)
 * @version 0.1
 * @date 2026-11-23
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
#include <string>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Integrators/TimeNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Time derivative with a relaxation time (e.g. LHS of the Allen-Cahn equation
 * tau dphi/dt = - dF/dphi):
 *   R_b(psi) = int tau du_b/dt psi,
 * with tau the coefficient of type Glossary::RelaxationTime of each block (one by default).
 * tau can be a constant, a function of the blocks or a variational coefficient tau(u, grad u)
 * ("variational": true in the JSON file) of the first blocks and of their gradients, e.g. an
 * anisotropic relaxation time tau(n), n = grad(phi) / |grad(phi)|. Implicit tau: value at the end
 * of the step u = u_old + dt k and derivatives with respect to the blocks and to their gradients
 * in the Jacobian. Explicit tau: values and gradients of u_old, no contribution to the Jacobian.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class RelaxationTimeNLFormIntegrator : public TimeNLFormIntegrator<VARS> {
 private:
  /// Base class of the nonlinear integrators of SLOTH
  using SlothNLBase = SlothNLFormIntegrator<VARS>;

  /// Variational relaxation time tau(u, grad u) at the integration points of the current element
  struct VariationalTau {
    bool active = false;  ///< tau of the current block variational
    int nb = 0;           ///< number of variables of tau (the first blocks)
    int dim = 0;          ///< dimension of the gradients
    bool with_derivatives = false;
    std::vector<mfem::Vector> dofs, dofs_old;
    mfem::Vector values, explicit_values, auxiliary_values;
    mfem::Vector gradients, explicit_gradients, auxiliary_gradients, coordinates;
    mfem::DenseMatrix points;
    /// tau, dtau/du_j [q nb + j] and dtau/d(grad u_j) [(q nb + j) dim + k] at the points q
    mfem::Vector F, dtau_du, dtau_dgrad;
  };
  VariationalTau tau_;
  mfem::DenseMatrix gradPsi_j_;
  mfem::Vector shape_j_, row_j_;

  void variational_tau_at_points(const Coefficient& tau,
                                 const mfem::Array<const mfem::FiniteElement*>& el,
                                 mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
                                 const mfem::Array<const mfem::Vector*>& elfun,
                                 const bool jacobian);

 protected:
  void get_coefficients() override;
  void eval_a_at_points(const int blk, const mfem::Array<const mfem::FiniteElement*>& el,
                        mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
                        const mfem::Array<const mfem::Vector*>& elfun,
                        const bool jacobian) override;
  double a_value(const int i) const override;
  double a_derivative(const int i, const int jblk) const override;
  void add_a_jacobian(const int blk, const int i, const mfem::Array<const mfem::FiniteElement*>& el,
                      mfem::ElementTransformation& Tr, const double factor, const mfem::Vector& psi,
                      const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;

 public:
  RelaxationTimeNLFormIntegrator(Geometry geometry, const double time_step,
                                 const std::vector<mfem::ParGridFunction> u_old,
                                 const std::vector<mfem::ParGridFunction>& aux_old,
                                 const Parameters& params, std::vector<VARS*> auxvars,
                                 const std::vector<Coefficients>& coefficients);
  virtual ~RelaxationTimeNLFormIntegrator() = default;
};

#include "Integrators/RelaxationTimeNLFormIntegrator.tpp"
