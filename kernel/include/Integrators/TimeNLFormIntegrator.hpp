/**
 * @file TimeNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative
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
#include <span>
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
 * @brief Class dedicated to the VF of the time-derivative (eg. LHS of Allen-Cahn  equation)
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class TimeNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  /// Base class of the nonlinear integrators of SLOTH
  using SlothNLBase = SlothNLFormIntegrator<VARS>;

  mfem::Vector Psi;
  /// Coefficients A and B at the integration points of the current element, at the end of the
  /// step (see SlothNLFormIntegrator::eval_coefficient_at_points)
  typename SlothNLBase::CoefficientAtPoints a_at_points_, b_at_points_;

 protected:
  /// Coefficient A of block blk at all the integration points of the element, at the end of the
  /// step (with its first derivatives for the Jacobian); see the derived classes for other forms
  virtual void eval_a_at_points(const int blk, const mfem::Array<const mfem::FiniteElement*>& el,
                                mfem::ElementTransformation& Tr, const mfem::IntegrationRule& ir,
                                const mfem::Array<const mfem::Vector*>& elfun, const bool jacobian);
  /// Value of A at the integration point i (see eval_a_at_points)
  virtual double a_value(const int i) const;
  /// Derivative of A with respect to block jblk at the integration point i
  virtual double a_derivative(const int i, const int jblk) const;
  /// Other terms of the Jacobian of (A B k, psi) due to A at the integration point i (none by
  /// default), multiplied by factor = dt B k w
  virtual void add_a_jacobian(const int /*blk*/, const int /*i*/,
                              const mfem::Array<const mfem::FiniteElement*>& /*el*/,
                              mfem::ElementTransformation& /*Tr*/, const double /*factor*/,
                              const mfem::Vector& /*psi*/,
                              const mfem::Array2D<mfem::DenseMatrix*>& /*elmats*/) {}

  std::vector<mfem::ParGridFunction> vaux_gf_;
  std::list<GlossaryType> expected_list_;
  Coefficients coefficient_A;
  Coefficients coefficient_B;
  void get_coefficients() override;

 public:
  void init() override;
  TimeNLFormIntegrator(Geometry geometry, const double time_step,
                       const std::vector<mfem::ParGridFunction> u_old,
                       const std::vector<mfem::ParGridFunction> aux_old, const Parameters& params,
                       std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
};

#include "Integrators/TimeNLFormIntegrator.tpp"
