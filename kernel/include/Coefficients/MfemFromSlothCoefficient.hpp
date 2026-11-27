/**
 * @file MfemFromSlothCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Wraps a SLOTH Coefficient, possibly depending on other variables, as an mfem::Coefficient
 * @version 0.1
 * @date 2026-11-26
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
#include <span>
#include <string>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Wraps a SLOTH ::Coefficient as an mfem::Coefficient, used to define the
 *        initial condition or the analytical solution of a Variable.
 *
 * @details Depending on the wrapped coefficient:
 *          - scalar: the constant value is returned;
 *          - file-based (Scheme::File): the physical point is forwarded as the
 *            first argument of ::Coefficient::compute() (same convention as
 *            MfemFileCoefficient);
 *          - otherwise: the values of the dependency grid functions at the
 *            integration point are forwarded as the first argument of
 *            ::Coefficient::compute(), in the order given to bind(). The physical
 *            point is passed through ::Coefficient::set_coordinates() when the
 *            coefficient uses coordinates.
 *
 *          Dependency grid functions are only stored as non-owning pointers, and
 *          only for the duration of a projection: call bind() right before the
 *          projection and unbind() right after it.
 */
class MfemFromSlothCoefficient : public mfem::Coefficient {
 private:
  ::Coefficient sloth_coefficient_;
  std::vector<const mfem::ParGridFunction*> dependencies_;

  mfem::Vector x_;                    ///< Physical point, reused across Eval calls
  std::vector<double> values_at_ip_;  ///< Values of the dependencies, reused across Eval calls

 public:
  explicit MfemFromSlothCoefficient(::Coefficient coef);

  bool is_scalar() const;
  bool is_file() const;
  bool requires_dependencies() const;

  void bind(std::vector<const mfem::ParGridFunction*> dependencies);
  void unbind();
  std::size_t get_dependencies_number() const;

  void SetTime(mfem::real_t t) override;
  mfem::real_t Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;

  ~MfemFromSlothCoefficient() override = default;
};
