/**
 * @file AxiCylindricalCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Built a MFEM Coefficient for AxiCylindricalCoefficient
 * @version 0.1
 * @date 2026-06-16
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
#include "Coefficients/AxiCylindricalCoefficient.hpp"

#include <cmath>
#include <cstdlib>
#include <span>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AxiCylindricalCoefficient::AxiCylindricalCoefficient object
 *
 */
AxiCylindricalCoefficient::AxiCylindricalCoefficient() : transip(3) {}

/**
 * @brief Eval AxiCylindricalCoefficient at the integration point
 *
 * @param T
 * @param ip
 * @return real_t
 */
mfem::real_t AxiCylindricalCoefficient::Eval(mfem::ElementTransformation& T,
                                             const mfem::IntegrationPoint& ip) {
  T.Transform(ip, transip);
  return abs(transip[0]);
}
