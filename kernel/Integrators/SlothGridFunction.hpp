/**
 * @file SlothGridFunction.hpp
 * @author clement.plumecocq@cea.fr
 * @brief Class dedicated to the calculation of the gradient of a GridFunction
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
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

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class SlothGridFunction : public mfem::ParGridFunction {
 private:
 protected:
 public:
  SlothGridFunction();
  explicit SlothGridFunction(const mfem::ParGridFunction& orig);
  void GetGradient(mfem::ElementTransformation& Tr, mfem::DenseMatrix& gradPsi, mfem::Vector& grad);
  ~SlothGridFunction();
};

/**
 * @brief Construct a new Sloth Grid Function:: Sloth Grid Function object
 *
 */
SlothGridFunction::SlothGridFunction() : mfem::ParGridFunction() {}

/**
 * @brief Construct a new Sloth Grid Function:: Sloth Grid Function object
 *
 * @param gf
 */
SlothGridFunction::SlothGridFunction(const mfem::ParGridFunction& gf) : mfem::ParGridFunction(gf) {}

/**
 * @brief Compute the gradient of a given GridFunction
 *
 * @param Tr
 * @param gradPsi
 * @param grad
 */
void SlothGridFunction::GetGradient(mfem::ElementTransformation& Tr, mfem::DenseMatrix& gradPsi,
                                    mfem::Vector& grad) {
  mfem::Vector elfun;
  const int& num = Tr.ElementNo;
  GetElementDofValues(num, elfun);
  gradPsi.MultTranspose(elfun, grad);
}

/**
 * @brief Destroy the Sloth Grid Function:: Sloth Grid Function object
 *
 */
SlothGridFunction::~SlothGridFunction() {}
