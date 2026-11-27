/**
 * @file FieldsAtPoints.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Values and physical gradients of scalar fields at the integration points of an element
 * @version 0.1
 * @date 2026-11-22
 *
 * Copyright CEA (C) 2026
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

#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Values and physical gradients of scalar fields at the integration points of an element.
 *
 * The fields share the finite element fe: element_dofs[f] are the dofs of the field f on the
 * element. Layout of the outputs (inputs of eval_variational_at_points, see
 * Coefficients/VariationalEval.hpp), with nf fields, nip points and dim = Tr.GetSpaceDim():
 *   values(f, q)       = values[q * nf + f]
 *   gradients(f, k, q) = gradients[(q * nf + f) * dim + k]
 *
 * @param fe Finite element shared by the fields
 * @param Tr Transformation of the element
 * @param ir Integration rule
 * @param element_dofs Dofs of each field on the element
 * @param values Values of the fields at the points (resized to nip * nf)
 * @param gradients Physical gradients of the fields at the points (resized to nip * nf * dim),
 *        or nullptr if they are not needed
 */
inline void fields_at_points(const mfem::FiniteElement& fe, mfem::ElementTransformation& Tr,
                             const mfem::IntegrationRule& ir,
                             const std::vector<mfem::Vector>& element_dofs, mfem::Vector& values,
                             mfem::Vector* gradients) {
  const int nip = ir.GetNPoints();
  const int nf = static_cast<int>(element_dofs.size());
  const int nd = fe.GetDof();
  const int dim = Tr.GetSpaceDim();
  mfem::Vector shape(nd);
  mfem::DenseMatrix dshape(nd, dim);

  values.SetSize(nip * nf);
  double* v = values.HostWrite();
  double* g = nullptr;
  if (gradients != nullptr) {
    gradients->SetSize(nip * nf * dim);
    g = gradients->HostWrite();
  }
  for (int q = 0; q < nip; ++q) {
    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
    Tr.SetIntPoint(&ip);
    fe.CalcPhysShape(Tr, shape);
    if (g != nullptr) {
      fe.CalcPhysDShape(Tr, dshape);
    }
    for (int f = 0; f < nf; ++f) {
      const mfem::Vector& dofs = element_dofs[static_cast<std::size_t>(f)];
      v[q * nf + f] = shape * dofs;
      if (g == nullptr) continue;
      for (int k = 0; k < dim; ++k) {
        double sum = 0.;
        for (int j = 0; j < nd; ++j) sum += dshape(j, k) * dofs(j);
        g[(q * nf + f) * dim + k] = sum;
      }
    }
  }
}
