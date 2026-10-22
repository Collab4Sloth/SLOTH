/**
 * @file test_SlothGridFunction.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for SlothGridFunction
 * @version 0.1
 * @date 2026-10-22
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
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Integrators/SlothGridFunction.hpp"
#include "mfem.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
double affine(const mfem::Vector& p) { return 1. + 2. * p[0] - 3. * p[1]; }
}  // namespace

TEST_CASE_METHOD(sloth_unit::UnitSquare, "SlothGridFunction: gradient from shape derivatives",
                 "[Integrators]") {
  CHECK_NOTHROW(SlothGridFunction());  // default constructor

  // grad(1 + 2x - 3y) = (2, -3) everywhere (exact in Q1)
  SlothGridFunction u(project(affine));
  // grad must be sized by the caller: SlothGridFunction::GetGradient does not
  // resize it (unlike mfem::GridFunction::GetGradient) and crashes otherwise
  mfem::Vector grad(2);
  for (int e = 0; e < pmesh->GetNE(); ++e) {
    const mfem::FiniteElement* fe = fes->GetFE(e);
    mfem::ElementTransformation* T = fes->GetElementTransformation(e);
    const mfem::IntegrationPoint& ip = mfem::Geometries.GetCenter(fe->GetGeomType());
    T->SetIntPoint(&ip);
    mfem::DenseMatrix dshape(fe->GetDof(), 2);
    fe->CalcPhysDShape(*T, dshape);  // derivatives of the shape functions
    u.GetGradient(*T, dshape, grad);
    CAPTURE(e);
    CHECK_THAT(grad(0), WithinAbs(2., 1.e-12));
    CHECK_THAT(grad(1), WithinAbs(-3., 1.e-12));
  }
}
