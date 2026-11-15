/**
 * @file test_SlothParBlockNonlinearForm.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests of SlothParBlockNonlinearForm: same assembly as mfem::ParBlockNonlinearForm,
 *        integrators prepared before each assembly
 * @version 0.1
 * @date 2026-11-14
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
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

namespace {
using VARS2D = Sloth2D::VARS;

/// @brief Fourier integrator recording the calls of prepare_quadrature_points
class PreparedFourier : public FourierNLFormIntegrator<VARS2D> {
 public:
  PreparedFourier(const std::vector<mfem::ParGridFunction>& u_old, const bool uses)
      : FourierNLFormIntegrator<VARS2D>(Geometry::Cartesian, 0.1, u_old, {},
                                        Parameters(std::vector<Parameter>{}), {},
                                        {Coefficients(Coefficient(Glossary::Conductivity, 2.))}),
        uses_(uses) {}

  bool uses_quadrature_points() const override { return uses_; }
  void prepare_quadrature_points(const mfem::BlockVector& x) override {
    ++calls;
    prepared = x;
  }

  int calls = 0;
  mfem::Vector prepared;  ///< values received by the last call

 private:
  bool uses_;
};

/// @brief Largest absolute difference between two vectors (all ranks)
double max_difference(const mfem::Vector& a, const mfem::Vector& b) {
  REQUIRE(a.Size() == b.Size());
  double local = 0., global = 0.;
  for (int i = 0; i < a.Size(); ++i) local = std::max(local, std::abs(a(i) - b(i)));
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

double phi_old(const mfem::Vector& p) { return 0.5 + 0.2 * std::sin(3. * p[0]); }
}  // namespace

// Residual and Jacobian identical (same operations) to those of mfem::ParBlockNonlinearForm
TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "SlothParBlockNonlinearForm: same residual and Jacobian as "
                 "mfem::ParBlockNonlinearForm",
                 "[Operators]") {
  const bool uses = GENERATE(false, true);
  INFO("uses_quadrature_points: " << uses);
  const std::vector<mfem::ParGridFunction> u_old{project(phi_old)};
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes.get();

  mfem::ParBlockNonlinearForm reference(spaces);
  auto* reference_integrator = new PreparedFourier(u_old, uses);
  reference_integrator->init();
  reference.AddDomainIntegrator(reference_integrator);

  SlothParBlockNonlinearForm<VARS2D> form(spaces);
  auto* integrator = new PreparedFourier(u_old, uses);
  integrator->init();
  form.AddDomainIntegrator(integrator);

  mfem::Vector x(fes->GetTrueVSize()), v(fes->GetTrueVSize());
  x.Randomize(5);
  v.Randomize(7);
  mfem::Vector y_reference(x.Size()), y(x.Size());
  reference.Mult(x, y_reference);
  form.Mult(x, y);
  CHECK(max_difference(y_reference, y) == 0.);

  mfem::Vector jv_reference(x.Size()), jv(x.Size());
  reference.GetGradient(x).Mult(v, jv_reference);
  form.GetGradient(x).Mult(v, jv);
  CHECK(max_difference(jv_reference, jv) == 0.);
}

// prepare_quadrature_points is called before each assembly with the local dofs of x, only for
// the integrators using the quadrature points
TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "SlothParBlockNonlinearForm: integrators prepared before each assembly",
                 "[Operators]") {
  const std::vector<mfem::ParGridFunction> u_old{project(phi_old)};
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes.get();
  SlothParBlockNonlinearForm<VARS2D> form(spaces);
  auto* prepared = new PreparedFourier(u_old, true);
  auto* not_prepared = new PreparedFourier(u_old, false);
  prepared->init();
  not_prepared->init();
  form.AddDomainIntegrator(prepared);
  form.AddDomainIntegrator(not_prepared);

  mfem::Vector x(fes->GetTrueVSize()), y(fes->GetTrueVSize());
  x.Randomize(5);
  mfem::ParGridFunction x_local(fes.get());
  x_local.SetFromTrueDofs(x);

  form.Mult(x, y);
  CHECK(prepared->calls == 1);
  CHECK(max_difference(prepared->prepared, x_local) == 0.);

  x *= 0.5;
  x_local.SetFromTrueDofs(x);
  form.GetGradient(x);
  CHECK(prepared->calls == 2);
  CHECK(max_difference(prepared->prepared, x_local) == 0.);

  CHECK(not_prepared->calls == 0);
}
