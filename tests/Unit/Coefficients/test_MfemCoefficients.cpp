/**
 * @file test_MfemCoefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the coefficients evaluated on a mesh (MfemCoefficient, MfemFileCoefficient,
 * AxiCylindricalCoefficient)
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
#include <cmath>
#include <functional>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Coefficients/Coefficient.hpp"
#include "Coefficients/FileCoefficient.hpp"
#include "Coefficients/MfemCoefficient.hpp"
#include "Coefficients/MfemFileCoefficient.hpp"
#include "mfem.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
constexpr double tol = 1.e-12;

/**
 * @brief Small parallel mesh built in memory: unit square, 4 x 4 quadrangles,
 * distributed over all the MPI ranks, with an H1 (Q1) finite element space.
 *
 * Used with TEST_CASE_METHOD: each test case gets a new MeshFixture, built
 * before the test and destroyed after it.
 *
 * Q1 elements interpolate exactly the affine functions (a + b x + c y): every
 * expected field below is affine, so the projections are exact (to round-off).
 */
struct MeshFixture {
  mfem::H1_FECollection fec{1, 2};
  std::unique_ptr<mfem::ParMesh> pmesh;
  std::unique_ptr<mfem::ParFiniteElementSpace> fes;

  MeshFixture() {
    auto mesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL, true, 1., 1.);
    pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
    fes = std::make_unique<mfem::ParFiniteElementSpace>(pmesh.get(), &fec);
  }

  /// @brief Grid function interpolating f
  mfem::ParGridFunction project(const std::function<double(const mfem::Vector&)>& f) const {
    mfem::ParGridFunction gf(fes.get());
    mfem::FunctionCoefficient c(f);
    gf.ProjectCoefficient(c);
    return gf;
  }

  /**
   * @brief Checks coef.Eval at the centre of every local element against
   * expected(x), x being the physical position of the centre.
   */
  void check_eval(mfem::Coefficient& coef,
                  const std::function<double(const mfem::Vector&)>& expected) const {
    mfem::Vector x;
    for (int e = 0; e < pmesh->GetNE(); ++e) {
      mfem::ElementTransformation* T = pmesh->GetElementTransformation(e);
      const mfem::IntegrationPoint& ip =
          mfem::Geometries.GetCenter(pmesh->GetElementBaseGeometry(e));
      T->SetIntPoint(&ip);
      T->Transform(ip, x);
      CAPTURE(e, x[0], x[1]);
      CHECK_THAT(coef.Eval(*T, ip), WithinAbs(expected(x), tol));
    }
  }

  /// @brief Global (all ranks) L2 distance between gf and f: collective
  double l2_error(const mfem::ParGridFunction& gf,
                  const std::function<double(const mfem::Vector&)>& f) const {
    mfem::FunctionCoefficient exact(f);
    return gf.ComputeL2Error(exact);
  }
};

/**
 * @brief Test coefficient that shows which values it receives:
 *   F = u + 10 u_explicit + 100 aux
 * (each term only if the corresponding array is not empty).
 */
class ProbeCoefficient final : public FunctionCoefficient {
 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [](const std::span<const double>& u, const std::span<const double>& u_exp,
              const std::span<const double>& aux, const unsigned int) {
      double f = u.empty() ? 0. : u[0];
      if (!u_exp.empty()) f += 10. * u_exp[0];
      if (!aux.empty()) f += 100. * aux[0];
      return f;
    };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  GradientF() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&,
              const unsigned int) { return std::vector<double>(u.size(), 0.); };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  HessianF() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&,
              const unsigned int) { return std::vector<double>(u.size() * u.size(), 0.); };
  }
};

double x_of(const mfem::Vector& p) { return p[0]; }
double y_of(const mfem::Vector& p) { return p[1]; }
double one(const mfem::Vector&) { return 1.; }
}  // namespace

// -----------------------------------------------------------------------------
// MfemCoefficient: a SLOTH Coefficient evaluated at the integration points,
// from the grid functions of the current time (u), previous time (un) and
// auxiliary variables (aux). Here u = x, un = y, aux = 1.
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(MeshFixture, "MfemCoefficient: values used by each time scheme",
                 "[Coefficients][mesh]") {
  std::vector<mfem::ParGridFunction> u{project(x_of)};
  std::vector<mfem::ParGridFunction> un{project(y_of)};
  std::vector<mfem::ParGridFunction> aux{project(one)};

  SECTION("Implicit: F(u, aux) = x + 100") {
    MfemCoefficient c(Coefficient(Glossary::Mobility, Scheme::Implicit, ProbeCoefficient()), u, un,
                      aux);
    auto expected = [](const mfem::Vector& p) { return p[0] + 100.; };
    check_eval(c, expected);
    CHECK(l2_error(c.ProjectToGridFunction(), expected) < 1.e-10);
  }
  SECTION("Explicit: F(un, aux) = y + 100") {
    MfemCoefficient c(Coefficient(Glossary::Mobility, Scheme::Explicit, ProbeCoefficient()), u, un,
                      aux);
    check_eval(c, [](const mfem::Vector& p) { return p[1] + 100.; });
  }
  SECTION("SemiImplicit: F(u, un, aux) = x + 10 y + 100") {
    MfemCoefficient c(Coefficient(Glossary::Mobility, Scheme::SemiImplicit, ProbeCoefficient()), u,
                      un, aux);
    check_eval(c, [](const mfem::Vector& p) { return p[0] + 10. * p[1] + 100.; });
  }
  SECTION("Constant: the grid functions are not used") {
    MfemCoefficient c(Coefficient(Glossary::Mobility, 5.), u, un, aux);
    check_eval(c, [](const mfem::Vector&) { return 5.; });
    CHECK(l2_error(c.ProjectToGridFunction(*fes), [](const mfem::Vector&) { return 5.; }) < 1.e-10);
  }
  SECTION("current and previous times must have the same number of grid functions") {
    std::vector<mfem::ParGridFunction> none;
    SLOTH_CHECK_MFEM_ERROR_WITH(MfemCoefficient(Coefficient(Glossary::Mobility, 5.), u, none, aux),
                                "consistency");
  }
}

// -----------------------------------------------------------------------------
// MfemFileCoefficient: a CSV field evaluated at the physical points
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(MeshFixture, "MfemFileCoefficient: CSV field on the mesh",
                 "[Coefficients][mesh][file]") {
  // F = t + 2x + 3y on a grid larger than the mesh
  std::string content = "t,x,y,F\n";
  for (const double t : {0., 1.}) {
    for (const double x : {0., 0.5, 1.}) {
      for (const double y : {0., 1.}) {
        content += std::to_string(t) + "," + std::to_string(x) + "," + std::to_string(y) + "," +
                   std::to_string(t + 2. * x + 3. * y) + "\n";
      }
    }
  }
  const sloth_unit::TemporaryFile csv("mesh_field.csv", content);
  Coefficient file_coef{FileCoefficient(csv.path())};
  file_coef.set_time(0.5);
  auto expected = [](const mfem::Vector& p) { return 0.5 + 2. * p[0] + 3. * p[1]; };

  SECTION("Eval and projection on a given space") {
    MfemFileCoefficient c(file_coef);
    check_eval(c, expected);
    CHECK(l2_error(c.ProjectToGridFunction(*fes), expected) < 1.e-10);
  }
  SECTION("projection on the space of the grid functions given at construction") {
    std::vector<mfem::ParGridFunction> u{project(x_of)};
    MfemFileCoefficient c(file_coef, u);
    CHECK(l2_error(c.ProjectToGridFunction(), expected) < 1.e-10);
  }
  SECTION("projection without grid function nor space is an error") {
    MfemFileCoefficient c(file_coef);
    SLOTH_CHECK_MFEM_ERROR_WITH(c.ProjectToGridFunction(), "no grid function");
  }
  SECTION("the wrapped coefficient must be read from a file") {
    SLOTH_CHECK_MFEM_ERROR_WITH(MfemFileCoefficient(Coefficient(Glossary::Mobility, 5.)),
                                "not file-based");
  }
}

// -----------------------------------------------------------------------------
// AxiCylindricalCoefficient: radius r = |x| (weight of axisymmetric integrals)
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(MeshFixture, "AxiCylindricalCoefficient: r = |x|", "[Coefficients][mesh]") {
  AxiCylindricalCoefficient r;
  // Element centres have non-integer x (0.125, 0.375...): an integer abs()
  // would truncate them to 0.
  check_eval(r, [](const mfem::Vector& p) { return std::abs(p[0]); });

  mfem::ParGridFunction gf(fes.get());
  gf.ProjectCoefficient(r);
  CHECK(l2_error(gf, [](const mfem::Vector& p) { return std::abs(p[0]); }) < 1.e-10);
}
