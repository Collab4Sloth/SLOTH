/**
 * @file test_ReducedOperators.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for PhaseFieldReducedOperator and SteadyPhaseFieldReducedOperator
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
#include <cmath>
#include <memory>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Operators/ReducedOperator.hpp"
#include "Operators/SteadyReducedOperator.hpp"
#include "mfem.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
/**
 * @brief Linear block integrator for the tests, for N variables on the same space:
 *   block (i, i) = a K + c M,   block (i, j != i) = b M
 * with K and M the element diffusion and mass matrices.
 */
class LinearBlockIntegrator : public mfem::BlockNonlinearFormIntegrator {
 public:
  LinearBlockIntegrator(double a, double b, double c) : a_(a), b_(b), c_(c) {}

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>&,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override {
    mfem::DenseMatrix K, M;
    diffusion_.AssembleElementMatrix(*el[0], Tr, K);
    mass_.AssembleElementMatrix(*el[0], Tr, M);
    const int n = el.Size();
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        mfem::DenseMatrix& Eij = *elmats(i, j);
        Eij.SetSize(K.Height());
        Eij = 0.;
        if (i == j) {
          Eij.Add(a_, K);
          Eij.Add(c_, M);
        } else {
          Eij.Add(b_, M);
        }
      }
    }
  }

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override {
    const int n = el.Size();
    mfem::Array2D<mfem::DenseMatrix*> elmats(n, n);
    std::vector<mfem::DenseMatrix> storage(static_cast<size_t>(n * n));
    for (int k = 0; k < n * n; ++k) elmats(k / n, k % n) = &storage[static_cast<size_t>(k)];
    AssembleElementGrad(el, Tr, elfun, elmats);
    for (int i = 0; i < n; ++i) {
      elvec[i]->SetSize(el[i]->GetDof());
      *elvec[i] = 0.;
      for (int j = 0; j < n; ++j) elmats(i, j)->AddMult(*elfun[j], *elvec[i]);
    }
  }

 private:
  double a_, b_, c_;
  mfem::DiffusionIntegrator diffusion_;
  mfem::MassIntegrator mass_;
};

/**
 * @brief N variables on the unit square, a block form "LHS" (mass on the
 * diagonal) and a coupled block form "RHS", Dirichlet on the whole boundary
 * for every variable.
 */
struct BlockProblem : sloth_unit::UnitSquare {
  int nvar;
  mfem::Array<mfem::ParFiniteElementSpace*> spaces;
  std::unique_ptr<mfem::ParBlockNonlinearForm> lhs, rhs;
  std::vector<mfem::Array<int>> ess;
  int tsize;

  explicit BlockProblem(int n) : nvar(n) {
    spaces.SetSize(nvar);
    for (int i = 0; i < nvar; ++i) spaces[i] = fes.get();
    lhs = std::make_unique<mfem::ParBlockNonlinearForm>(spaces);
    lhs->AddDomainIntegrator(new LinearBlockIntegrator(0., 0., 1.));  // M on the diagonal
    rhs = std::make_unique<mfem::ParBlockNonlinearForm>(spaces);
    rhs->AddDomainIntegrator(new LinearBlockIntegrator(1., 0.5, 0.));  // K + 0.5 M coupling
    ess.assign(static_cast<size_t>(nvar), ess_tdof);
    tsize = fes->GetTrueVSize();
  }

  mfem::Vector random(int seed) const {
    mfem::Vector v(nvar * tsize);
    v.Randomize(seed);
    return v;
  }
};

/// @brief Marks the rows of the Dirichlet dofs of every variable
std::vector<bool> dirichlet_rows(const BlockProblem& p) {
  std::vector<bool> is_ess(static_cast<size_t>(p.nvar * p.tsize), false);
  for (int i = 0; i < p.nvar; ++i) {
    for (const int d : p.ess[static_cast<size_t>(i)]) {
      is_ess[static_cast<size_t>(i * p.tsize + d)] = true;
    }
  }
  return is_ess;
}

/**
 * @brief Rows of the non-Dirichlet dofs: J v = (F(k + v) - F(k - v)) / 2
 * (exact: F is affine), for v = 0 on the Dirichlet dofs. The Jacobian also
 * eliminates the Dirichlet columns, which is consistent with a Newton
 * increment that vanishes on the Dirichlet dofs.
 */
void check_gradient_free_rows(const mfem::Operator& F, const BlockProblem& p,
                              const mfem::Vector& k) {
  const auto is_ess = dirichlet_rows(p);
  mfem::Vector v = p.random(7);
  for (int r = 0; r < v.Size(); ++r) {
    if (is_ess[static_cast<size_t>(r)]) v(r) = 0.;
  }
  mfem::Vector kp(k), km(k), fp(k.Size()), fm(k.Size()), fd(k.Size()), jv(k.Size());
  kp += v;
  km -= v;
  F.Mult(kp, fp);
  F.Mult(km, fm);
  subtract(0.5, fp, fm, fd);
  F.GetGradient(k).Mult(v, jv);
  for (int r = 0; r < k.Size(); ++r) {
    if (is_ess[static_cast<size_t>(r)]) continue;
    CAPTURE(r);
    CHECK_THAT(jv(r), WithinAbs(fd(r), 1.e-10 * (1. + std::abs(fd(r)))));
  }
}

/**
 * @brief Rows of the Dirichlet dofs: identity rows, (J v)_r = v_r, for any v.
 * A 1 in an off-diagonal block (i, j) would add v of the other variable.
 */
void check_gradient_dirichlet_rows(const mfem::Operator& F, const BlockProblem& p,
                                   const mfem::Vector& k) {
  const auto is_ess = dirichlet_rows(p);
  const mfem::Vector v = p.random(7);
  mfem::Vector jv(k.Size());
  F.GetGradient(k).Mult(v, jv);
  for (int r = 0; r < k.Size(); ++r) {
    if (!is_ess[static_cast<size_t>(r)]) continue;
    CAPTURE(r);
    CHECK_THAT(jv(r), WithinAbs(v(r), 1.e-12));
  }
}
}  // namespace

TEST_CASE("SteadyPhaseFieldReducedOperator", "[Operators]") {
  const int nvar = GENERATE(1, 2);
  CAPTURE(nvar);
  BlockProblem p(nvar);
  SteadyPhaseFieldReducedOperator F(p.rhs.get(), p.ess);
  REQUIRE(F.Height() == nvar * p.tsize);

  SECTION("Mult: residual of RHS, zero on Dirichlet dofs") {
    const mfem::Vector k = p.random(3);
    mfem::Vector y(k.Size()), expected(k.Size());
    F.Mult(k, y);
    p.rhs->Mult(k, expected);
    for (int i = 0; i < nvar; ++i) {
      for (const int d : p.ess[static_cast<size_t>(i)]) expected(i * p.tsize + d) = 0.;
    }
    expected -= y;
    CHECK(expected.Normlinf() < 1.e-12);
  }
  SECTION("GetGradient, non-Dirichlet rows") { check_gradient_free_rows(F, p, p.random(3)); }
}

TEST_CASE("PhaseFieldReducedOperator", "[Operators]") {
  // F(k) = LHS(k) + RHS(u + dt k), Dirichlet rows set to zero
  const int nvar = GENERATE(1, 2);
  CAPTURE(nvar);
  BlockProblem p(nvar);
  PhaseFieldReducedOperator F(p.lhs.get(), p.rhs.get(), p.ess);
  const double dt = 0.1;
  const mfem::Vector u = p.random(1);
  F.SetParameters(dt, &u);

  SECTION("Mult") {
    const mfem::Vector k = p.random(3);
    mfem::Vector y(k.Size()), expected(k.Size()), z(u);
    F.Mult(k, y);
    z.Add(dt, k);
    p.rhs->Mult(z, expected);
    p.lhs->AddMult(k, expected);
    for (int i = 0; i < nvar; ++i) {
      for (const int d : p.ess[static_cast<size_t>(i)]) expected(i * p.tsize + d) = 0.;
    }
    expected -= y;
    CHECK(expected.Normlinf() < 1.e-12);
  }
  SECTION("GetGradient, non-Dirichlet rows") { check_gradient_free_rows(F, p, p.random(3)); }
}

// With several variables, only the diagonal blocks get a 1 on their diagonal
// when the Dirichlet conditions are eliminated: the Dirichlet row of variable
// i reads v_i (a 1 in an off-diagonal block would make it read v_i + v_j).
TEST_CASE("Reduced operators: Dirichlet rows of the Jacobian are identity rows",
          "[Operators]") {
  const int nvar = GENERATE(1, 2);
  CAPTURE(nvar);
  BlockProblem p(nvar);
  const mfem::Vector k = p.random(3);

  SECTION("SteadyPhaseFieldReducedOperator") {
    SteadyPhaseFieldReducedOperator F(p.rhs.get(), p.ess);
    check_gradient_dirichlet_rows(F, p, k);
  }
  SECTION("PhaseFieldReducedOperator") {
    PhaseFieldReducedOperator F(p.lhs.get(), p.rhs.get(), p.ess);
    const mfem::Vector u = p.random(1);
    F.SetParameters(0.1, &u);
    check_gradient_dirichlet_rows(F, p, k);
  }
}

// Same property in a Newton solve. Two coupled variables, Dirichlet on
// variable 0 only: the time derivative k0 stays 0 on the Dirichlet dofs (the
// boundary value of u0 is fixed). A 1 in an off-diagonal block would make the
// Dirichlet row read dk0 + dk1 = 0, hence k0 != 0 on the boundary.
TEST_CASE("PhaseFieldReducedOperator: Newton keeps the Dirichlet values",
          "[Operators]") {
  BlockProblem p(2);
  p.ess[1].DeleteAll();  // no Dirichlet condition on variable 1
  PhaseFieldReducedOperator F(p.lhs.get(), p.rhs.get(), p.ess);
  const mfem::Vector u = p.random(1);
  F.SetParameters(0.1, &u);

  mfem::GMRESSolver gmres(MPI_COMM_WORLD);  // the Jacobian is not symmetric
  gmres.SetRelTol(1.e-12);
  gmres.SetMaxIter(500);
  mfem::NewtonSolver newton(MPI_COMM_WORLD);
  newton.SetSolver(gmres);
  newton.SetOperator(F);
  newton.SetRelTol(1.e-10);
  newton.SetMaxIter(20);

  mfem::Vector k(F.Height());
  k = 0.;
  newton.Mult(mfem::Vector(), k);
  REQUIRE(newton.GetConverged());

  for (const int d : p.ess[0]) {
    CAPTURE(d);
    CHECK_THAT(k(d), WithinAbs(0., 1.e-10));  // k0 on the Dirichlet dofs of variable 0
  }
}
