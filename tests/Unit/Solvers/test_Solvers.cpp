/**
 * @file test_Solvers.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the linear and nonlinear solvers (SlothSolver, LSolver, NLSolver)
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
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "Parameters/Parameters.hpp"
#include "Solvers/HSolverBase.hpp"
#include "Solvers/IPrecondBase.hpp"
#include "Solvers/ISolverBase.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Solvers/SlothSolver.hpp"
#include "mfem.hpp"
#include "sloth_unit.hpp"

namespace {
/**
 * @brief Laplace problem on the unit square (4 x 4 quadrangles, Q1, all MPI
 * ranks): -Delta u = 0, u = x on the boundary. The exact solution u = x is
 * affine, hence exactly represented in Q1: every solver must recover it up to
 * its tolerance.
 */
struct PoissonFixture {
  mfem::H1_FECollection fec{1, 2};
  std::unique_ptr<mfem::ParMesh> pmesh;
  std::unique_ptr<mfem::ParFiniteElementSpace> fes;
  mfem::Array<int> ess_bdr, ess_tdof;
  mfem::FunctionCoefficient exact{[](const mfem::Vector& p) { return p[0]; }};

  PoissonFixture() {
    auto mesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL, true, 1., 1.);
    pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
    fes = std::make_unique<mfem::ParFiniteElementSpace>(pmesh.get(), &fec);
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  // Dirichlet on the whole boundary
    fes->GetEssentialTrueDofs(ess_bdr, ess_tdof);
  }
};

/// @brief Applies the solver held by the variant (nothing for "no solver")
void apply(VSharedMFEMSolver& solver, const mfem::Vector& b, mfem::Vector& x) {
  std::visit(
      [&](auto& s) {
        if constexpr (!std::is_same_v<std::decay_t<decltype(s)>, sptr<std::monostate>>) {
          s->Mult(b, x);
        }
      },
      solver);
}

/// @brief Convergence flag of an MFEM iterative solver held by the variant
bool converged(VSharedMFEMSolver& solver) {
  bool result = false;
  std::visit(
      [&](auto& s) {
        using S = typename std::decay_t<decltype(s)>::element_type;
        if constexpr (std::is_base_of_v<mfem::IterativeSolver, S>) result = s->GetConverged();
      },
      solver);
  return result;
}

/// @brief A solver / preconditioner pair, by its SLOTH enums
struct Combo {
  std::string name;
  VSolverType solver;
  VSolverType precond;
};

/// @brief Valid solver / preconditioner pairs
std::vector<Combo> all_combos() {
  std::vector<Combo> combos{
      {"CG", IterativeSolverType::CG, PreconditionerType::NO},
      {"CG + HypreSmoother", IterativeSolverType::CG, HyprePreconditionerType::HYPRE_SMOOTHER},
      {"BiCGSTAB + HypreSmoother", IterativeSolverType::BICGSTAB,
       HyprePreconditionerType::HYPRE_SMOOTHER},
      {"GMRES", IterativeSolverType::GMRES, PreconditionerType::NO},
      {"MINRES", IterativeSolverType::MINRES, PreconditionerType::NO},
      {"HyprePCG", HypreSolverType::HYPRE_PCG, HyprePreconditionerType::NO},
      {"HyprePCG + BoomerAMG", HypreSolverType::HYPRE_PCG,
       HyprePreconditionerType::HYPRE_BOOMER_AMG},
      {"HyprePCG + DiagScale", HypreSolverType::HYPRE_PCG,
       HyprePreconditionerType::HYPRE_DIAG_SCALE},
      {"HypreGMRES + ILU", HypreSolverType::HYPRE_GMRES, HyprePreconditionerType::HYPRE_ILU},
      {"HypreFGMRES + BoomerAMG", HypreSolverType::HYPRE_FGMRES,
       HyprePreconditionerType::HYPRE_BOOMER_AMG}};
#ifdef SLOTH_USE_MUMPS
  combos.push_back({"MUMPS", DirectSolverType::MUMPS, PreconditionerType::NO});
#endif
  return combos;
}
}  // namespace

// -----------------------------------------------------------------------------
// SlothSolver: enum -> MFEM object
// -----------------------------------------------------------------------------
TEST_CASE("SlothSolver: each enum value creates the corresponding MFEM object", "[Solvers]") {
  const Parameters params(std::vector<Parameter>{});
  auto make = [&](VSolverType type) { return SlothSolver(type, params).get_value(); };

  CHECK(std::holds_alternative<sptr<mfem::CGSolver>>(make(IterativeSolverType::CG)));
  CHECK(std::holds_alternative<sptr<mfem::GMRESSolver>>(make(IterativeSolverType::GMRES)));
  CHECK(std::holds_alternative<sptr<mfem::BiCGSTABSolver>>(make(IterativeSolverType::BICGSTAB)));
  CHECK(std::holds_alternative<sptr<mfem::MINRESSolver>>(make(IterativeSolverType::MINRES)));

  CHECK(std::holds_alternative<sptr<mfem::HyprePCG>>(make(HypreSolverType::HYPRE_PCG)));
  CHECK(std::holds_alternative<sptr<mfem::HypreGMRES>>(make(HypreSolverType::HYPRE_GMRES)));
  CHECK(std::holds_alternative<sptr<mfem::HypreFGMRES>>(make(HypreSolverType::HYPRE_FGMRES)));

  CHECK(std::holds_alternative<sptr<mfem::NewtonSolver>>(make(NLSolverType::NEWTON)));
  CHECK(std::holds_alternative<sptr<mfem::LBFGSSolver>>(make(NLSolverType::LBFGS)));

  CHECK(std::holds_alternative<sptr<mfem::DSmoother>>(make(PreconditionerType::SMOOTHER)));
  CHECK(std::holds_alternative<sptr<std::monostate>>(make(PreconditionerType::NO)));

  CHECK(std::holds_alternative<sptr<mfem::HypreILU>>(make(HyprePreconditionerType::HYPRE_ILU)));
  CHECK(std::holds_alternative<sptr<mfem::HypreBoomerAMG>>(
      make(HyprePreconditionerType::HYPRE_BOOMER_AMG)));
  CHECK(std::holds_alternative<sptr<mfem::HypreDiagScale>>(
      make(HyprePreconditionerType::HYPRE_DIAG_SCALE)));
  CHECK(std::holds_alternative<sptr<mfem::HypreSmoother>>(
      make(HyprePreconditionerType::HYPRE_SMOOTHER)));
  CHECK(std::holds_alternative<sptr<std::monostate>>(make(HyprePreconditionerType::NO)));

#ifdef SLOTH_USE_MUMPS
  CHECK(std::holds_alternative<sptr<mfem::MUMPSSolver>>(make(DirectSolverType::MUMPS)));
#else
  SLOTH_CHECK_MFEM_ERROR_WITH(make(DirectSolverType::MUMPS), "DirectSolverType");
#endif
}

TEST_CASE("SolverBase: description", "[Solvers]") {
  SolverCG cg;
  cg.create_solver(Parameters(Parameter("description", std::string("my CG"))));
  CHECK(cg.get_description() == "my CG");

  SECTION("the description is required") {
    SolverCG other;
    SLOTH_CHECK_MFEM_ERROR_WITH(other.create_solver(Parameters(std::vector<Parameter>{})),
                                "Not Found");
  }
}

// -----------------------------------------------------------------------------
// LSolver: every valid solver / preconditioner pair solves the problem
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(PoissonFixture, "LSolver: solver / preconditioner pairs", "[Solvers]") {
  const Combo combo = GENERATE(Catch::Generators::from_range(all_combos()));
  INFO(combo.name);  // printed if a check fails

  // Linear system A X = B, Dirichlet values u = x on the boundary
  mfem::ParBilinearForm a(fes.get());
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.Assemble();
  mfem::ParLinearForm b(fes.get());
  b.Assemble();  // zero right-hand side
  mfem::ParGridFunction u(fes.get());
  u = 0.;
  u.ProjectBdrCoefficient(exact, ess_bdr);
  mfem::HypreParMatrix A;
  mfem::Vector B, X;
  a.FormLinearSystem(ess_tdof, u, b, A, X, B);

  const Parameters no_params(std::vector<Parameter>{});
  LSolver lsolver(combo.solver, no_params, combo.precond, no_params, A);
  auto solver = lsolver.get_solver();
  apply(solver, B, X);
  a.RecoverFEMSolution(X, b, u);

  CHECK(u.ComputeL2Error(exact) < 1.e-8);
}

TEST_CASE_METHOD(PoissonFixture, "LSolver: solver parameters are used", "[Solvers]") {
  mfem::ParBilinearForm a(fes.get());
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.Assemble();
  mfem::ParLinearForm b(fes.get());
  b.Assemble();
  mfem::ParGridFunction u(fes.get());
  u = 0.;
  u.ProjectBdrCoefficient(exact, ess_bdr);
  mfem::HypreParMatrix A;
  mfem::Vector B, X;
  a.FormLinearSystem(ess_tdof, u, b, A, X, B);

  SECTION("default parameters: converged (constructor without preconditioner)") {
    LSolver lsolver(IterativeSolverType::CG, Parameters(std::vector<Parameter>{}), A);
    auto solver = lsolver.get_solver();
    apply(solver, B, X);
    CHECK(converged(solver));
  }
  SECTION("iter_max = 1: not converged") {
    LSolver lsolver(IterativeSolverType::CG, Parameters(Parameter("iter_max", 1)), A);
    auto solver = lsolver.get_solver();
    apply(solver, B, X);
    CHECK_FALSE(converged(solver));
  }
}

TEST_CASE_METHOD(PoissonFixture, "LSolver: incompatible pairs are rejected", "[Solvers]") {
  mfem::ParBilinearForm a(fes.get());
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  const Parameters no_params(std::vector<Parameter>{});

  SECTION("MFEM iterative solver with a Hypre preconditioner") {
    SLOTH_CHECK_MFEM_ERROR_WITH(LSolver(IterativeSolverType::CG, no_params,
                                        HyprePreconditionerType::HYPRE_BOOMER_AMG, no_params, *A),
                                "can only be associated");
  }
  SECTION("Hypre solver with an MFEM preconditioner") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        LSolver(HypreSolverType::HYPRE_PCG, no_params, PreconditionerType::SMOOTHER, no_params, *A),
        "can only be associated");
  }
}

// -----------------------------------------------------------------------------
// NLSolver: Newton on the same problem written as a nonlinear form
// (linear problem: Newton converges in one iteration, up to the linear solver)
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(PoissonFixture, "NLSolver: Newton", "[Solvers]") {
  mfem::ParNonlinearForm n(fes.get());
  n.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  n.SetEssentialTrueDofs(ess_tdof);

  // Initial guess: exact values on the boundary, 0 inside
  mfem::ParGridFunction u_exact(fes.get());
  u_exact.ProjectCoefficient(exact);
  mfem::Vector X_exact;
  u_exact.GetTrueDofs(X_exact);
  mfem::Vector X(fes->GetTrueVSize());
  X = 0.;
  for (const int i : ess_tdof) X(i) = X_exact(i);

  // iterative_mode: Newton starts from X (and keeps the boundary values)
  const Parameters nl_params(Parameter("iterative_mode", true));
  const Parameters no_params(std::vector<Parameter>{});

  SECTION("with CG") {
    NLSolver nl(NLSolverType::NEWTON, nl_params, IterativeSolverType::CG, no_params, n);
    nl.get_nl_solver()->Mult(mfem::Vector(), X);
    CHECK(nl.get_nl_solver()->GetConverged());
  }
  SECTION("with HyprePCG + BoomerAMG") {
    NLSolver nl(NLSolverType::NEWTON, nl_params, HypreSolverType::HYPRE_PCG, no_params,
                HyprePreconditionerType::HYPRE_BOOMER_AMG, no_params, n);
    nl.get_nl_solver()->Mult(mfem::Vector(), X);
    CHECK(nl.get_nl_solver()->GetConverged());
  }

  X -= X_exact;
  CHECK(X.Normlinf() < 1.e-8);  // on this rank
}
