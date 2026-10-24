/**
 * @file test_CahnHilliardDirichlet.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit test of a complete coupled problem: Cahn-Hilliard (phi, mu) with a Dirichlet condition on mu only
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
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

/**
 * Small transient Cahn-Hilliard problem, split form (phi, mu), solved with the
 * whole SLOTH chain (Variables, TransientOperator, Problem, Coupling,
 * TimeDiscretization) on 8 x 8 quadrangles, 3 time steps, no output file.
 *
 * Boundary conditions: Neumann for phi everywhere; mu = 0 (Dirichlet) on the
 * left and right sides, Neumann elsewhere.
 *
 * phi starts from 0.5 + small perturbation (spinodal decomposition): it
 * evolves everywhere, including on the boundaries where mu is imposed. This
 * makes the test discriminating: a Jacobian coupling mu to phi on the
 * Dirichlet dofs of mu would move mu away from 0 there as soon as phi changes.
 */
TEST_CASE("Problems: coupled Cahn-Hilliard keeps mu = 0 on its Dirichlet boundary",
          "[Problems][Operators]") {
  using namespace Sloth2D;
  const auto dims = std::make_tuple(8, 8, 1., 1.);
  SPA spatial("InlineSquareWithQuadrangles", 1, 0, dims);

  auto boundaries_phi = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                         Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  auto boundaries_mu = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 0.),
                        Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
  auto bcs_phi = BCS(&spatial, boundaries_phi);
  auto bcs_mu = BCS(&spatial, boundaries_mu);

  // Same model parameters as tests/CahnHilliard/2D/test4
  const double epsilon = 0.1, sigma = 1., mob = 1.;
  const double lambda = 1.5 * sigma * epsilon;
  const double omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                           Parameter("lambda", lambda), Parameter("omega", omega),
                           Parameter("power", 2.0));
  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);

  auto phi0 = AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& x, double) {
        return 0.5 + 0.05 * std::cos(3. * M_PI * x[0]) * std::cos(2. * M_PI * x[1]);
      }));
  auto v1 = VAR(&spatial, bcs_phi, std::string("phi"), Glossary::PhaseField, 2, phi0);
  auto v2 = VAR(&spatial, bcs_mu, std::string("mu"), Glossary::ChemicalPotential, 2, 0.);
  auto vars = VARS(v1, v2);

  Coefficients coefs(double_well, capillary, mobility, grad_energy);
  std::vector<SPA*> spatials{&spatial, &spatial};
  TransientOPE oper(spatials, {"CahnHilliard"}, params, TimeScheme::EulerImplicit,
                    "SplitTimeDerivative");
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", std::string("Newton")),
                                     Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12),
                                     Parameter("iter_max", 50)));
  oper.overload_solver(HypreSolverType::HYPRE_GMRES);
  oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU);

  TransientPB problem(oper, vars, {coefs, coefs});
  auto coupling = Coupling("CahnHilliard", problem);
  auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                            Parameter("final_time", 3.e-3),
                                            Parameter("time_step", 1.e-3)),
                                 coupling);

  const mfem::Vector phi_initial = vars.getVariables()[0].get_last();
  time.solve();
  const mfem::Vector phi_final = vars.getVariables()[0].get_last();
  const mfem::Vector mu_final = vars.getVariables()[1].get_last();

  // Dirichlet dofs of mu (true dofs, local to this rank)
  const mfem::Array<int> mu_dirichlet = bcs_mu.GetEssentialDofs();

  // The test is only meaningful if phi has changed on these dofs
  double local_change = 0.;
  for (const int d : mu_dirichlet) {
    local_change = std::max(local_change, std::abs(phi_final(d) - phi_initial(d)));
  }
  double change = 0.;
  MPI_Allreduce(&local_change, &change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  REQUIRE(change > 1.e-6);

  // mu must still be 0 on its Dirichlet dofs
  for (const int d : mu_dirichlet) {
    CAPTURE(d);
    CHECK(std::abs(mu_final(d)) < 1.e-12);
  }
}
