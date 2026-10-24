/**
 * @file test_AMR.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the adaptive mesh refinement (SingleVariableAMR, MultiVariableMaxAMR)
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
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
using namespace Sloth2D;

/// @brief Sharp interface x = x0 (or y = x0 if along_y), thickness 0.02
std::function<double(const mfem::Vector&, double)> interface(double x0, bool along_y = false) {
  return [x0, along_y](const mfem::Vector& p, double) {
    return 0.5 + 0.5 * std::tanh(((along_y ? p[1] : p[0]) - x0) / 0.02);
  };
}

/**
 * @brief Unit square, 8 x 8 quadrangles (non-conforming refinement allowed),
 * Kelly estimator. Variables are created by the test cases.
 */
struct AMRSetup {
  static constexpr double kCoarseArea = 1. / 64.;
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(8, 8, 1., 1.), true};
  BCS bcs{&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  mfem::ConstantCoefficient one{1.};
  mfem::DiffusionIntegrator integ{one};
  SlothErrorEstimators estimator{ErrorEstimatorType::KELLY, &integ};
  Parameters params{Parameter("max_elem_error", 1.e-4), Parameter("amr_max_level", 3),
                    Parameter("nc_limit", 1),           Parameter("max_preref_cycles", 3),
                    Parameter("scale_down_factor", 0.25), Parameter("rebalance", true)};

  mfem::ParMesh& mesh() { return *spatial.get_mesh(); }
  long long global_ne() { return mesh().GetGlobalNE(); }

  VAR variable(const std::string& name, std::function<double(const mfem::Vector&, double)> f) {
    return VAR(&spatial, bcs, name, Glossary::PhaseField, 2, AnalyticalFunctions<DIM>(f));
  }

  /// @brief Smallest element area over all the ranks
  double min_area() {
    double local = 1., global = 1.;
    for (int e = 0; e < mesh().GetNE(); ++e) local = std::min(local, mesh().GetElementVolume(e));
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global;
  }

  /// @brief Largest area of the elements whose centre satisfies pred (all ranks)
  double max_area_where(const std::function<bool(const mfem::Vector&)>& pred) {
    double local = 0., global = 0.;
    mfem::Vector c;
    for (int e = 0; e < mesh().GetNE(); ++e) {
      mesh().GetElementCenter(e, c);
      if (pred(c)) local = std::max(local, mesh().GetElementVolume(e));
    }
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
  }

  /// @brief L2 error of a variable against a function (all ranks)
  static double l2_error(VAR& v, std::function<double(const mfem::Vector&, double)> f) {
    mfem::FunctionCoefficient exact([f](const mfem::Vector& p) { return f(p, 0.); });
    return v.get_ref_gf().ComputeL2Error(exact);
  }
};

double x_of(const mfem::Vector& p, double) { return p[0]; }
}  // namespace

TEST_CASE_METHOD(AMRSetup, "SingleVariableAMR: initial refinement at the interface",
                 "[AMR]") {
  auto vars = VARS(variable("phi", interface(0.5)));
  auto aux = VARS(VAR(&spatial, bcs, "c", Glossary::Concentration, 2,
                      AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
                          x_of))));
  SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
  amr.SetCriteria(&estimator, params);

  REQUIRE(global_ne() == 64);
  const double coarse_error = l2_error(vars[0], interface(0.5));

  amr.InitialRefine(vars, {&aux});

  SECTION("the mesh is refined, up to amr_max_level, at the interface only") {
    CHECK(global_ne() > 64);
    CHECK(min_area() < kCoarseArea);
    CHECK(min_area() >= kCoarseArea / std::pow(4., 3) * (1. - 1.e-12));  // 3 levels at most
    CHECK_THAT(max_area_where([](const mfem::Vector& c) { return std::abs(c[0] - 0.5) > 0.3; }),
               WithinAbs(kCoarseArea, 1.e-12));  // far from the interface: not refined
  }
  SECTION("the initial condition is projected on the refined mesh: smaller error") {
    CHECK(l2_error(vars[0], interface(0.5)) < 0.5 * coarse_error);
  }
  SECTION("variables and auxiliary variables follow the mesh (sizes, values)") {
    CHECK(vars[0].get_unknown().Size() == vars[0].get_fespace()->GetTrueVSize());
    CHECK(aux[0].get_unknown().Size() == aux[0].get_fespace()->GetTrueVSize());
    // c = x is exact in Q1: its transfer to the refined (and rebalanced) mesh is exact
    CHECK(l2_error(aux[0], x_of) < 1.e-12);
  }
  SECTION("load balance after rebalancing") {
    const long long local = mesh().GetNE();
    long long max_ne = 0, min_ne = 0;
    MPI_Allreduce(&local, &max_ne, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local, &min_ne, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    CHECK(static_cast<double>(max_ne) <= 1.1 * static_cast<double>(min_ne) + 4.);
  }
}

TEST_CASE_METHOD(AMRSetup, "SingleVariableAMR: derefinement of a flat field", "[AMR]") {
  auto vars = VARS(variable("phi", interface(0.5)));
  SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
  amr.SetCriteria(&estimator, params);
  amr.InitialRefine(vars, {});
  const long long refined_ne = global_ne();
  REQUIRE(refined_ne > 64);

  // Uniform field: zero error everywhere, every refined element can be derefined
  mfem::Vector u = vars[0].get_unknown();
  u = 0.3;
  vars[0].update(u);
  amr.StepDerefine(vars, {});
  amr.StepRebalance(vars, {});

  CHECK(global_ne() < refined_ne);
  CHECK(vars[0].get_unknown().Size() == vars[0].get_fespace()->GetTrueVSize());
  CHECK(vars[0].get_unknown().Normlinf() <= 0.3 + 1.e-12);  // still 0.3 (if any dof)
}

TEST_CASE_METHOD(AMRSetup, "SingleVariableAMR: step refinement", "[AMR]") {
  auto vars = VARS(variable("phi", interface(0.5)));
  SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
  amr.SetCriteria(&estimator, params);
  amr.StepRefine(vars, {});  // one level, no re-projection of the initial condition
  amr.StepRebalance(vars, {});
  CHECK(global_ne() > 64);
  CHECK(vars[0].get_unknown().Size() == vars[0].get_fespace()->GetTrueVSize());
}

TEST_CASE_METHOD(AMRSetup, "MultiVariableMaxAMR: refinement driven by several variables",
                 "[AMR]") {
  // phi: interface x = 0.5, psi: interface y = 0.5
  auto vars = VARS(variable("phi", interface(0.5)), variable("psi", interface(0.5, true)));
  MultiVariableMaxAMR<VARS> amr(mesh(), spatial.is_nc_simplices(),
                                std::vector<std::string>{"phi", "psi"});
  amr.SetCriteria(&estimator, params);
  amr.InitialRefine(vars, {});

  auto near = [](double cx, double cy) {
    return [cx, cy](const mfem::Vector& c) {
      return std::abs(c[0] - cx) < 0.07 && std::abs(c[1] - cy) < 0.07;
    };
  };
  CHECK(max_area_where(near(0.5, 0.15)) < kCoarseArea);   // refined for phi only
  CHECK(max_area_where(near(0.15, 0.5)) < kCoarseArea);   // refined for psi only
  CHECK_THAT(max_area_where(near(0.1, 0.1)), WithinAbs(kCoarseArea, 1.e-12));  // neither
}

TEST_CASE_METHOD(AMRSetup, "AMR: invalid use", "[AMR][errors]") {
  auto vars = VARS(variable("phi", interface(0.5)));
  SECTION("criteria not set") {
    SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
    SLOTH_CHECK_MFEM_ERROR_WITH(amr.StepRefine(vars, {}), "AMR criteria not set");
  }
  SECTION("unknown pilot variable") {
    SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "psi");
    amr.SetCriteria(&estimator, params);
    SLOTH_CHECK_MFEM_ERROR_WITH(amr.StepRefine(vars, {}), "pilot variable 'psi' not found");
  }
  SECTION("max_elem_error is required") {
    SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
    SLOTH_CHECK_MFEM_ERROR_WITH(amr.SetCriteria(&estimator, Parameters(std::vector<Parameter>{})),
                                "Not Found");
  }
}

// -----------------------------------------------------------------------------
// AMR attached to a problem (set_amr): Allen-Cahn front at x = 0.5, Dirichlet
// phi = 0 on x = 0 and phi = 1 on x = 1. The mesh is refined before the first
// time step and at each step; the operator and the Dirichlet dofs follow the
// new mesh, and the Dirichlet values are kept exactly.
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(AMRSetup, "AMR: attached to a problem, during the time loop",
                 "[AMR][Problems]") {
  BCS bcs_d{&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 1.),
            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
  auto vars = VARS(VAR(&spatial, bcs_d, "phi", Glossary::PhaseField, 2,
                       AnalyticalFunctions<DIM>(interface(0.5))));
  const double lambda = 0.15, omega = 120.;
  Coefficients coefs(Coefficient(Glossary::FreeEnergy, Scheme::Explicit, W(omega)),
                     Coefficient(Glossary::Capillary, lambda), Coefficient(Glossary::Mobility, 1.),
                     Coefficient(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda)));
  std::vector<SPA*> spatials{&spatial};
  TransientOPE oper(spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", std::string("Newton")),
                                     Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12)));
  TransientPB problem("phi", oper, vars, {coefs});

  SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
  amr.SetCriteria(&estimator, params);
  problem.set_amr(&amr);

  auto coupling = Coupling("AMR", problem);
  auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                            Parameter("final_time", 2.e-3),
                                            Parameter("time_step", 1.e-3)),
                                 coupling);
  time.solve();

  CHECK(global_ne() > 64);
  mfem::ParFiniteElementSpace* fes = vars[0].get_fespace();
  const mfem::Vector u = vars[0].get_unknown();
  REQUIRE(u.Size() == fes->GetTrueVSize());

  // Dirichlet values on the refined mesh
  double local_err = 0.;
  for (const auto& [bdr, value] : {std::pair{1, 1.}, std::pair{3, 0.}}) {
    mfem::Array<int> marker(fes->GetParMesh()->bdr_attributes.Max());
    marker = 0;
    marker[bdr] = 1;
    mfem::Array<int> dofs;
    fes->GetEssentialTrueDofs(marker, dofs);
    for (int k = 0; k < dofs.Size(); ++k) {
      local_err = std::max(local_err, std::abs(u(dofs[k]) - value));
    }
  }
  double err = 0.;
  MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  CHECK(err < 1.e-10);
}

TEST_CASE_METHOD(AMRSetup, "MultiVariableMaxAMR: derefinement of flat fields", "[AMR]") {
  auto vars = VARS(variable("phi", interface(0.5)), variable("psi", interface(0.5, true)));
  MultiVariableMaxAMR<VARS> amr(mesh(), spatial.is_nc_simplices(),
                                std::vector<std::string>{"phi", "psi"});
  amr.SetCriteria(&estimator, params);
  amr.InitialRefine(vars, {});
  const long long refined_ne = global_ne();
  REQUIRE(refined_ne > 64);

  // Uniform fields: zero error everywhere, every refined element can be derefined
  for (std::size_t k = 0; k < 2; ++k) {
    mfem::Vector u = vars[k].get_unknown();
    u = 0.3;
    vars[k].update(u);
  }
  amr.StepDerefine(vars, {});
  CHECK(global_ne() < refined_ne);
}
