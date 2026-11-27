/**
 * @file test_Variable.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for Variable: initial conditions, analytical solutions, dependencies
 * @version 0.1
 * @date 2026-11-26
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
#include <set>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
using namespace Sloth2D;

/// @brief f(v_0, v_1, ...) = sum_k (k + 1) v_k, or 7 without value
class WeightedSum final : public FunctionCoefficient {
 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&, const unsigned int) {
      if (u.empty()) {
        return 7.;
      }
      double sum = 0.;
      for (std::size_t k = 0; k < u.size(); ++k) sum += static_cast<double>(k + 1) * u[k];
      return sum;
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

/// @brief f(x, y) = x + y, given through the coordinates
class CoordinateSum final : public FunctionCoefficient {
 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [this](const std::span<const double>&, const std::span<const double>&,
                  const std::span<const double>&,
                  const unsigned int) { return this->coordinates_[0] + this->coordinates_[1]; };
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

 public:
  bool uses_coordinates() const override { return true; }
};

double x_of(const mfem::Vector& p, double) { return p[0]; }

/**
 * @brief Unit square, 8 x 8 quadrangles, Q1. Variables are created by the test cases.
 */
struct VariableSetup {
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(8, 8, 1., 1.), true};
  BCS bcs{&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};

  mfem::ParMesh& mesh() { return *spatial.get_mesh(); }

  VAR x_variable(const std::string& name) {
    return VAR(&spatial, bcs, name, Glossary::Concentration, 2,
               AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(x_of)));
  }

  /// @brief L2 error of a variable against a function (all ranks)
  static double l2_error(VAR& v, std::function<double(const mfem::Vector&)> f) {
    mfem::FunctionCoefficient exact(f);
    return v.get_ref_gf().ComputeL2Error(exact);
  }
};
}  // namespace

// -----------------------------------------------------------------------------
// The historical ways of giving an initial condition and an analytical solution
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(VariableSetup, "Variable: initial condition and analytical solution",
                 "[Variables]") {
  std::function<double(const mfem::Vector&, double)> f = x_of;
  const mfem::FunctionCoefficient fc(f);
  const auto af = AnalyticalFunctions<DIM>(f);

  SECTION("double, AnalyticalFunctions, mfem::FunctionCoefficient, SLOTH Coefficient") {
    auto v_double = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, 1.5);
    auto v_int = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, 0);
    auto v_af = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, af);
    auto v_fc = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, fc);
    auto v_coef =
        VAR(&spatial, bcs, "v", Glossary::Concentration, 2, Coefficient(Glossary::Default, 2.5));

    CHECK(l2_error(v_double, [](const mfem::Vector&) { return 1.5; }) < 1.e-12);
    CHECK(l2_error(v_int, [](const mfem::Vector&) { return 0.; }) < 1.e-12);
    CHECK(l2_error(v_af, [](const mfem::Vector& p) { return p[0]; }) < 1.e-12);
    CHECK(l2_error(v_fc, [](const mfem::Vector& p) { return p[0]; }) < 1.e-12);
    CHECK(l2_error(v_coef, [](const mfem::Vector&) { return 2.5; }) < 1.e-12);
    CHECK(v_af.isInitialConditionReprojectable());
    CHECK(v_af.get_dependencies().empty());
  }
  SECTION("the history is initialized to the initial condition") {
    auto v = VAR(&spatial, bcs, "v", Glossary::Concentration, 3, af);
    for (const auto& [key, vec] : v.get_map_unknown()) {
      mfem::Vector diff(vec);
      diff -= v.get_unknown();
      CHECK(diff.Normlinf() < 1.e-14);
    }
  }
  SECTION("analytical solution given as AnalyticalFunctions") {
    auto v = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, 1.5, af);
    CHECK(v.get_analytical_solution() != nullptr);
    REQUIRE(v.get_analytical_solution_coefficient() != nullptr);
    CHECK(v.get_ref_gf().ComputeL2Error(*v.get_analytical_solution_coefficient()) ==
          v.get_ref_gf().ComputeL2Error(*std::make_shared<mfem::FunctionCoefficient>(f)));
  }
  SECTION("analytical solution given as mfem::FunctionCoefficient or SLOTH Coefficient") {
    auto v_fc = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, af, fc);
    auto v_coef = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, af, std::set<std::string>{},
                      Coefficient(Glossary::Default, Scheme::Implicit, CoordinateSum()));
    // Only an AnalyticalFunctions can be returned as a std::function
    CHECK(v_fc.get_analytical_solution() == nullptr);
    REQUIRE(v_fc.get_analytical_solution_coefficient() != nullptr);
    CHECK(v_fc.get_ref_gf().ComputeL2Error(*v_fc.get_analytical_solution_coefficient()) < 1.e-12);
    REQUIRE(v_coef.get_analytical_solution_coefficient() != nullptr);
    mfem::FunctionCoefficient y_of([](const mfem::Vector& p) { return p[1]; });
    mfem::ParGridFunction diff(v_coef.get_fespace());
    diff.ProjectCoefficient(*v_coef.get_analytical_solution_coefficient());
    diff -= v_coef.get_ref_gf();
    CHECK(diff.ComputeL2Error(y_of) < 1.e-12);
  }
}

// -----------------------------------------------------------------------------
// Initial condition given by a SLOTH Coefficient depending on other variables
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(VariableSetup, "Variable: initial condition depending on other variables",
                 "[Variables]") {
  const Coefficient weighted_sum(Glossary::Default, Scheme::Implicit, WeightedSum());

  SECTION("set_dependencies(): values of the dependencies, in the order of the arguments") {
    auto c = VAR(&spatial, bcs, "c", Glossary::Concentration, 2, 0.5);
    auto x = x_variable("x");
    auto v = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, weighted_sum);
    v.set_dependencies(c, x);

    CHECK(v.get_dependencies() == std::vector<std::string>{"c", "x"});
    CHECK(l2_error(v, [](const mfem::Vector& p) { return 0.5 + 2. * p[0]; }) < 1.e-12);
    mfem::Vector diff(v.get_last());
    diff -= v.get_unknown();
    CHECK(diff.Normlinf() < 1.e-14);
  }
  SECTION("without set_dependencies(): projected at the first use") {
    auto v = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, weighted_sum);
    auto w = VAR(&spatial, bcs, "w", Glossary::Concentration, 2,
                 Coefficient(Glossary::Default, Scheme::Implicit, CoordinateSum()));
    CHECK(v.get_unknown().Size() == v.get_fespace()->GetTrueVSize());
    CHECK(l2_error(v, [](const mfem::Vector&) { return 7.; }) < 1.e-12);
    CHECK(l2_error(w, [](const mfem::Vector& p) { return p[0] + p[1]; }) < 1.e-12);
  }
  SECTION("the dependencies are kept by the copies stored in a Variables container") {
    auto c = VAR(&spatial, bcs, "c", Glossary::Concentration, 2, 0.5);
    auto v = VAR(&spatial, bcs, "v", Glossary::Concentration, 2, weighted_sum);
    v.set_dependencies(c);
    auto vars = VARS(v, c);
    CHECK(vars[0].get_dependencies() == std::vector<std::string>{"c"});
    CHECK(l2_error(vars[0], [](const mfem::Vector&) { return 0.5; }) < 1.e-12);
  }
}

// -----------------------------------------------------------------------------
// Reprojection of the initial conditions during the initial refinement
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(VariableSetup, "Variable: initial refinement with dependent initial conditions",
                 "[Variables][AMR]") {
  const auto interface = [](const mfem::Vector& p, double) {
    return 0.5 + 0.5 * std::tanh((p[0] - 0.5) / 0.02);
  };
  auto vars = VARS(
      VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2,
          AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(interface))));

  // v depends on c and x, and is stored first: the dependencies must be reset before it
  auto c = VAR(&spatial, bcs, "c", Glossary::Concentration, 2, 0.5);
  auto x = x_variable("x");
  auto v = VAR(&spatial, bcs, "v", Glossary::Concentration, 2,
               Coefficient(Glossary::Default, Scheme::Implicit, WeightedSum()));
  v.set_dependencies(c, x);
  auto aux = VARS(v, c, x);

  mfem::ConstantCoefficient one{1.};
  mfem::DiffusionIntegrator integ{one};
  SlothErrorEstimators estimator{ErrorEstimatorType::KELLY, &integ};
  Parameters params{Parameter("max_elem_error", 1.e-4),
                    Parameter("amr_max_level", 2),
                    Parameter("nc_limit", 1),
                    Parameter("max_preref_cycles", 2),
                    Parameter("scale_down_factor", 0.25),
                    Parameter("rebalance", true)};
  SingleVariableAMR<VARS> amr(mesh(), spatial.is_nc_simplices(), "phi");
  amr.SetCriteria(&estimator, params);

  REQUIRE(mesh().GetGlobalNE() == 64);
  amr.InitialRefine(vars, {&aux});
  REQUIRE(mesh().GetGlobalNE() > 64);

  CHECK(aux[0].get_unknown().Size() == aux[0].get_fespace()->GetTrueVSize());
  CHECK(l2_error(aux[0], [](const mfem::Vector& p) { return 0.5 + 2. * p[0]; }) < 1.e-12);
  mfem::Vector diff(aux[0].get_last());
  diff -= aux[0].get_unknown();
  CHECK(diff.Normlinf() < 1.e-14);
}
