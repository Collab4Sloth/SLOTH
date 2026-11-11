/**
 * @file test_CustomQuantity.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the quantities defined outside the kernel (CustomQuantity), read by
 *        custom integrators
 * @version 0.1
 * @date 2026-11-09
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
#include <string>
#include <tuple>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

namespace {
using namespace Sloth2D;

/**
 * @brief Custom diffusion integrator whose diffusion coefficient is the custom quantity
 *        "MyConductivity".
 */
template <class VARS>
class NamedConductivityIntegrator : public DiffusionNLFormIntegrator<VARS> {
 public:
  NamedConductivityIntegrator(Geometry geometry, const double time_step,
                              const std::vector<mfem::ParGridFunction>& u_old,
                              const std::vector<mfem::ParGridFunction>& aux_old,
                              const Parameters& params, std::vector<VARS*> auxvars,
                              const std::vector<Coefficients>& coefficients)
      : DiffusionNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                        coefficients) {
    this->integrator_name_ = "NamedConductivity";
    this->expected_list_.push_back(GlossaryType::Custom);
  }

 protected:
  void get_coefficients() override {
    this->check_custom_coefficient_names({"MyConductivity"});
    for (unsigned int i = 0; i < this->nb_blk_; i++) {
      this->diffusion.add(*(this->get_coefficient(i, "MyConductivity", 0)));
    }
  }
};

/**
 * @brief Steady heat conduction on the unit square (4 x 4 quadrangles, Q1), T = 1 on x = 0 and
 * T = 2 on x = 1: exact solution T = 1 + x (exact in Q1).
 * Boundaries: 0 lower (y = 0), 1 right (x = 1), 2 upper (y = 1), 3 left (x = 0).
 */
struct HeatSquare {
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(4, 4, 1., 1.)};
  std::vector<SPA*> spatials{&spatial};
  BCS bcs{&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 2.),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.)};
  VARS vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 1.5));
  Parameters newton{Parameter("description", std::string("Newton")), Parameter("rel_tol", 1.e-12),
                    Parameter("abs_tol", 1.e-12)};
};
}  // namespace

TEST_CASE("CustomQuantity: quantity defined outside the kernel", "[Glossary][Extensions]") {
  const GlossaryQuantity qty = CustomQuantity("MyConductivity", "W.m-1.K-1", "My conductivity", 2);
  CHECK(qty.type == GlossaryType::Custom);
  CHECK(qty.name == "MyConductivity");
  CHECK(qty.id == 2);
  CHECK(qty.description == "My conductivity in W.m-1.K-1");

  const Coefficient coef(CustomQuantity("MyConductivity"), 1.);
  CHECK(coef.get_type() == GlossaryType::Custom);
  CHECK(coef.get_custom_name() == "MyConductivity");
  CHECK(Coefficient(Glossary::Conductivity, 1.).get_custom_name().empty());
}

// A custom integrator reads its coefficient by name: T = 1 + x
TEST_CASE_METHOD(HeatSquare, "CustomQuantity: read by a custom integrator",
                 "[Operators][Extensions]") {
  SteadyOPE oper(spatials, {CustomIntegrator<NamedConductivityIntegrator>("NamedFourier")});
  oper.overload_nl_solver(NLSolverType::NEWTON, newton);
  SteadyPB problem("Heat", oper, vars,
                   {Coefficients(Coefficient(CustomQuantity("MyConductivity"), 1.))});
  auto coupling = Coupling("Heat", problem);
  auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                            Parameter("final_time", 1.),
                                            Parameter("time_step", 1.)),
                                 coupling);
  time.solve();

  mfem::FunctionCoefficient exact([](const mfem::Vector& p) { return 1. + p[0]; });
  CHECK(vars[0].get_ref_gf().ComputeL2Error(exact) < 1.e-10);
}

// The custom quantity expected by the custom integrator is missing
TEST_CASE_METHOD(HeatSquare, "CustomQuantity: missing custom quantity",
                 "[Operators][Extensions]") {
  const IntegratorEntry<VARS> named = CustomIntegrator<NamedConductivityIntegrator>("NamedFourier");
  SteadyOPE oper(spatials, {named});
  oper.set_coefficients({Coefficients(Coefficient(CustomQuantity("Other"), 1.))}, false);
  const std::vector<mfem::Vector> u{vars[0].get_unknown()};
  auto build = [&]() { return oper.set_nlfi_ptr(1., named, u); };
  SLOTH_CHECK_MFEM_ERROR_WITH(build(), "MyConductivity");
}

// One coefficient set shared by an operator with built-in integrators only and an operator
// mixing built-in and custom integrators: the custom quantity is ignored by the first one and
// read by the second one. T = 1 + x in both cases.
TEST_CASE_METHOD(HeatSquare, "CustomQuantity: coefficient set shared by built-in and custom",
                 "[Operators][Extensions]") {
  const Coefficients shared(Coefficient(Glossary::Conductivity, 1.),
                            Coefficient(CustomQuantity("MyConductivity"), 1.));
  mfem::FunctionCoefficient exact([](const mfem::Vector& p) { return 1. + p[0]; });

  SECTION("built-in integrators only") {
    SteadyOPE oper(spatials, {"Fourier"});
    oper.overload_nl_solver(NLSolverType::NEWTON, newton);
    SteadyPB problem("Heat", oper, vars, {shared});
    auto coupling = Coupling("Heat", problem);
    auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                              Parameter("final_time", 1.),
                                              Parameter("time_step", 1.)),
                                   coupling);
    time.solve();
    CHECK(vars[0].get_ref_gf().ComputeL2Error(exact) < 1.e-10);
  }

  SECTION("built-in and custom integrators") {
    SteadyOPE oper(spatials,
                   {"Fourier", CustomIntegrator<NamedConductivityIntegrator>("NamedFourier")});
    oper.overload_nl_solver(NLSolverType::NEWTON, newton);
    SteadyPB problem("Heat", oper, vars, {shared});
    auto coupling = Coupling("Heat", problem);
    auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                              Parameter("final_time", 1.),
                                              Parameter("time_step", 1.)),
                                   coupling);
    time.solve();
    CHECK(vars[0].get_ref_gf().ComputeL2Error(exact) < 1.e-10);
  }
}
