/**
 * @file test_InterDiffusionCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for InterDiffusionCoefficient (and PropertyBase)
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
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "Parameters/Parameters.hpp"
#include "Property/InterDiffusionCoefficient.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;
using sloth_unit::point;

namespace {
// Outputs: [name, element, "inter_mob"] -> vector filled by the property
using Output =
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>;
// Inputs: [element, "x"], [phase, element, "mob"], [phase, "phi"] -> values at the nodes
using Input = std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>;

// Binary U-O system on 3 nodes
const mfem::Vector xU = point({0.2, 0.5, 0.8});
const mfem::Vector xO = point({0.8, 0.5, 0.2});
constexpr double MU = 1.e-3;  // mobility of U
constexpr double MO = 4.e-3;  // mobility of O

/// @brief Inter-diffusion mobility of a binary system: x_U x_O (x_O M_U + x_U M_O)
double binary_L(double x_u, double x_o, double m_u, double m_o) {
  return x_u * x_o * (x_o * m_u + x_u * m_o);
}

mfem::Vector constant(double v) { return point({v, v, v}); }

/// @brief Parameters of the U-O system (U first, O last). Must outlive the
/// property: PropertyBase keeps a reference to them.
Parameters uo_parameters() {
  return Parameters(Parameter("first_component", std::string("U")),
                    Parameter("last_component", std::string("O")),
                    Parameter("primary_phase", std::string("FLUORITE")),
                    Parameter("secondary_phase", std::string("LIQUID")));
}

/// @brief Single-phase inputs: [element, "x"] and [phase, element, "mob"]
Input single_phase_input() {
  return {{{"U", "x"}, xU},
          {{"O", "x"}, xO},
          {{"FLUORITE", "U", "mob"}, constant(MU)},
          {{"FLUORITE", "O", "mob"}, constant(MO)}};
}
}  // namespace

TEST_CASE("InterDiffusionCoefficient: binary system, one phase", "[Property]") {
  const Parameters params = uo_parameters();
  InterDiffusionCoefficient prop(params);

  mfem::Vector L_U, L_O;
  Output output{{{"L_U", "U", "inter_mob"}, std::ref(L_U)},
                {{"L_O", "O", "inter_mob"}, std::ref(L_O)}};
  prop.compute(output, single_phase_input());

  REQUIRE(L_U.Size() == 3);
  REQUIRE(L_O.Size() == 3);
  for (int k = 0; k < 3; ++k) {
    CAPTURE(k);
    const double L = binary_L(xU(k), xO(k), MU, MO);
    CHECK_THAT(L_U(k), WithinRel(L, 1.e-12));
    CHECK_THAT(L_O(k), WithinRel(-L, 1.e-12));  // opposite fluxes
  }

  SECTION("a second call uses the new inputs (no state kept between calls)") {
    Input input = single_phase_input();
    std::get<1>(input[2]) = constant(2. * MU);
    prop.compute(output, input);
    CHECK_THAT(L_U(1), WithinRel(binary_L(0.5, 0.5, 2. * MU, MO), 1.e-12));
  }
  SECTION("set_times") { CHECK_NOTHROW(prop.set_times(1., 0.1)); }
  SECTION("phase and element names are case-insensitive") {
    const Input input{{{"u", "x"}, xU},
                      {{"o", "x"}, xO},
                      {{"Fluorite", "u", "mob"}, constant(MU)},
                      {{"fluorite", "O", "mob"}, constant(MO)}};
    prop.compute(output, input);
    CHECK_THAT(L_U(0), WithinRel(binary_L(xU(0), xO(0), MU, MO), 1.e-12));
  }
}

TEST_CASE("InterDiffusionCoefficient: binary system, two phases", "[Property]") {
  const Parameters params = uo_parameters();
  InterDiffusionCoefficient prop(params);

  // phi = 1: fluorite only; phi = 0: liquid only; phi = 0.5: mixture.
  // phi is clamped to [0, 1] (1.5 -> 1).
  const mfem::Vector phi = point({1.5, 0., 0.5});
  constexpr double MU_liq = 5.e-3, MO_liq = 2.e-3;
  Input input = single_phase_input();
  input.push_back({{"LIQUID", "U", "mob"}, constant(MU_liq)});
  input.push_back({{"LIQUID", "O", "mob"}, constant(MO_liq)});
  input.push_back({{"FLUORITE", "phi"}, phi});

  mfem::Vector L_U, L_O;
  Output output{{{"L_U", "U", "inter_mob"}, std::ref(L_U)},
                {{"L_O", "O", "inter_mob"}, std::ref(L_O)}};
  prop.compute(output, input);

  for (int k = 0; k < 3; ++k) {
    CAPTURE(k);
    const double p = std::clamp(phi(k), 0., 1.);
    const double mu = p * MU + (1. - p) * MU_liq;
    const double mo = p * MO + (1. - p) * MO_liq;
    CHECK_THAT(L_U(k), WithinRel(binary_L(xU(k), xO(k), mu, mo), 1.e-12));
  }
}

TEST_CASE("InterDiffusionCoefficient: invalid data", "[Property][errors]") {
  const Parameters params = uo_parameters();
  mfem::Vector L_U;

  SECTION("missing parameter") {
    const Parameters incomplete(Parameter("first_component", std::string("U")));
    SLOTH_CHECK_MFEM_ERROR_WITH(InterDiffusionCoefficient{incomplete}, "Not Found");
  }
  SECTION("no inter_mob output") {
    InterDiffusionCoefficient prop(params);
    Output output{{{"L_U", "U", "other"}, std::ref(L_U)}};
    SLOTH_CHECK_MFEM_ERROR_WITH(prop.compute(output, single_phase_input()), "inter_mob");
  }
  SECTION("inter_mob output without element") {
    InterDiffusionCoefficient prop(params);
    Output output{{{"L_U", "inter_mob"}, std::ref(L_U)}};
    SLOTH_CHECK_MFEM_ERROR_WITH(prop.compute(output, single_phase_input()), "inter_mob");
  }
  SECTION("no molar fraction") {
    InterDiffusionCoefficient prop(params);
    Output output{{{"L_U", "U", "inter_mob"}, std::ref(L_U)}};
    Input input{{{"FLUORITE", "U", "mob"}, constant(MU)}};
    SLOTH_CHECK_MFEM_ERROR_WITH(prop.compute(output, input), "x auxiliary");
  }
  SECTION("no mobility") {
    InterDiffusionCoefficient prop(params);
    Output output{{{"L_U", "U", "inter_mob"}, std::ref(L_U)}};
    Input input{{{"U", "x"}, xU}};
    SLOTH_CHECK_MFEM_ERROR_WITH(prop.compute(output, input), "mob auxiliary");
  }
  SECTION("mobility without phase") {
    InterDiffusionCoefficient prop(params);
    Output output{{{"L_U", "U", "inter_mob"}, std::ref(L_U)}};
    Input input{{{"U", "x"}, xU}, {{"U", "mob"}, constant(MU)}};
    SLOTH_CHECK_MFEM_ERROR_WITH(prop.compute(output, input), "mobilities");
  }
  SECTION("inter-diffusion mobility >= 1 is rejected") {
    InterDiffusionCoefficient prop(params);
    mfem::Vector L_O;
    Output output{{{"L_U", "U", "inter_mob"}, std::ref(L_U)},
                  {{"L_O", "O", "inter_mob"}, std::ref(L_O)}};
    Input input = single_phase_input();
    std::get<1>(input[2]) = constant(100.);
    SLOTH_CHECK_MFEM_ERROR_WITH(prop.compute(output, input), "Error for j");
  }
}
