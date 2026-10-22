/**
 * @file test_Options.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the keyword -> enum conversions (from) of Options
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
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "Options/Options.hpp"
#include "Utils/UtilsForData.hpp"
#include "sloth_unit.hpp"

namespace {
/**
 * @brief Checks Enum::from on a table {keyword, expected value}:
 *  - every keyword gives the expected value,
 *  - whatever the case (keywords are case-insensitive),
 *  - an unknown keyword is an error whose message names the option.
 *
 * The table is written here from the documentation of the options, not copied
 * from the code: a keyword renamed or removed in SLOTH makes the test fail.
 */
template <typename Enum>
void check_from(const std::string& option_name,
                const std::vector<std::pair<std::string, typename Enum::value>>& table) {
  for (const auto& [keyword, expected] : table) {
    CAPTURE(keyword);  // printed with the failure, to know which keyword failed
    CHECK(Enum::from(keyword) == expected);
    CHECK(Enum::from(toUpperCase(keyword)) == expected);
    CHECK(Enum::from(toLowerCase(keyword)) == expected);
  }
  SLOTH_CHECK_MFEM_ERROR_WITH(Enum::from("NotAKeyword"), option_name);
}
}  // namespace

TEST_CASE("Options: TimeScheme::from", "[Options]") {
  check_from<TimeScheme>("TimeScheme", {{"EulerImplicit", TimeScheme::EulerImplicit},
                                        {"SDIRK23", TimeScheme::SDIRK23},
                                        {"SDIRK33", TimeScheme::SDIRK33},
                                        {"EulerExplicit", TimeScheme::EulerExplicit},
                                        {"RungeKutta4", TimeScheme::RungeKutta4}});
}

TEST_CASE("Options: Problems::from", "[Options]") {
  check_from<Problems>("Problems", {{"Diffusion", Problems::Diffusion},
                                    {"AllenCahn", Problems::AllenCahn},
                                    {"Calphad", Problems::Calphad}});
}

TEST_CASE("Options: Meshes::from", "[Options]") {
  check_from<Meshes>("Meshes",
                     {{"InlineLineWithSegments", Meshes::InlineLineWithSegments},
                      {"InlineSquareWithTriangles", Meshes::InlineSquareWithTriangles},
                      {"InlineSquareWithQuadrangles", Meshes::InlineSquareWithQuadrangles},
                      {"InlineSquareWithTetraedres", Meshes::InlineSquareWithTetraedres},
                      {"InlineSquareWithHexaedres", Meshes::InlineSquareWithHexaedres},
                      {"GMSH", Meshes::GMSH}});
}

TEST_CASE("Options: BoundaryConditionType::from", "[Options]") {
  check_from<BoundaryConditionType>("BoundaryConditionType",
                                    {{"Dirichlet", BoundaryConditionType::Dirichlet},
                                     {"Neumann", BoundaryConditionType::Neumann},
                                     {"Robin", BoundaryConditionType::Robin},
                                     {"Periodic", BoundaryConditionType::Periodic}});
}

TEST_CASE("Options: Integrators::from", "[Options]") {
  check_from<Integrators>("Integrators",
                          {{"TimeDerivative", Integrators::TimeDerivative},
                           {"HeatTimeDerivative", Integrators::HeatTimeDerivative},
                           {"SplitTimeDerivative", Integrators::SplitTimeDerivative},
                           {"Fick", Integrators::Fick},
                           {"Fourier", Integrators::Fourier},
                           {"CahnHilliard", Integrators::CahnHilliard},
                           {"AllenCahn", Integrators::AllenCahn},
                           {"SplitAllenCahn", Integrators::SplitAllenCahn},
                           {"MeltingTemperature", Integrators::MeltingTemperature},
                           {"MeltingCalphad", Integrators::MeltingCalphad},
                           {"MeltingConstant", Integrators::MeltingConstant},
                           {"MassFlux", Integrators::MassFlux},
                           {"LatentHeat", Integrators::LatentHeat},
                           {"Neumann", Integrators::Neumann},
                           {"Robin", Integrators::Robin},
                           {"SourceTerm", Integrators::SourceTerm},
                           {"Advection", Integrators::Advection},
                           {"SplitAdvection", Integrators::SplitAdvection},
                           {"AdvectionBoundary", Integrators::AdvectionBoundary},
                           {"AdvectionBoundaryCH", Integrators::AdvectionBoundaryCH}});
}

TEST_CASE("Options: Calphad options", "[Options][Calphad]") {
  SECTION("calphad_outputs::from") {
    check_from<calphad_outputs>(
        "calphad_outputs",
        {{"mu", calphad_outputs::mu},   {"dmu", calphad_outputs::dmu},
         {"x", calphad_outputs::x},     {"xp", calphad_outputs::xp},
         {"xph", calphad_outputs::xph}, {"y", calphad_outputs::y},
         {"g", calphad_outputs::g},     {"gm", calphad_outputs::gm},
         {"h", calphad_outputs::h},     {"hm", calphad_outputs::hm},
         {"dgm", calphad_outputs::dgm}, {"cp", calphad_outputs::cp},
         {"mob", calphad_outputs::mob}, {"nucleus", calphad_outputs::nucleus},
         {"error", calphad_outputs::error}});
  }
  SECTION("temperature_sort_method::from") {
    check_from<temperature_sort_method>("temperature_sort_method",
                                        {{"Ascending", temperature_sort_method::ascending},
                                         {"Descending", temperature_sort_method::descending},
                                         {"No", temperature_sort_method::no}});
  }
  SECTION("pressure_sort_method::from") {
    check_from<pressure_sort_method>("pressure_sort_method",
                                     {{"Ascending", pressure_sort_method::ascending},
                                      {"Descending", pressure_sort_method::descending},
                                      {"No", pressure_sort_method::no}});
  }
  SECTION("KKS_nucleation_strategy::from") {
    check_from<KKS_nucleation_strategy>(
        "KKS_nucleation_strategy",
        {{"LiquidFraction", KKS_nucleation_strategy::liquid_fraction},
         {"GivenMeltingTemperature", KKS_nucleation_strategy::given_melting_temperature}});
  }
}
