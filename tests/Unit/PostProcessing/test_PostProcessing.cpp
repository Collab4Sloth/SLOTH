/**
 * @file test_PostProcessing.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for PostProcessing (output schedule, ParaView files, CSV of errors,
 *        integrals and iso-values)
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
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
using namespace Sloth2D;

/// @brief Splits a CSV line
std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

/// @brief First (header) and last lines of a CSV file
std::pair<std::vector<std::string>, std::vector<std::string>> header_and_last(
    const std::filesystem::path& file) {
  std::ifstream in(file);
  std::string line, header, last;
  std::getline(in, header);
  while (std::getline(in, line)) {
    if (!line.empty()) last = line;
  }
  return {split(header), split(last)};
}

/// @brief Value of the column named `name` in the last line of a CSV file
double column(const std::filesystem::path& file, const std::string& name) {
  const auto [header, last] = header_and_last(file);
  for (std::size_t k = 0; k < header.size() && k < last.size(); ++k) {
    if (header[k] == name) return std::stod(last[k]);
  }
  FAIL("column " << name << " not found in " << file);
  return 0.;
}

/// @brief True if a file with the given extension exists below dir
bool has_file_with_extension(const std::filesystem::path& dir, const std::string& ext) {
  for (const auto& entry : std::filesystem::recursive_directory_iterator(dir)) {
    if (entry.path().extension() == ext) return true;
  }
  return false;
}

struct PstSetup {
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(8, 8, 1., 1.)};
  sloth_unit::TemporaryDirectory dir{"pst"};
  PstSetup() { Profiling::getInstance().disable(); }

  Parameters base_parameters() const {
    return Parameters(Parameter("main_folder_path", dir.path()),
                      Parameter("calculation_path", std::string("Calc")));
  }
};
}  // namespace

// -----------------------------------------------------------------------------
// Output schedule
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(PstSetup, "PostProcessing: when are the fields saved", "[PostProcessing]") {
  SECTION("every 2 iterations") {
    const Parameters p = base_parameters() + Parameter("frequency", 2);
    PST pst(&spatial, p);
    CHECK(pst.get_frequency() == 2);
    CHECK(pst.need_to_be_saved(0, 0.));
    CHECK_FALSE(pst.need_to_be_saved(1, 0.1));
    CHECK(pst.need_to_be_saved(2, 0.2));
  }
  SECTION("given iterations") {
    const Parameters p = base_parameters() + Parameter("iterations_list", vInt{1, 3});
    PST pst(&spatial, p);
    CHECK(pst.need_to_be_saved(1, 0.1));
    CHECK_FALSE(pst.need_to_be_saved(2, 0.2));
    CHECK(pst.need_to_be_saved(3, 0.3));
  }
  SECTION("given times") {
    const Parameters p = base_parameters() + Parameter("times_list", vDouble{0.5});
    PST pst(&spatial, p);
    CHECK(pst.need_to_be_saved(7, 0.5));
    CHECK_FALSE(pst.need_to_be_saved(8, 0.6));
  }
  SECTION("specialized outputs (CSV) every 3 iterations") {
    const Parameters p =
        base_parameters() + Parameter("frequency", 1) + Parameter("specialized_frequency", 3);
    PST pst(&spatial, p);
    CHECK(pst.get_specialized_frequency() == 3);
    CHECK(pst.need_specialized(3));
    CHECK_FALSE(pst.need_specialized(4));
  }
}

TEST_CASE_METHOD(PstSetup, "PostProcessing: invalid parameters", "[PostProcessing][errors]") {
  SECTION("no output schedule") {
    const Parameters p = base_parameters();
    SLOTH_CHECK_MFEM_ERROR_WITH(PST(&spatial, p), "At least one");
  }
  SECTION("specialized_frequency < 1") {
    const Parameters p =
        base_parameters() + Parameter("frequency", 1) + Parameter("specialized_frequency", 0);
    SLOTH_CHECK_MFEM_ERROR_WITH(PST(&spatial, p), "specialized_frequency");
  }
  SECTION("integral bounds in the wrong order") {
    const Parameters p = base_parameters() + Parameter("frequency", 1) +
                         Parameter("integral_to_compute", MapString2Double{{"T", {1., 0.}}});
    SLOTH_CHECK_MFEM_ERROR_WITH(PST(&spatial, p), "Lower bound");
  }
}

TEST_CASE_METHOD(PstSetup, "PostProcessing: force_clean_output_dir", "[PostProcessing]") {
  const std::filesystem::path old_file = std::filesystem::path(dir.path()) / "old.txt";
  if (mfem::Mpi::Root()) std::ofstream(old_file) << "old";
  MPI_Barrier(MPI_COMM_WORLD);
  const Parameters p = base_parameters() + Parameter("frequency", 1) +
                       Parameter("force_clean_output_dir", true);
  PST pst(&spatial, p);
  MPI_Barrier(MPI_COMM_WORLD);
  CHECK_FALSE(std::filesystem::exists(old_file));
}

// -----------------------------------------------------------------------------
// Complete problem with post-processing: steady heat conduction, T = 1 + x
// (T = 1 on x = 0, T = 2 on x = 1), exact in Q1.
//   - L2 error against the analytical solution: 0
//   - integral and average of T over the unit square: 3/2
//   - iso-value T = 3/2: the line x = 1/2
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(PstSetup, "PostProcessing: complete problem (ParaView, errors, integral, iso)",
                 "[PostProcessing][Problems]") {
  auto bcs = BCS(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 2.),
                 Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.));
  auto exact = AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& p, double) { return 1. + p[0]; }));
  auto vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2,
                       AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Uniform, 1.5), exact));

  const Parameters p_pst = base_parameters() + Parameter("frequency", 1) +
                           Parameter("enable_compute_energies", false) +
                           Parameter("integral_to_compute", MapString2Double{{"T", {0., 10.}}}) +
                           Parameter("iso_val_to_compute", MapStringDouble{{"T", 1.5}});
  PST pst(&spatial, p_pst);

  Coefficients coefs(Coefficient(Glossary::Conductivity, 3.));
  std::vector<SPA*> spatials{&spatial};
  SteadyOPE oper(spatials, {"Fourier"});
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", std::string("Newton")),
                                     Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12)));
  SteadyPB problem("Heat", oper, vars, {coefs}, pst);
  auto coupling = Coupling("Heat", problem);
  auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                            Parameter("final_time", 0.1),
                                            Parameter("time_step", 0.1)),
                                 coupling);
  time.solve();
  MPI_Barrier(MPI_COMM_WORLD);

  // Files are written by rank 0
  if (mfem::Mpi::Root()) {
    const std::filesystem::path out = pst.get_post_processing_directory();
    const auto csv = out / "time_specialized.csv";
    REQUIRE(std::filesystem::exists(csv));
    CHECK(column(csv, "T_L2-error[-]") < 1.e-10);
    CHECK_THAT(column(csv, "T_integral[-]"), WithinAbs(1.5, 1.e-12));
    CHECK_THAT(column(csv, "T_average[-]"), WithinAbs(1.5, 1.e-12));

    // Iso-value T = 1.5: every point found has x = 0.5
    const auto iso = out / "T_iso_computation.csv";
    REQUIRE(std::filesystem::exists(iso));
    const auto [header, last] = header_and_last(iso);
    REQUIRE(last.size() > 3);
    for (std::size_t k = 3; k < last.size(); k += 2) {  // after iter, dt, t: (x, y) pairs
      CAPTURE(k);
      CHECK_THAT(std::stod(last[k]), WithinAbs(0.5, 1.e-12));
    }

    CHECK(has_file_with_extension(dir.path(), ".pvd"));  // ParaView collection
  }
}
