/**
 * @file test_MPIProblem.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for MPI_Problem (PB_MPI) and the default behaviour of ProblemBase
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
#include <tuple>

#include <catch2/catch_test_macros.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

// MPI_Problem fills its variable with the MPI rank of each process. It has no
// operator: initialization, CSV output and finalization are those of
// ProblemBase.
TEST_CASE("MPI_Problem: each rank writes its rank", "[Problems]") {
  using namespace Sloth2D;
  SPA spatial("InlineSquareWithQuadrangles", 1, 0, std::make_tuple(4, 4, 1., 1.));
  BCS bcs(&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann"));
  auto vars = VARS(VAR(&spatial, bcs, "rank", Glossary::Concentration, 2, -1.));

  sloth_unit::TemporaryDirectory dir("mpi");
  PST pst(&spatial, Parameters(Parameter("main_folder_path", dir.path()),
                               Parameter("calculation_path", std::string("MPI")),
                               Parameter("frequency", 1),
                               Parameter("enable_compute_energies", false)));
  PB_MPI problem(vars, pst);
  auto coupling = Coupling("MPI", problem);
  auto time = TimeDiscretization(Parameters(Parameter("initial_time", 0.),
                                            Parameter("final_time", 1.),
                                            Parameter("time_step", 1.)),
                                 coupling);
  time.solve();

  const mfem::Vector u = vars[0].get_unknown();
  const double rank = static_cast<double>(mfem::Mpi::WorldRank());
  for (int i = 0; i < u.Size(); ++i) {
    CAPTURE(i);
    CHECK(u(i) == rank);
  }
}
