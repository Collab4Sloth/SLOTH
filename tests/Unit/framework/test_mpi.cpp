/**
 * @file test_mpi.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Sanity checks of the unit-test harness itself (MPI initialisation)
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

#include "mfem.hpp"

TEST_CASE("Harness: MPI is initialised by the common main", "[Framework][MPI]") {
  REQUIRE(mfem::Mpi::IsInitialized());

  const int size = mfem::Mpi::WorldSize();
  const int rank = mfem::Mpi::WorldRank();
  CHECK(size >= 1);
  CHECK(rank >= 0);
  CHECK(rank < size);

  // Every rank contributes 1: the sum must be the number of ranks
  int one = 1;
  int sum = 0;
  MPI_Allreduce(&one, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  CHECK(sum == size);
}
