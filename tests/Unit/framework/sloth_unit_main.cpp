/**
 * @file sloth_unit_main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Common main for all SLOTH unit tests (Catch2 v3 + MPI)
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
 * What this main does:
 *  1. initialises MPI and Hypre through MFEM, as the other SLOTH tests do
 *     (mfem::Mpi::Init also finalises MPI at exit);
 *  2. runs Catch2 on every rank, with the usual Catch2 command-line options;
 *  3. returns the same non-zero code on every rank if any rank failed, so that
 *     mpirun and ctest see a failure even when only rank 1 fails.
 */
#include <catch2/catch_session.hpp>

#include "mfem.hpp"

int main(int argc, char* argv[]) {
  // 1. MPI + Hypre
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();

  // 2. Catch2 session (same behaviour as Catch2::Catch2WithMain)
  Catch::Session session;
  const int cli_status = session.applyCommandLine(argc, argv);
  if (cli_status != 0) {
    return cli_status;  // invalid command line
  }
  // Number of failed test cases on this rank (0 if all passed)
  const int local_failures = session.run();

  // 3. Global result: the worst result over all ranks
  int global_failures = 0;
  MPI_Allreduce(&local_failures, &global_failures, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (mfem::Mpi::Root() && global_failures != 0 && local_failures == 0) {
    std::cerr << "[SLOTH unit tests] failure(s) on other MPI rank(s)" << std::endl;
  }
  return global_failures;
}
