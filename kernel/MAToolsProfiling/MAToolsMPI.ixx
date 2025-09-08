/**
 * @file MAToolsMPI.ixx
 * @author Raphaël Prat (raphael.prat@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-09-08
 *
 * Copyright CEA (C) 2025
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
#pragma once
#include <cassert>
#include <iostream>
#include <mfem.hpp>

/**
 * @namespace MATools
 * @brief Namespace containing utility tools for various purposes.
 */
namespace MATools {
/**
 * @namespace MPI
 * @brief Namespace containing MPI-related utilities.
 */
namespace MPI {
constexpr int master = 0;

void mpi_initialize([[maybe_unused]] int* argc, [[maybe_unused]] char*** argv) {
  const auto is_init = check_mpi_initialized();
  if (!is_init) {
    MPI_Init(argc, argv);
  }
}

void mpi_finalize() {
  const auto is_final = check_mpi_finalized();
  if (!is_final) {
    MPI_Finalize();
  }
}

bool check_mpi_initialized() {
  int val = -1;
  MPI_Initialized(&val);
  assert(val != -1 && "error in check mpi init");
  bool ret = val == 1 ? true : false;
  return ret;
}

bool check_mpi_finalized() {
  int val = -1;
  MPI_Finalized(&val);
  assert(val != -1 && "error in check mpi finalize");
  bool ret = val == 1 ? true : false;
  return ret;
}

int get_rank() { return mfem::Mpi::WorldRank(); }

int get_mpi_size() { return mfem::Mpi::WorldSize(); }

bool is_master() { return mfem::Mpi::Root(); }

template <typename T>
T reduce(T a_in, MPI_Op a_op) {
  std::cout << "error" << std::endl;
  std::exit(EXIT_FAILURE);
  return -666;
}

template <>
double reduce(double a_in, MPI_Op a_op) {
  double global(0.0);
  MPI_Reduce(&a_in, &global, 1, MPI_DOUBLE, a_op, master, MPI_COMM_WORLD);  // master rank is 0
  return global;
}

template <>
int reduce(int a_in, MPI_Op a_op) {
  int global(0.0);
  MPI_Reduce(&a_in, &global, 1, MPI_INT, a_op, master, MPI_COMM_WORLD);  // master rank is 0
  return global;
}

double reduce_max(double a_duration) {
  double ret = reduce(a_duration, MPI_MAX);
  return ret;
}

double reduce_mean(double a_duration) {
  double ret = reduce(a_duration, MPI_SUM);
  auto mpi_size = get_mpi_size();
  assert(mpi_size > 0);
  return ret / mpi_size;
}
};  // namespace MPI
};  // namespace MATools
