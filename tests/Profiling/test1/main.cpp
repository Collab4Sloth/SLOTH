/**
 * @file ex1.cpp
 * @author MB278596 (mouad.bakhkhakh@cea.fr)
 * @brief This example code compute the sum of numbers from 1 to n using MPI
 * @version 0.1
 * @date 2024-06-10
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <mpi.h>

#include <chrono>  // NOLINT [build/c++11]
#include <iomanip>
#include <iostream>

#include "MAToolsProfiling/MATimersAPI.hxx"

int main(int argc, char** argv) {
  //---------------------------------------
  // Initialize MPI
  //---------------------------------------
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  {
    Catch_Time_Section("test1");
    //---------------------------------------

    const int& n = 1e6;
    int somme_locale = 0;
    int somme_globale = 0;

    // Compute the local sum for each proc
    for (int i = rank + 1; i <= n; i += size) {
      Catch_Time_Section("Local sum");
      somme_locale += i;
    }

    // Reduction of the local sum
    MPI_Reduce(&somme_locale, &somme_globale, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      std::cout << "The sum of numbers from 1 to  " << n << " is : " << somme_globale << std::endl;
    }
  }

  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------

  return 0;
}
