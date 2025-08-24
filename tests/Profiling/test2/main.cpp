/**
 * @file main.cpp
 * @author MB278596 (mouad.bakhkhakh@cea.fr)
 * @brief This example prints a Hello message from each proc and computes the time elapsed using
 * MPI_Wtime and the profiling classes
 * @version 0.1
 * @date 2024-06-10
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <mpi.h>

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
    Catch_Time_Section("test2");
    //---------------------------------------

    // Start Chrono using MPI_Wtime
    double start_time = MPI_Wtime();

    // Print the Hello message and the rank for each proc
    std::cout << "Hello, World! from process " << rank << " out of " << size << std::endl;

    // Stop chrono
    double end_time = MPI_Wtime();

    // Print elapsed time using MPI_Wtime
    std::cout << "(MPI_Wtime)time spent for rank " << rank << " : " << (end_time - start_time)
              << std::endl;
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
