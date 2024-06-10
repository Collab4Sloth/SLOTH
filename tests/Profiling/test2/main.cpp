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

#include "Profiling/UtilsforOutput.hpp"
#include "Profiling/output.hpp"
#include "Profiling/timers.hpp"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // Initialize MPI
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ------Start profiling-------------------------
  Output output2("output2");

  // --Enable profiling--
  UtilsForOutput::getInstance().get_enableOutput();

  // --Disable profiling--
  // UtilsForOutput::getInstance().get_disableOutput();

  Timers timer_ex2("ex2");
  timer_ex2.start();
  //---------------------------------------------

  // Start Chrono using MPI_Wtime
  double start_time = MPI_Wtime();

  // Print the Hello message and the rank for each proc
  std::cout << "Hello, World! from process " << rank << " out of " << size << std::endl;

  // Stop chrono
  double end_time = MPI_Wtime();
  timer_ex2.stop();

  // Print elapsed time using MPI_Wtime
  std::cout << "(MPI_Wtime)temps écoulé pour le rang " << rank << " : " << (end_time - start_time)
            << std::endl;

  // ---------End Profiling------------------
  UtilsForOutput::getInstance().update_timer("ex2", timer_ex2);
  UtilsForOutput::getInstance().print_timetable();
  UtilsForOutput::getInstance().savefiles();
  //----------------------------------------

  MPI_Finalize();

  return 0;
}
