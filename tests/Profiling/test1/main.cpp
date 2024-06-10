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

#include <chrono>
#include <iomanip>
#include <iostream>

#include "Profiling/UtilsforOutput.hpp"
#include "Profiling/output.hpp"
#include "Profiling/timers.hpp"

int main(int argc, char** argv) {
  // Initialize MPI
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ------Start profiling------------------------
  Output output1("output1");

  // --Enable profiling--
  UtilsForOutput::getInstance().get_enableOutput();
  // --Disable profiling--
  // UtilsForOutput::getInstance().get_disableOutput();

  Timers timer_ex1("ex1");
  Timers timer_somme_locale("somme_locale");
  timer_ex1.start();
  //---------------------------------------------

  int n = 1e6;
  int somme_locale = 0;
  int somme_globale = 0;

  // Compute the local sum for each proc
  for (int i = rank + 1; i <= n; i += size) {
    timer_somme_locale.start();
    somme_locale += i;
    timer_somme_locale.stop();
    UtilsForOutput::getInstance().update_timer("somme_locale", timer_somme_locale);
  }

  // Reduction of the local sum
  MPI_Reduce(&somme_locale, &somme_globale, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "the sum of numbers from 1 to  " << n << " is : " << somme_globale << std::endl;
  }

  // ---------End Profiling------------------
  timer_ex1.stop();
  UtilsForOutput::getInstance().update_timer("ex1", timer_ex1);
  UtilsForOutput::getInstance().print_timetable();
  UtilsForOutput::getInstance().savefiles();
  //----------------------------------------
  MPI_Finalize();

  return 0;
}