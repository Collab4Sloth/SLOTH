/**
 * @file main.cpp
 * @author MB278596 (mouad.bakhkhakh@cea.fr)
 * @brief This example calculates the execution time of Timers functions using a test code.
 * @version 0.1
 * @date 2024-06-10
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>

#include "Profiling/Profiling.hpp"

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
    Catch_Time_Section("test3");
    //---------------------------------------

    int n = 0;
    const int nmax = 1e3;  // 1e8 initially but reduced for integration
    for (int i = 0; i < nmax; i++) {
      Catch_Time_Section("somme");
      n++;
    }
    std::cout << "Test code executed successfully" << std::endl;
    std::cout << "n = " << n << std::endl;
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
