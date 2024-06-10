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

#include "Profiling/UtilsforOutput.hpp"
#include "Profiling/output.hpp"
#include "Profiling/timers.hpp"

int main(int argc, char** argv) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // ------Start profiling-------------------------

  Timers timer0("timer0");
  timer0.start();

  ///////////// Test code /////////////////////////
  // ------Start profiling of the test code-------------------------
  Output output3("output3");

  // --Enable profiling--
  UtilsForOutput::getInstance().get_enableOutput();

  // --Disable profiling--
  // UtilsForOutput::getInstance().get_disableOutput();
  // ---------------------------------------------
  Timers timer_ex3("ex3");
  timer_ex3.start();

  int n = 0;
  for (int i = 0; i < 100000000; i++) {
    n++;
  }
  std::cout << "Test code executed successfully" << std::endl;
  std::cout << "n = " << n << std::endl;

  timer_ex3.stop();

  // ---------End Profiling of the test code-------------------
  UtilsForOutput::getInstance().update_timer("ex3", timer_ex3);
  UtilsForOutput::getInstance().print_timetable();
  UtilsForOutput::getInstance().savefiles();
  //----------------------------------------

  // End Profiling and calculate the elapsed time of Timers functions

  timer0.stop();
  std::cout << "temps d'exécution du timer est : " << timer0.stop(true) - timer_ex3.stop(true)
            << std::endl;

  MPI_Finalize();
  return 0;
}