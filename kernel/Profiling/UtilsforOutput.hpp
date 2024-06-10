/**
 * @file UtilsforOutput.hpp
 * @author (mouad.Bakhkhakh@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-07
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <mpi.h>

#include <Profiling/output.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#pragma once

/**
 * @brief Useful methods for managing outputs
 *
 */
class UtilsForOutput {
 public:
  static UtilsForOutput& getInstance() {
    static UtilsForOutput instance;
    return instance;
  }

  Output* get_timers();
  void get_enableOutput();
  void get_disableOutput();
  void update_timer(std::string name, Timers& timer);
  void savefiles();
  void print_timetable();
  void get_indentation();
  void get_lose_indentation(std::string func);
};

/**
 * @brief Get the timers object
 *
 * @return Output*
 */
Output* UtilsForOutput::get_timers() {
  static Output timers;
  return &timers;
}

/**
 * @brief Method allows to enable profiling
 *
 */
void UtilsForOutput::get_enableOutput() {
  auto timers = get_timers();
  timers->enableOutput();
}

/**
* @brief Method allows to enable profiling
* 
*/
void UtilsForOutput::get_disableOutput()
{
auto timers = get_timers();
timers->disableOutput();
}

/**
 * @brief Method used to save the results of a timer in the tuple type list : timerlist
 *
 * @param name
 * @param timer
 */
void UtilsForOutput::update_timer(std::string name, Timers& timer) {
  auto timers = get_timers();
  timers->save(name, timer);
}

/**
 * @brief Method used to save timerlist and files
 *
 */
void UtilsForOutput::savefiles() {
  auto timers = get_timers();
  timers->savetimerfile();
}

/**
 * @brief Method used to print timetable in terminal
 *
 */
void UtilsForOutput::print_timetable() {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  auto timers = get_timers();
  timers->timetable(rank, size);
}

/**
 * @brief Method used to add identation to a line in the timetable in order to have hierarchy
 *
 */
void UtilsForOutput::get_indentation() {
  auto timers = get_timers();
  timers->indentation();
}

/**
 * @brief Method used to remove identation from a line in the timetable in order to have hierarchy
 *
 * @param func
 */
void UtilsForOutput::get_lose_indentation(std::string func) {
  auto timers = get_timers();
  timers->lose_indentation(func);
}
