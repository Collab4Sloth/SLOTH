/**
 * @file output.hpp
 * @author (mouad.Bakhkhakh@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-07
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <mpi.h>

#include <Profiling/timers.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#pragma once

/**
 * @brief TimerTuple struct to store profiling information for a function
 *
 */
struct TimerTuple {
  int order;
  std::string function_name;
  int call_nmbrCalls;
  double total_time;

  /**
   * @brief Comparison operator to allow sorting TimerTuple objects by function name
   *
   * @param other
   * @return true
   * @return false
   */
  bool operator<(const TimerTuple& other) const { return function_name < other.function_name; }
};

/**
 * @brief Class for defining the the object Output used to print and print the profiling
 *
 */
class Output {
 private:
  std::string name;
  int order_;
  int level_;

 public:
  Output();
  explicit Output(std::string name);

  std::vector<TimerTuple> timerlist;
  std::vector<std::tuple<int, std::string>> levellist;
  std::vector<std::tuple<int, std::string>> indentationlist;
  int indentation_level_;
  bool Output_enabled_;

  void enableOutput();
  void disableOutput();

  void save(std::string name_func, Timers& func);
  void timetable(int rank, int size);
  void savetimerfile();
  void cleartimerfile();
  void indentation();
  void lose_indentation(std::string func);

  ~Output();
};

/**
 * @brief Method used to print lines in timetable
 *
 * @param c
 * @param length
 * @return std::string
 */
std::string line(char c, size_t length) { return std::string(length, c); }

/**
 * @brief Construct a new Output:: Output object
 *
 */
Output::Output() {}

/**
 * @brief Construct a new Output:: Output object
 *
 * @param name
 * @param rank
 */
Output::Output(std::string name) : order_(0), level_(0), name(name), Output_enabled_(true) {
  cleartimerfile();
}

/**
 * @brief Method to enable profiling
 *
 */
void Output::enableOutput() { 
  this->Output_enabled_ = true;
  get_enableTimers();
}

/**
 * @brief Method to disable profiling
 *
 */
void Output::disableOutput() { 
  this->Output_enabled_ = false; 
  get_disableTimers();
}

/**
 * @brief Method used to save the results of a timer in the tuple type list : timerlist
 *
 * @param name_func
 * @param func
 */
void Output::save(std::string name_func, Timers& func) {
  if (this->Output_enabled_ == false) {
    return;
  }

  bool inside = false;
  for (auto& timer : timerlist) {
    if (timer.function_name == name_func) {
      // increment the number of calls
      timer.call_nmbrCalls++;
      // update the total time
      timer.total_time += func.stop(true);
      inside = true;
      break;
    }
  }

  // If the timer does not exist, then add it to timerlist
  if (!inside) {
    TimerTuple new_timer{this->order_++, name_func, 1, func.stop(true)};
    timerlist.push_back(new_timer);
  }

  get_lose_indent(name_func);
}

/**
 * @brief Increases the indentation level if output is enabled
 *
 */
void Output::indentation() {
  if (this->Output_enabled_ == false) {
    return;
  }

  this->indentation_level_++;
}

/**
 * @brief Decrease the indentation level if output is disabled
 *
 * @param func
 */
void Output::lose_indentation(std::string func) {
  if (this->Output_enabled_ == false) {
    return;
  }

  bool inside = false;
  for (int i = 0; i < indentationlist.size(); i++) {
    if (std::get<1>(indentationlist[i]) == func) {
      inside = true;
      break;
    }
  }
  if (inside == false) {
    indentationlist.push_back(std::make_tuple(this->indentation_level_, func));
  }
  this->indentation_level_--;
}

/**
 * @brief Method used to save timerlist in the file timerlist.txt
 *
 */
void Output::savetimerfile() {
  if (this->Output_enabled_ == false) {
    return;
  }
  // name of the file
  std::string filename = "timerlist.txt";

  // open file to add data
  std::ofstream outputFile(filename, std::ios::app);
  if (!outputFile.is_open()) {
    std::cerr << "Error opening the file timerlist.txt for writing" << std::endl;
    return;
  }

  // Traverse the vector of timers and write to the file.
  for (const auto& timer : timerlist) {
    outputFile << timer.function_name << " " << timer.call_nmbrCalls << " " << timer.total_time
               << std::endl;
  }

  // close the file
  outputFile.close();
}

/**
 * @brief Method used to clear timerlist.txt
 *
 */
void Output::cleartimerfile() {
  if (this->Output_enabled_ == false) {
    return;
  }

  std::string filename = "timerlist.txt";

  std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
  if (!ofs.is_open()) {
    std::cerr << "Error: Unable to open the file timerlist.txt" << std::endl;
    return;
  }
  ofs.close();
}

/**
 * @brief Method used to print timetable in terminal
 *
 * @param rank
 * @param size
 */
void Output::timetable(int rank = 0, int size = 1) {
  if (this->Output_enabled_ == false) {
    return;
  }

  std::string filename = "timerlist.txt";  // name of the saving file

  std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);
  if (!outputFile.is_open()) {
    std::cerr << "Erreur lors de l'ouverture du fichier pour écriture!" << std::endl;
    return;
  }

  // looking for the maximum time in order to be used in time ratio
  double glob_time, max_timer;

  // sequentiel case
  if (size == 1) {
    max_timer = (timerlist[0]).total_time;
    for (int i = 1; i < timerlist.size(); i++) {
      if (max_timer < (timerlist[i]).total_time) {
        max_timer = (timerlist[i]).total_time;
      }
    }
  } else {
    // parallel case
    MPI_Reduce(&(timerlist[0]).total_time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      max_timer = glob_time / size;
    }
    for (int i = 1; i < timerlist.size(); i++) {
      MPI_Reduce(&(timerlist[i]).total_time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (rank == 0) {
        if (max_timer < (glob_time / size)) {
          max_timer = glob_time / size;
        }
      }
    }
  }

  // sequentiel case
  if (size == 1) {
    std::cout << " |" << line('-', 4) << "Start Timetable" << line('-', 123) << "|" << '\n';
    std::cout << " |" << line(' ', 38) << "Name" << line(' ', 40)
              << "| Number of Calls   |      Time(s)      |  time ratio (%)   |" << std::endl;
    std::cout << " |" << line('-', 142) << "|" << '\n';

    for (int j = timerlist.size() - 1; j >= 0; j--) {
      const TimerTuple& timer = timerlist[j];

      std::vector<std::tuple<int, std::string>>* indent_list = get_indentation_list();

      int index = 0;
      for (const auto& indent : *indent_list) {
        if (std::get<1>(indent) == timer.function_name) {
          index = std::get<0>(indent);
          break;
        }
      }

      int space = 80;
      std::cout << " | ";

      // add indentation
      for (int k = 0; k < index; k++) {
        std::cout << "-->";
        space -= 3;
      }

      // print results
      std::cout << std::setw(space) << std::left << timer.function_name << " | ";
      std::cout << std::setw(17) << std::left << timer.call_nmbrCalls << " | ";
      std::cout << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                << timer.total_time << " | ";
      std::cout << std::setw(17) << std::left << std::fixed << std::setprecision(2)
                << (timer.total_time * 100) / max_timer << " | " << std::endl;
    }
    std::cout << " |-- End Timetable " << line('-', 125) << "|" << std::endl;

  } else {
    // parallel case
    if (rank == 0) {
      std::cout << " |" << line('-', 4) << "Start Timetable" << line('-', 163) << "|" << '\n';
      std::cout << " |" << line(' ', 28) << "name" << line(' ', 50)
                << "| number of Calls   |    time Max (s)   |     mean time(s)  |     Imbalance    "
                   " | time ratio (%)    |"
                << std::endl;
      std::cout << " |" << line('-', 182) << "|" << '\n';
    }

    for (int j = timerlist.size() - 1; j >= 0; j--) {
      const TimerTuple& timer = timerlist[j];

      // Trouver l'index d'indentation
      std::vector<std::tuple<int, std::string>>* indent_list = get_indentation_list();

      int index = 0;
      for (const auto& indent : *indent_list) {
        if (std::get<1>(indent) == timer.function_name) {
          index = std::get<0>(indent);
          break;
        }
      }
      int global_calls;
      MPI_Reduce(&timer.call_nmbrCalls, &global_calls, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      double global_time;
      MPI_Reduce(&timer.total_time, &global_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      double Max_time;
      MPI_Reduce(&timer.total_time, &Max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      double Min_time;
      MPI_Reduce(&timer.total_time, &Min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

      if (rank == 0) {
        int space = 80;
        std::cout << " | ";

        for (int k = 0; k < index; k++) {
          std::cout << "-->";
          space -= 3;
        }
        // Name
        std::cout << std::setw(space) << std::left << timer.function_name << " | ";
        // Number of Calls
        std::cout << std::setw(17) << std::left << global_calls << " | ";
        // time_Max(s)
        std::cout << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                  << Max_time << " | ";
        // time_Mean(s)
        std::cout << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                  << global_time / size << " | ";
        // Imbalance
        std::cout << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                  << ((Max_time - (global_time / size)) / (global_time / size)) * 100 << " | ";
        // time ratio (%)
        std::cout << std::setw(17) << std::left << std::fixed << std::setprecision(2)
                  << ((global_time / size) * 100) / max_timer << " | " << std::endl;

        // save the results in the file timerlist.txt
        for (int k = 0; k < index; k++) {
          outputFile << "   ";
          space -= 3;
        }
        outputFile << timer.function_name << " " << global_calls << " " << global_time / size
                   << std::endl;
      }
    }
    if (rank == 0) {
      outputFile.close();
      std::cout << " |-- End Timetable " << line('-', 165) << "|" << std::endl;
    }
  }
}

Output::~Output() {}
