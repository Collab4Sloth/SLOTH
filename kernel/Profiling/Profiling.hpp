
/**
 * @file Profiling.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Macro used to catch timer section
 * @version 0.1
 * @date 2024-06-14
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <mpi.h>

#include <filesystem>  // NOLINT [avoid  <filesystem> is an unapproved C++17 header.]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "Profiling/Timer.hpp"

#pragma once

/**
 * @brief Method used to print lines in timetable
 *
 * @param c
 * @param length
 * @return std::string
 */
std::string line(char c, size_t length) { return std::string(length, c); }

/**
 * @brief TimerTuple struct to store profiling information for a function
 *
 */
struct TimerTuple {
  int order;
  std::string function_name;
  std::uint64_t call_nmbrCalls;
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
/////////////////////////////////////
////////////////////////////////////
/**
 * @brief Useful methods to manage profiling
 *
 */
class Profiling {
 private:
  std::string name;
  int order_{0};

  void save(std::string name_func, Timer& func);
  void get_indentation();
  void get_lose_indentation(std::string func);
  bool Profiling_enabled_{true};

  std::vector<TimerTuple> timerlist;
  std::vector<std::tuple<int, std::string>> levellist;
  std::vector<std::tuple<int, std::string>> indentationlist;

 public:
  static Profiling& getInstance() {
    static Profiling instance;
    return instance;
  }
  int indentation_level_;
  void enable();
  void update_timer(std::string name, Timer& timer);
  void print();
};

/**
 * @brief Method allows to enable profiling
 *
 */
void Profiling::enable() { Timer::Timer_enabled_ = true; }

/**
 * @brief Method used to save the results of a timer in the tuple type list : timerlist
 *
 * @param name_func
 * @param func
 */
void Profiling::save(std::string name_func, Timer& func) {
  if (this->Profiling_enabled_ == false) {
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
 * @brief Method used to save the results of a timer in the tuple type list : timerlist
 *
 * @param name
 * @param timer
 */
void Profiling::update_timer(std::string name, Timer& timer) { this->save(name, timer); }

/**
 * @brief Method used to print timetable in terminal
 *
 */
void Profiling::print() {
  if (this->Profiling_enabled_ == false) {
    return;
  }
  int rank = mfem::Mpi::WorldRank();
  int size = mfem::Mpi::WorldSize();

  std::filesystem::path filename = "profiling.txt";
  std::ofstream outputFile;

  // looking for the maximum time in order to be used in time ratio
  double glob_time, max_timer;

  // serial case
  if (size == 1) {
    max_timer = (timerlist[0]).total_time;
    for (std::size_t i = 1; i < timerlist.size(); i++) {
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
    for (std::size_t i = 1; i < timerlist.size(); i++) {
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
    outputFile.open(filename, std::ios::out | std::ios::trunc);
    if (!outputFile.is_open()) {
      std::cerr << "Error while opening " << filename << std::endl;
      return;
    }
    outputFile << " |" << line('-', 4) << "Start Timetable" << line('-', 123) << "|" << '\n';
    outputFile << " |" << line(' ', 38) << "Name" << line(' ', 40)
               << "| Number of Calls   |      Time(s)      |  time ratio (%)   |" << std::endl;
    outputFile << " |" << line('-', 142) << "|" << '\n';

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
      outputFile << " | ";

      // Add indentation
      for (int k = 0; k < index; k++) {
        outputFile << "---";
        space -= 3;
      }

      // print results
      outputFile << std::setw(space) << std::left << timer.function_name << " | ";
      outputFile << std::setw(17) << std::left << timer.call_nmbrCalls << " | ";
      outputFile << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                 << timer.total_time << " | ";
      outputFile << std::setw(17) << std::left << std::fixed << std::setprecision(2)
                 << (timer.total_time * 100) / max_timer << " | " << std::endl;
    }
    outputFile << " |-- End Timetable " << line('-', 125) << "|" << std::endl;
    outputFile.close();

  } else {
    // parallel case
    if (rank == 0) {
      outputFile.open(filename, std::ios::out | std::ios::trunc);
      if (!outputFile.is_open()) {
        std::cerr << "Error while opening " << filename << std::endl;
        return;
      }
      outputFile << " |" << line('-', 4) << "Start Timetable" << line('-', 163) << "|" << '\n';
      outputFile
          << " |" << line(' ', 28) << "name" << line(' ', 50)
          << "| number of Calls   |    time Max (s)   |     mean time(s)  |     Imbalance    "
             " | time ratio (%)    |"
          << std::endl;
      outputFile << " |" << line('-', 182) << "|" << '\n';
    }

    for (int j = timerlist.size() - 1; j >= 0; j--) {
      const TimerTuple& timer = timerlist[j];

      // Find index of indentation
      std::vector<std::tuple<int, std::string>>* indent_list = get_indentation_list();

      int index = 0;
      for (const auto& indent : *indent_list) {
        if (std::get<1>(indent) == timer.function_name) {
          index = std::get<0>(indent);
          break;
        }
      }
      double global_time;
      MPI_Reduce(&timer.total_time, &global_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      double Max_time;
      MPI_Reduce(&timer.total_time, &Max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      double Min_time;
      MPI_Reduce(&timer.total_time, &Min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      if (rank == 0) {
        int space = 80;
        outputFile << " | ";

        for (int k = 0; k < index - 1; k++) {
          outputFile << "---";
          space -= 3;
        }
        outputFile << "-->";
        space -= 3;
        // Name
        outputFile << std::setw(space) << std::left << timer.function_name << " | ";
        // Number of Calls
        outputFile << std::setw(17) << std::left << timer.call_nmbrCalls << " | ";
        // time_Max(s)
        outputFile << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                   << Max_time << " | ";
        // time_Mean(s)
        outputFile << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                   << global_time / size << " | ";
        // Imbalance
        outputFile << std::setw(17) << std::left << std::scientific << std::setprecision(10)
                   << ((Max_time - (global_time / size)) / (global_time / size)) * 100 << " | ";
        // time ratio (%)
        outputFile << std::setw(17) << std::left << std::fixed << std::setprecision(2)
                   << ((global_time / size) * 100) / max_timer << " | " << std::endl;
      }
    }
    if (rank == 0) {
      outputFile << " |-- End Timetable " << line('-', 165) << "|" << std::endl;
      outputFile.close();
    }
  }

  // Print on screen
  if (rank == 0) {
    std::ifstream inputFile(filename);

    if (!inputFile) {
      std::cerr << "Error while opening " << filename << " for reading." << std::endl;
      return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
      std::cout << line << std::endl;
    }
  }
}

/**
 * @brief Method used to add identation to a line in the timetable in order to have hierarchy
 *
 */
void Profiling::get_indentation() {
  if (this->Profiling_enabled_ == false) {
    return;
  }

  this->indentation_level_++;
}

/**
 * @brief Method used to remove identation from a line in the timetable in order to have hierarchy
 *
 * @param func
 */
void Profiling::get_lose_indentation(std::string func) {
  if (this->Profiling_enabled_ == false) {
    return;
  }

  bool inside = false;
  for (std::size_t i = 0; i < indentationlist.size(); i++) {
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
 * @brief Macro used to catch time of a section
 *
 */
#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a##b

#define CONCAT_LINE(x) CONCAT(x, __LINE__)
#define CONCAT_FILE(x) CONCAT(x, __FILE__)
#define CONCAT_COUNTER(x) CONCAT(x, __COUNTER__)
#define VARNAME() CONCAT_COUNTER(BASE)

struct APITimer {
  Timer m_timer;
  std::string m_name;

  explicit APITimer(std::string name) : m_timer(name), m_name(name) { m_timer.start(); }

  ~APITimer() {
    m_timer.stop();
    Profiling::getInstance().update_timer(m_name, m_timer);
  }
};

#define Catch_Time_Section(XNAME) APITimer VARNAME()(XNAME);
