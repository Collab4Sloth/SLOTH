/**
 * @file MATimersAPI.hxx
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
#include <MAToolsProfiling/MATimers.hxx>
#include <MAToolsProfiling/MATimersVerbosity.hxx>
#include <string>
// One level of macro indirection is required in order to resolve __COUNTER__,
// and get varname1 instead of varname__COUNTER__.
#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a##b

#define CONCAT_LINE(x) CONCAT(x, __LINE__)
#define CONCAT_FILE(x) CONCAT(x, __FILE__)
#define CONCAT_COUNTER(x) CONCAT(x, __COUNTER__)
#define VARNAME() CONCAT_COUNTER(BASE)

inline std::chrono::duration<double>* get_duration(std::string a_name) {
  auto& ptr = MATools::MATimer::get_MATimer_node<CURRENT>();
  if (MATools::MATimer::is_enable()) {
    assert(ptr != nullptr && "do not use an undefined MATimerNode");
    ptr = ptr->find(a_name);
    ptr->update_count();
    MATools::MATimer::print_verbosity_level_1(a_name, ptr->get_level());
    return ptr->get_ptr_duration();
  }
  return nullptr;
}

#define Catch_Time_Section(XNAME) MATools::MATimer::Timer VARNAME()(get_duration(XNAME));

// alias
#define Catch_Nested_Time_Section(XNAME) Catch_Time_Section(XNAME)
#define START_TIMER(XNAME) Catch_Time_Section(XNAME)

/**
 * @brief This function captures the runtime of a given section.
 * @param [in] lambda section that the user wants to measure.
 * @return The runtime of lambda
 */
template <typename Lambda>
double chrono_section(Lambda&& lambda_function) {
  using namespace MATools::MATimer;
  double ret;
  BasicTimer timer;
  timer.start();
  lambda_function();
  timer.end();
  ret = timer.get_duration();
  return ret;
}

/**
 * @brief This function captures the runtime of a given section and add it in the current
 * MATimerNode named a_name.
 * @param [in] a_name of the chrono section measured.
 * @param [in] a_lambda is the section captured.
 */
template <typename Lambda>
void add_capture_chrono_section(std::string a_name, Lambda&& a_lambda_function) {
  Catch_Time_Section(a_name);
  a_lambda_function();
}

/** Create an object MATimersManager:
 * To enhance your application's functionality, you need to create an object called
 * "MATimersManager" and place it at the beginning of your application. The MATimersManager serves
 * as a tool for monitoring and analyzing the performance of your application. By initializing it at
 * the outset, you can ensure that all relevant metrics and data are captured from the start and
 * throughout the execution of your application. This allows for comprehensive profiling and
 * optimization of your application's performance.
 */
class MATimersManager {
 public:
  /**
   * @brief Constructor for MATimersManager.
   * Initializes the MATimersManager by initializing the timer.
   */
  MATimersManager() { MATools::MATimer::initialize(); }

  /**
   * @brief Constructor for MATimersManager.
   * Initializes the MATimersManager by initializing the timer.
   */
  MATimersManager([[maybe_unused]] int* argc, [[maybe_unused]] char*** argv) {
    // MATools::MPI::mpi_initialize(argc,argv);
    MATools::MATimer::initialize();
  }

  /**
   * @brief Disables printing of the timetable.
   * This function disables the printing of the timetable during profiling.
   */
  void disable_timetable() { MATools::MATimer::Optional::disable_print_timetable(); }

  /**
   * @brief Disables writing data to a file.
   * This function disables writing profiling data to a file.
   */
  void disable_write_file() { MATools::MATimer::Optional::disable_write_file(); }

  void enable() {
    MATools::MATimer::finalize();  // reset if necessary
    MATools::MATimer::initialize();
  }

  void print() { MATools::MATimer::print(); }

  /**
   * @brief Destructor for MATimersManager.
   * Finalizes the MATimersManager by cleaning up the timer resources.
   */
  ~MATimersManager() {
    MATools::MATimer::finalize();
    // MATools::MPI::mpi_finalize();
  }
};

namespace Profiling {
MATimersManager& getInstance() {
  static MATimersManager manager;
  return manager;
}
}  // namespace Profiling
