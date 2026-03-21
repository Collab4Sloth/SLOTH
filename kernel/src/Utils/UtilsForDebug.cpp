/**
 * @file UtilsForDebug.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Usefull methods for debug
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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

#include "Utils/UtilsForDebug.hpp"

#include <sys/resource.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief This function creates a rusage variable (by rp)
 * @return rusage data type that get memory information
 */

rusage UtilsForDebug::make_memory_checkpoint() {
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto res = getrusage(who, &obj);
  assert((res = -1) && "error: getrusage has failed");
  return obj;
}

/**
 * @brief Get and Print memory information
 *
 * @param msg
 */
void UtilsForDebug::memory_checkpoint(const std::string& msg) {
  auto mem = UtilsForDebug::make_memory_checkpoint();
  std::cout << "========= MEMORY CHECKPOINT ====== " << std::endl;
  std::cout << "<<<" << msg << ">>>" << std::endl;
#if defined(__APPLE__)
  // macOS : ru_maxrss est en octets
  std::cout << "Memory footprint "
            << static_cast<double>(mem.ru_maxrss) / (1024.0 * 1024.0 * 1024.0) << " GB"
            << std::endl;
#else
  // Linux : ru_maxrss est en kB
  std::cout << "Memory footprint " << static_cast<double>(mem.ru_maxrss) / (1024.0 * 1024.0)
            << " GB" << std::endl;
#endif
}
