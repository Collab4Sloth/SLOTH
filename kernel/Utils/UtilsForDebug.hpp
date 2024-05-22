/**
 * @file UtilsForDebug.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-05-21
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <sys/resource.h>

#include <string>

#include "mfem.hpp"

#ifndef UTILSFORDEBUG_H_
#define UTILSFORDEBUG_H_

/**
 * @brief
 *
 */
class UtilsForDebug {
 public:
  static rusage make_memory_checkpoint();

  static void memory_checkpoint(const std::string& msg);
};

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
  std::cout << "Memory footprint " << mem.ru_maxrss * 1e-6 << " GB" << std::endl;
}
#endif

/* UTILSFORDEBUG_H_ */
