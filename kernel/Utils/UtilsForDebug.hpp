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
#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

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

/*!
 * \brief Lambda expression used to manage exception
 * \param[in] b: boolean variable used in a conditional test
 * \param[in] method: string variable specifying the method concerned by the
 * exception
 * \param[in] msg: string variable specifying the error message written on the
 * screen
 */
static auto throw_if = [](const bool b, const std::string& method, const std::string& msg) {
  if (b) {
    SlothInfo::error(method, ": ", msg);
    mfem::mfem_error("Error message");
  }
};
