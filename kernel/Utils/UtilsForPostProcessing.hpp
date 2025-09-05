/**
 * @file UtilsForPostProcessing.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Usefull methods for Postprocessing
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#pragma once

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Custom IterationKey for specialized map
 *
 */
struct IterationKey {
  std::pair<std::string, int> iter_;
  std::pair<std::string, double> time_step_;
  std::pair<std::string, double> time_;

  /**
   * @brief Construct a new Iteration Key object
   *
   * @param iter
   * @param time_step
   * @param time
   * @param s_iter
   * @param s_time_step
   * @param s_time
   */
  IterationKey(int iter, double time_step, double time, std::string s_iter = "Iter[-]",
               std::string s_time_step = "Dt[s]", std::string s_time = "Time[s]")
      : iter_(s_iter, iter), time_step_(s_time_step, time_step), time_(s_time, time) {}

  /**
   * @brief Comparison operator : mandatory for using IterationKey as a key in map
   *
   * @param user_key
   * @return true
   * @return false
   */
  bool operator<(const IterationKey& user_key) const {
    return std::tie(iter_, time_step_, time_) <
           std::tie(user_key.iter_, user_key.time_step_, user_key.time_);
  }
};
