/**
 * @file UtilsForPostProcessing.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Usefull methods for post-processing
 * @version 0.1
 * @date 2025-09-05
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
