/**
 * @file Timer.hpp
 * @author (mouad.Bakhkhakh@cea.fr)
 * @brief Class Timer
 * @version 0.1
 * @date 2024-06-07
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <chrono>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#pragma once
/**
 * @brief Timer Class
 *
 */
class Timer {
 private:
  std::string name_;  // name of the function

  std::chrono::time_point<std::chrono::high_resolution_clock> start_utime_;
  std::chrono::duration<double> accum_time_;

 public:
  explicit Timer(std::string name);

  static bool Timer_enabled_;
  std::vector<std::tuple<int, std::string>> indentation_list_;
  int indentation_level_;

  void indentation();
  void lose_indentation(std::string func);
  inline void start();
  inline double stop(bool result_time);
  std::vector<std::tuple<int, std::string>>* return_indentation_list();

  ~Timer();
};

/**
 * @brief Global method used to get a Timer object
 *
 * @return Timer*
 */
inline Timer* get_indentation() {
  static Timer indent("");
  return &indent;
}

/**
 * @brief Global method used to increase the level of indentation by 1
 *
 */
void get_indent() {
  auto indent = get_indentation();
  indent->indentation();
}
/**
 * @brief Construct a new Timer:: Timer object
 *
 * @param name
 */
Timer::Timer(std::string name) : accum_time_(0.0), name_(name) {}

/**
 * @brief Increase the indentation level
 *
 */
void Timer::indentation() { this->indentation_level_++; }

/**
 * @brief Decrease indentation level
 *
 * @param func
 */
void Timer::lose_indentation(std::string func) {
  bool inside = false;
  for (int i = 0; i < this->indentation_list_.size(); i++) {
    if (std::get<1>(this->indentation_list_[i]) == func) {
      inside = true;
      break;
    }
  }
  if (inside == false) {
    this->indentation_list_.push_back(std::make_tuple(this->indentation_level_, func));
  }
  this->indentation_level_--;
}

/**
 * @brief Return list of indentation
 *
 * @return std::vector<std::tuple<int, std::string>>*
 */
std::vector<std::tuple<int, std::string>>* Timer::return_indentation_list() {
  auto indent = get_indentation();
  return &indentation_list_;
}

/**
 * @brief Method used to start chrono
 *
 */
inline void Timer::start() {
  if (this->Timer_enabled_ == false) {
    return;
  }
  get_indent();
  this->start_utime_ = std::chrono::high_resolution_clock::now();
}

/**
 * @brief Method used to get the value of the elapsed time
 *
 * @param result_time
 * @return double
 */
inline double Timer::stop(bool result_time = false) {
  if (this->Timer_enabled_ == false) {
    return 0.0;
  }
  const auto& utime = std::chrono::high_resolution_clock::now();

  if (result_time == false) {
    this->accum_time_ =
        std::chrono::duration_cast<std::chrono::duration<double>>(utime - this->start_utime_);
  }

  return this->accum_time_.count();
}

/**
 * @brief Destroy the Timer:: Timer object
 *
 */
Timer::~Timer() {}

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

/**
 * @brief Global method used to decrease the level of indentation by 1
 *
 * @param name_func
 */
void get_lose_indent(std::string name_func) {
  auto indent = get_indentation();
  indent->lose_indentation(name_func);
}

/**
 * @brief Get the indentation list
 *
 * @return std::vector<std::tuple<int, std::string>>*
 */
static std::vector<std::tuple<int, std::string>>* get_indentation_list() {
  auto indent = get_indentation();
  return indent->return_indentation_list();
}

bool Timer::Timer_enabled_ = false;
