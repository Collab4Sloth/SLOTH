/**
 * @file timers.hpp
 * @author (mouad.Bakhkhakh@cea.fr)
 * @brief Class Timers
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
 * @brief Timers Class
 *
 */
class Timers {
 private:
  std::string name_;  // name of the function

  std::chrono::time_point<std::chrono::high_resolution_clock> start_utime_;
  std::chrono::duration<double> accum_time_;

 public:
  explicit Timers(std::string name);

  static bool timers_enabled_;
  std::vector<std::tuple<int, std::string>> indentation_list_;
  int indentation_level_;

  void indentation();
  void lose_indentation(std::string func);
  inline void start();
  inline double stop(bool result_time);
  std::vector<std::tuple<int, std::string>>* return_indentation_list();

  ~Timers();
};

/**
 * @brief Global method used to get a timers object
 *
 * @return Timers*
 */
inline Timers* get_indentation() {
  static Timers indent("");
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
 * @brief Construct a new Timers:: Timers object
 *
 * @param name
 */
Timers::Timers(std::string name) : accum_time_(0.0), name_(name) {}

/**
 * @brief Increase the indentation level
 *
 */
void Timers::indentation() { this->indentation_level_++; }

/**
 * @brief Decrease indentation level
 *
 * @param func
 */
void Timers::lose_indentation(std::string func) {
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
std::vector<std::tuple<int, std::string>>* Timers::return_indentation_list() {
  auto indent = get_indentation();
  return &indentation_list_;
}

/**
 * @brief Method used to start chrono
 *
 */
inline void Timers::start() {
  if (this->timers_enabled_ == false) {
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
inline double Timers::stop(bool result_time = false) {
  if (this->timers_enabled_ == false) {
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
 * @brief Destroy the Timers:: Timers object
 *
 */
Timers::~Timers() {}

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

bool Timers::timers_enabled_ = false;

/**
 * @brief Method allows to enable Timers
 *
 */
void get_enableTimers() { Timers::timers_enabled_ = true; }

/**
 * @brief Method allows to disable Timers
 *
 */
void get_disableTimers() { Timers::timers_enabled_ = false; }
