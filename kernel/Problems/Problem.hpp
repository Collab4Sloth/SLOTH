/*
 * Copyright © CEA 2023
 *
 * \brief Phase-Field problems
 *
 * \file Problem.hpp
 * \author ci230846
 * \date 31/03/2023
 */

#pragma once
#include "../BCs/BoundaryConditions.hpp"
#include "../Operators/PhaseFieldOperator.hpp"
#include "../Parameters/Parameter.hpp"
#include "../PostProcessing/postprocessing.hpp"
#include "../Time/Time.hpp"
#include "../Utils/PhaseFieldOptions.hpp"
#include "../Variables/Variable.hpp"
#include "../Variables/Variables.hpp"
#include "mfem.hpp"

class Problem {
 private:
  SpatialDiscretization<T, DIM> spatial_;
  PhaseFieldOperator<T, DIM> operator_;
  PostProcessing<T, DC, DIM> postprocessing_;
  TimeDiscretization<T, DC, DIM> time_;

 public:
  Problem();
  virtual void initialize() = 0;
  virtual void solve() = 0;
  void check_before_solve();
  void solve();
  ~Problem();
};

/**
 * @brief Construct a new Problem:: Problem object
 *
 */
template <class T, class DC, int DIM>
Problem<T, DC, DIM>::Problem() {}

/**
 * @brief check the existance of mandatory objects
 *
 */
template <class T, class DC, int DIM>
void Problem<T, DC, DIM>::check_before_solve() {
  this->spatial_ = this->load_spatial();
  if (this->spatial_ == NULL) {
    throw std::runtime_error("Error while loading spatial object");
  }
  this->operator_ = this->load_operator();
  if (this->operator_ == NULL) {
    throw std::runtime_error("Error while loading operator object");
  }
  this->postprocessing_ = this->load_postprocessing();
  if (this->postprocessing_ == NULL) {
    throw std::runtime_error("Error while loading postprocessing object");
  }
  this->time_ = this->load_time();
  if (this->time_ == NULL) {
    throw std::runtime_error("Error while loading time object");
  }
}

template <class T, class DC, int DIM>
SpatialDiscretization<T, DIM>* Problem<T, DC, DIM>::load_spatial() {
  return nullptr;
}

template <class T, class DC, int DIM>
PhaseFieldOperator<T, DIM>* Problem<T, DC, DIM>::load_operator() {
  return nullptr;
}
template <class T, class DC, int DIM>
PostProcessing<T, DC, DIM>* Problem<T, DC, DIM>::load_postprocessing() {
  return nullptr;
}
template <class T, class DC, int DIM>
TimeDiscretization<T, DC, DIM>* Problem<T, DC, DIM>::load_time() {
  return nullptr;
}

/**
 * @brief check the existance of mandatory objects and solve the problem
 *
 */
template <class T, class DC, int DIM>
void Problem<T, DC, DIM>::solve() {
  this->check_before_solve();
  this->time_.execute();
}

/**
 * @brief Destroy the Problem:: Problem object
 *
 */
template <class T, class DC, int DIM>
Problem<T, DC, DIM>::~Problem() {}
