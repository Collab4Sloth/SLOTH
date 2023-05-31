/*
 * Copyright © CEA 2023
 *
 * \brief Analytical Mobilities
 *
 * \file PhaseFieldMobilities.hpp
 * \author ci230846
 * \date 15/05/2023
 */

#ifndef PHASEFIELDMOBILITIES_HPP_
#define PHASEFIELDMOBILITIES_HPP_
#pragma once
#include <functional>
#include <string>
#include <vector>
#include "Utils/PhaseFieldOptions.hpp"
// TODO(ci): voir si on fait une mobilité degénérée implicite...

template <Mobility MOBI>
class MobilityFunctions {
 private:
  template <class... Args>
  std::function<double(const double&)> getConstant(Args... args);

  template <class... Args>
  std::function<double(const double&)> getDegenerated(Args... args);

 public:
  MobilityFunctions();

  template <class... Args>
  std::function<double(const double&)> getMobilityFunction(Args... args);

  ~MobilityFunctions();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

//  explicit instantiation of AnalyticalFunctions class template/
template class MobilityFunctions<Mobility::Constant>;
template class MobilityFunctions<Mobility::Degenerated>;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new potential function:: potential function object
 *
 */
template <Mobility MOBI>
MobilityFunctions<MOBI>::MobilityFunctions() {}

template <Mobility MOBI>
template <typename... Args>
std::function<double(const double&)> MobilityFunctions<MOBI>::getConstant(Args... args) {
  return std::function<double(const double&)>([](double x) { return x; });
}

/**
 * @brief (explicit) degenerated mobility : mob_cst xn*xn * (1-xn)*(1-xn)
 *
 * @tparam MOBI
 * @tparam Args
 * @param args
 * @return std::function<double(const double&)>
 */
template <Mobility MOBI>
template <typename... Args>
std::function<double(const double&)> MobilityFunctions<MOBI>::getDegenerated(Args... args) {
  auto v = std::vector<double>{args...};

  if (v.size() == 1) {
    const auto xn = v[0];
    return std::function<double(const double&)>([xn](double mob_cst) {
      const auto mob = mob_cst * xn * xn * (1.0 - xn) * (1.0 - xn);
      return mob;
    });
  } else {
    throw std::runtime_error(
        "MobilityFunctions::getDegenerated: only one argument is expected for explicit mobility");
  }
}

/**
 * @brief return the function associated with the potential_name and its ORDER of
 * derivative
 *
 * @param potential_name
 * @return const double
 */
template <Mobility MOBI>
template <class... Args>
std::function<double(const double&)> MobilityFunctions<MOBI>::getMobilityFunction(Args... args) {
  switch (MOBI) {
    case Mobility::Constant:
      return this->getConstant(args...);
    case Mobility::Degenerated:
      return this->getDegenerated(args...);
    default:
      throw std::runtime_error(
          "MobilityFunctions::getMobilityFunctions: only constant and degenerated mobilities  are "
          "available");
      break;
  }
}

/**
 * @brief Destroy the potential function :: potential function  object
 *
 */
template <Mobility MOBI>
MobilityFunctions<MOBI>::~MobilityFunctions() {}

#endif /* PHASEFIELDMobilities_HPP_ */
