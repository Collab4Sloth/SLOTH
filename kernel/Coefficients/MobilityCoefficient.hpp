/*
 * Copyright © CEA 2023
 *
 * MobilityCoefficient.hpp
 *
 *  Created on: 7 fev. 2023
 *      Author: ci230846
 */
#include <algorithm>

#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]
#pragma once

template <Mobility MOBI>
class MobilityCoefficient : public mfem::Coefficient {
 private:
  double mobility_;
  const mfem::GridFunction gf;
  int degenerate_order_;

 public:
  explicit MobilityCoefficient(const double &mob_c);
  MobilityCoefficient(const double &mob_c, mfem::GridFunction mob_gf, const int &order);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

/**
 * @brief Construct a new Mobility Coefficient< MOBI>:: Mobility Coefficient object
 *
 * @tparam MOBI
 * @param mob_c
 */
template <Mobility MOBI>
MobilityCoefficient<MOBI>::MobilityCoefficient(const double &mob_c) : mobility_(mob_c) {}

/**
 * @brief Construct a new Mobility Coefficient< MOBI>:: Mobility Coefficient object
 *
 * @tparam MOBI
 * @param mob_c
 * @param mob_gf
 * @param order
 */
template <Mobility MOBI>
MobilityCoefficient<MOBI>::MobilityCoefficient(const double &mob_c, mfem::GridFunction mob_gf,
                                               const int &order)
    : mobility_(mob_c), gf(mob_gf), degenerate_order_(order) {}

/**
 * @brief Evaluation of the mobility coefficient at integration point
 *
 * @tparam MOBI
 * @param T
 * @param ip
 * @return double
 */
template <Mobility MOBI>
double MobilityCoefficient<MOBI>::Eval(mfem::ElementTransformation &T,
                                       const mfem::IntegrationPoint &ip) {
  switch (MOBI) {
    case Mobility::Constant:
      return this->mobility_;

    case Mobility::Degenerated: {
      const auto xx = gf.GetValue(T.ElementNo, ip);
      const auto a_xx = std::max(std::min(xx, 1.), 0.);
      const auto degenerated_mob_func = this->mobility_ *
                                        std::pow(1. - a_xx, this->degenerate_order_) *
                                        std::pow(a_xx, this->degenerate_order_);
      return degenerated_mob_func;
    }

    default:
      throw std::runtime_error(
          "MobilityCoefficient::Eval: only constant and degenerated mobilities  are available");
      break;
  }
}
