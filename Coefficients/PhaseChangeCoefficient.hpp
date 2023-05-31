/*
 * Copyright © CEA 2023
 *
 * PhaseChangeCoefficient.hpp
 *
 *  Created on: 15 may 2023
 *      Author: ci230846
 */
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"
#pragma once

// TODO(ci) : il faut prévoir plusieurs constructeur pour anticiper un changement de phase à la
// Welland, Calphad, métamodèle...., constant

class PhaseChangeCoefficient : public mfem::Coefficient {
 private:
  template <class... Args>
  std::function<double(const double &)> getConstant(Args... args);

  template <class... Args>
  std::function<double(const double &)> getCalphad(Args... args);

  const mfem::GridFunction gf_;

 public:
  PhaseChangeCoefficient();

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~PhaseChangeCoefficient();
};

/**
 * @brief Construct a new  PhaseChange Coefficient:: Source Term Coefficient object
 *
 * @tparam PHASE
 */

PhaseChangeCoefficient::PhaseChangeCoefficient() {}

double PhaseChangeCoefficient::Eval(mfem::ElementTransformation &T,
                                    const mfem::IntegrationPoint &ip) {
  // const auto xx = this->gf_.GetValue(T.ElementNo, ip);

  const auto phase_change_value = 0.;  // this->phase_change_ * xx;

  return phase_change_value;
  // switch (PHASE) {
  //   case PhaseChange::Null: {
  //     const auto phase_change_value = 0.;
  //     return phase_change_value;
  //   }

  //   case PhaseChange::Constant: {
  //     const auto phase_change_value = 0.;
  //     return phase_change_value;
  //   }
  //   case PhaseChange::Calphad: {
  //     throw std::runtime_error("PhaseChangeCoefficient::Eval: Calphad is not implemented yet");
  //     break;
  //   }
  //   default:
  //     throw std::runtime_error(
  //         "SourceTermCoefficient::Eval: only Null, Sinusoide2D source term are available");
  //     break;
  // }
}

/**
 * @brief Destroy the Phase Change Coefficient:: Phase Change Coefficient object
 *
 * @tparam PHASE
 */

PhaseChangeCoefficient::~PhaseChangeCoefficient() {}
