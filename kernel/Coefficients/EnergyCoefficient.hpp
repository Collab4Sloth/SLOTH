/**
 * @file EnergyCoefficient.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-18
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <numeric>

#include "mfem.hpp"
#pragma once

//--------------------------
//--------------------------
//--------------------------
class InterfacialCoefficient : public mfem::Coefficient {
 private:
  mfem::GridFunction *gfu;
  double lambda;

 public:
  InterfacialCoefficient(mfem::GridFunction *gfu_, const double &lambda_)
      : gfu(gfu_), lambda(lambda_) {}
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

double InterfacialCoefficient::Eval(mfem::ElementTransformation &T,
                                    const mfem::IntegrationPoint &ip) {
  mfem::Vector gradu;
  gfu->GetGradient(T, gradu);
  const auto &value = std::inner_product(gradu.begin(), gradu.end(), gradu.begin(), 0.);
  return lambda * value;
}

//--------------------------
//--------------------------
//--------------------------
class EnergyCoefficient : public mfem::Coefficient {
 private:
  mfem::GridFunction *gfu;
  double lambda;
  double omega;

 public:
  EnergyCoefficient(mfem::GridFunction *gfu_, const double &lambda_, const double &omega_)
      : gfu(gfu_), lambda(lambda_), omega(omega_) {}
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

double EnergyCoefficient::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) {
  mfem::Vector gradu;
  gfu->GetGradient(T, gradu);
  const auto phi = gfu->GetValue(T.ElementNo, ip);

  const auto &value = std::inner_product(gradu.begin(), gradu.end(), gradu.begin(), 0.);
  const auto f_int = lambda * value;
  const auto f_o = omega * phi * phi * (1. - phi) * (1. - phi);
  const auto energy = f_o + f_int;
  return energy;
}