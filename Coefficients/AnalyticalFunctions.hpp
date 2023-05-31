/*
 * Copyright © CEA 2022
 *
 * AnalyticalFunctions.h
 *
 *  Created on: 31 déc. 2021
 *      Author: ci230846
 */
#include <algorithm>
#include "mfem.hpp"
namespace analytical {

// Heaviside function centered at 0.5*r
double heavisideFunction(const mfem::Vector &x, const double &r) {
  if (x.Norml2() < 0.5 * r) {
    return 0.0;
  } else {
    return 1.0;
  }
}

// Heaviside function centered at 0.5*r
double heavisideFunction2(const mfem::Vector &x) {
  if (x[0] < 0.5) {
    return 0.0;
  } else {
    return 1.0;
  }
}

// Heaviside function centered at 0.5*r
double nullRHS(const mfem::Vector &x) {
  auto rhs = 0.;
  return rhs;
}

// Heaviside function centered at 0.5*r
double hyperbolicTangent(const mfem::Vector &x) {
  const auto a = 2.;
  const auto epsilon = 1.e-4;  // 4.e-4;
  const auto L = 1.e-4;        // 2.e-3;        // 1.e-3;
  const auto xx = x[0] - 1.e-3;
  const auto yy = x[1] - 4.e-4;
  const auto r = std::sqrt(xx * xx + yy * yy);
  return 0.5 + 0.5 * std::tanh(a * (r - 0.5 * L) / epsilon);
}
// Heaviside function centered at 0.5*r
double hyperbolicTangent2(const mfem::Vector &x) {
  auto a = 2.;
  const auto epsilon = 2.e-4;
  const auto L = 2.e-3;  // 1.e-3;
  const auto r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
  return 0.5 + 0.5 * std::tanh(a * (r - 0.5 * L) / epsilon);
}

// Analytical solution for -Delta u = f(x,y) in [0,1]x[0,1]
// with f(x,y) = 6x*y*(1-y)-2x**3
// solution :  = y*(1-y)*x**3
double exactPoissonFunctionT(const mfem::Vector &x) {
  if (x[0] < 1.) {
    return 0;

  } else {
    return x[1] * (1. - x[1]);
  }
}
// Analytical solution for -Delta u = f(x,y) in [0,1]x[0,1]
// with f(x,y) = 6x*y*(1-y)-2x**3
// solution :  = y*(1-y)*x**3
double exactPoissonFunctionA(const mfem::Vector &x) {
  return x[1] * (1. - x[1]) * std::pow(x[0], 3.);
}
// RHS for -Delta u = f(x,y) in [0,1]x[0,1]
// RHS: = 6x*y*(1-y)-2x**3
double PoissonFunctionA(const mfem::Vector &x) {
  return 6. * x[0] * x[1] * (1. - x[1]) - 2. * std::pow(x[0], 3.);
}
// Analytical solution for -Delta u = f(x,y) in [0,1]x[0,1]
// with f(x,y) = 1.
double exactPoissonFunctionB(const mfem::Vector &x) {
  auto exa = 1.;
  return exa;
}

}  // namespace analytical
