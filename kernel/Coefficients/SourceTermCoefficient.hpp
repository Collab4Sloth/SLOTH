/*
 * Copyright © CEA 2023
 *
 * SourceTermCoefficient.hpp
 *
 *  Created on: 15 may 2023
 *      Author: ci230846
 */
#include <string>

#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

class SourceTermCoefficient : public mfem::Coefficient {
 private:
  std::string source_term_name_;

 public:
  explicit SourceTermCoefficient(const std::string &source_term_name);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
  ~SourceTermCoefficient();
};

/**
 * @brief Construct a new Source Term Coefficient< SRC>:: Source Term Coefficient object
 *
 */
SourceTermCoefficient::SourceTermCoefficient(const std::string &source_term_name)
    : source_term_name_(source_term_name) {}

/**
 * @brief Evaluation of the source term coefficient
 *
 * @param ip
 * @return double
 */
double SourceTermCoefficient::Eval(mfem::ElementTransformation &T,
                                   const mfem::IntegrationPoint &ip) {
  switch (SourceTerm::from(this->source_term_name_)) {
    case SourceTerm::Null: {
      const auto src_value = 0.;
      return src_value;
    }

    case SourceTerm::Sinusoide2D: {
      const auto time = this->GetTime();
      const auto u = std::exp(-2. * time) * std::sin(ip.x + ip.y);
      const auto src_value = u * u * u - u;
      return -1. * src_value / (0.09);
    }
    default:
      mfem::mfem_error(
          "SourceTermCoefficient::Eval: only Null, Sinusoide2D source term are available");
      break;
  }
}

/**
 * @brief Destroy the Source Term Coefficient< SRC>:: Source Term Coefficient object
 *
 * @tparam SRC
 */
SourceTermCoefficient::~SourceTermCoefficient() {}
