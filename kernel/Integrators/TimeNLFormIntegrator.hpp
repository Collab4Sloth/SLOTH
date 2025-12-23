/**
 * @file TimeNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative
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
#include <string>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the Allen Cahn equation
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS>
class TimeNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  void check_variables_consistency();

 protected:
  std::vector<mfem::ParGridFunction> temp_gf_;
  std::list<GlossaryType> expected_list_;
  Coefficients coefficient_A;
  Coefficients coefficient_B;
  void get_coefficients() override;
  virtual double compute_coefficient(Coefficient coef, const std::vector<double>& values);

 public:
  void init() override;
  TimeNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old, const Parameters& params,
                       std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmats) override;
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new TimeNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 */
template <class VARS>
TimeNLFormIntegrator<VARS>::TimeNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                                                 const Parameters& params,
                                                 std::vector<VARS*> auxvars,
                                                 const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {
  this->check_variables_consistency();
}
template <class VARS>
void TimeNLFormIntegrator<VARS>::init() {
  if (this->expected_list_.empty()) {
    this->get_coefficients();
  } else {
    this->check_coefficient_types(this->expected_list_);
    this->get_coefficients();
  }
}

/**
 * @brief Check variables consistency
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "TimeNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      temperature_found = true;
      break;
    }
  }
}

/**
 * @brief Get default coefficients
 * @remark Could be overridden by child classes
 *
 * @tparam VARS
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    this->coefficient_A.add(Coefficient(Glossary::Default, 1.0));
    this->coefficient_B.add(Coefficient(Glossary::Default, 1.0));
  }
}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Catch_Time_Section("TimeNLFormIntegrator:AssembleElementVector");
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    // elvect.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      const auto& u = *elfun[blk] * Psi;
      const auto& un = this->u_old_[blk].GetValue(Tr, ip);

      // el[blk]->CalcPhysDShape(Tr, gradPsi);
      // gradPsi.MultTranspose(*elfun[blk], gradU);
      double coef_a = this->compute_coefficient(this->coefficient_A[blk], {u, un});
      double coef_b = this->compute_coefficient(this->coefficient_B[blk], {u, un});
      const double ww = coef_a * coef_b * u * ip.weight * Tr.Weight();
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("TimeNLFormIntegrator::AssembleElementGrad");
  int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    // int nd = el.GetDof();
    // int dim = el.GetDim();
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    // elmat.SetSize(nd);
    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    // elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      Tr.SetIntPoint(&ip);
      const auto& u = *elfun[blk] * Psi;
      const auto& un = this->u_old_[blk].GetValue(Tr, ip);

      double coef_a = this->compute_coefficient(this->coefficient_A[blk], {u, un});
      double coef_b = this->compute_coefficient(this->coefficient_B[blk], {u, un});
      double fun_val = coef_a * coef_b * ip.weight * Tr.Weight();
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));
    }
  }

  // Off-diagonal
  for (int blk = 0; blk < num_blocks; ++blk) {
    for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
      if (off_blk != blk) {
        int nd = el[blk]->GetDof();

        // elmat.SetSize(nd);
        elmats(blk, off_blk)->SetSize(nd);
        *elmats(blk, off_blk) = 0.0;
      }
    }
  }
}

/**
 * @brief  Return the value of a coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 * @tparam VARS
 * @param coef
 * @param values
 * @return double
 */
template <class VARS>
double TimeNLFormIntegrator<VARS>::compute_coefficient(Coefficient coef,
                                                       const std::vector<double>& values) {
  const double u = values[0];
  const double un = values[1];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    // coef_value = coef.compute({u});
    MFEM_VERIFY(false,
                "Implicit coefficient for TimeDerivative integrator not implemented yet. Please "
                "check your data.");
  } else if (coef.is_explicit()) {
    coef_value = coef.compute({un});
  } else if (coef.is_scalar()) {
    coef_value = coef.compute();
  }
  return coef_value;
}