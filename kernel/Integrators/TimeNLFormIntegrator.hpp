/**
 * @file TimeNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief VF of the time derivative
 * @version 0.1
 * @date 2025-06-04
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/LambdaCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/OmegaCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
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
class TimeNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator,
                             public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  void check_variables_consistency();

 protected:
  mfem::ParGridFunction u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  std::vector<mfem::Vector> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_gf_infos_;
  std::vector<mfem::ParGridFunction> temp_gf_;

 public:
  TimeNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                       std::vector<VARS*> auxvars);
  ~TimeNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvec);

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmats);
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
TimeNLFormIntegrator<VARS>::TimeNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                                 const Parameters& params,
                                                 std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->check_variables_consistency();
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
  this->aux_gf_ = this->get_aux_gf();
  this->aux_old_gf_ = this->get_aux_old_gf();
  this->aux_gf_infos_ = this->get_aux_infos();

  // Temperature scaling for mobility
  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
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
  int blk = 0;
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

    el[blk]->CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(*elfun[blk], gradU);

    const double ww = u * ip.weight * Tr.Weight();
    add(*elvect[blk], ww, Psi, *elvect[blk]);
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
  int blk = 0;
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
    double fun_val = ip.weight * Tr.Weight();
    AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));
  }
}

/**
 * @brief Destroy the TimeNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS>
TimeNLFormIntegrator<VARS>::~TimeNLFormIntegrator() {}
