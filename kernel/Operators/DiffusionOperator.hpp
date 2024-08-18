/**
 * @file DiffusionOperator.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Diffusion timedependent operator
 * @version 0.1
 * @date 2024-06-17
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Operators/PhaseFieldOperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Operators/SteadyPhaseFieldOperatorBase.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief DiffusionOperator class
 *
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
class DiffusionOperator final : public OPEBASE<T, DIM, NLFI> {
 private:
  const Parameters &params_;

 public:
  template <typename... Args>
  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Args &&...args)
      : OPEBASE<T, DIM, NLFI>(spatial, params, std::forward<Args>(args)...), params_(params) {
    this->get_parameters(params);
  }

  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;
  void get_parameters(const Parameters &vectr_param) override;
  void ComputeEnergies(const int &it, const double &dt, const double &t,
                       const mfem::Vector &u) override;
  ~DiffusionOperator();
};

/**
 * @brief Destroy the Diffusion Operator< T,  DIM,  NLFI,  OPEBASE>:: Diffusion Operator
 * object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
DiffusionOperator<T, DIM, NLFI, OPEBASE>::~DiffusionOperator() {}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief  Set the NonLinearFormIntegrator dedicated to diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
NLFI *DiffusionOperator<T, DIM, NLFI, OPEBASE>::set_nlfi_ptr(const double dt,
                                                             const mfem::Vector &u) {
  Catch_Time_Section("DiffusionOperator::set_nlfi_ptr");

  mfem::ParGridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  NLFI *nlfi_ptr = new NLFI(un, this->params_);
  return nlfi_ptr;
}

/**
 * @brief Get parameters for use in the current operator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param params
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
void DiffusionOperator<T, DIM, NLFI, OPEBASE>::get_parameters(const Parameters &params) {
  // this->alpha_ = params.get_param_value<double>("alpha");
  // this->kappa_ = params.get_param_value<double>("kappa");
}

/**
 * @brief  Compute Energy contribution at the end of time step
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param it
 * @param dt
 * @param t
 * @param u
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
void DiffusionOperator<T, DIM, NLFI, OPEBASE>::ComputeEnergies(const int &it, const double &dt,
                                                               const double &t,
                                                               const mfem::Vector &u) {
  Catch_Time_Section("DiffusionOperator::ComputeEnergies");
  mfem::ParGridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::ParGridFunction gf(this->fespace_);
  EnergyCoefficient g(&un_gf, 0., 1.);
  gf.ProjectCoefficient(g);

  // Calcul de l'intégrale de l'objet FunctionCoefficient sur le domaine
  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeL1Error(zero);
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", energy));
}
