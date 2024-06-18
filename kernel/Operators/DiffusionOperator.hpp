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
#include <vector>

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Operators/PhaseFieldOperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"

#pragma once

/**
 * @brief DiffusionOperator class
 *
 */
template <class T, int DIM, class NLFI>
class DiffusionOperator final : public PhaseFieldOperatorBase<T, DIM, NLFI> {
 private:
  double alpha_, kappa_;

 public:
  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Variables<T, DIM> &vars);
  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Variables<T, DIM> &vars, Variables<T, DIM> &auxvars);

  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Variables<T, DIM> &vars, AnalyticalFunctions<DIM> source_term_name);

  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Variables<T, DIM> &vars, Variables<T, DIM> &auxvars,
                    AnalyticalFunctions<DIM> source_term_name);

  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;
  void get_parameters(const Parameters &vectr_param) override;
  void ComputeEnergies(const int &it, const double &dt, const double &t,
                       const mfem::Vector &u) override;
  ~DiffusionOperator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Phase Field Operator< T,  DIM,  NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 * @param auxvars
 */
template <class T, int DIM, class NLFI>
DiffusionOperator<T, DIM, NLFI>::DiffusionOperator(SpatialDiscretization<T, DIM> *spatial,
                                                   const Parameters &params,
                                                   Variables<T, DIM> &vars,
                                                   Variables<T, DIM> &auxvars)
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars, auxvars) {
  this->get_parameters(params);
}

/**
 * @brief Construct a new Phase Field Operator< T,  DIM,  NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 */
template <class T, int DIM, class NLFI>
DiffusionOperator<T, DIM, NLFI>::DiffusionOperator(SpatialDiscretization<T, DIM> *spatial,
                                                   const Parameters &params,
                                                   Variables<T, DIM> &vars)
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars) {
  this->get_parameters(params);
}

/**
 * @brief Construct a new Phase Field Operator< T,  DIM,  NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 */
template <class T, int DIM, class NLFI>
DiffusionOperator<T, DIM, NLFI>::DiffusionOperator(SpatialDiscretization<T, DIM> *spatial,
                                                   const Parameters &params,
                                                   Variables<T, DIM> &vars,
                                                   AnalyticalFunctions<DIM> source_term_name)
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars, source_term_name) {
  this->get_parameters(params);
}

/**
 * @brief Construct a new Phase Field Operator<T, DIM, NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 * @param auxvars
 * @param source_term_name
 */
template <class T, int DIM, class NLFI>
DiffusionOperator<T, DIM, NLFI>::DiffusionOperator(SpatialDiscretization<T, DIM> *spatial,
                                                   const Parameters &params,
                                                   Variables<T, DIM> &vars,
                                                   Variables<T, DIM> &auxvars,
                                                   AnalyticalFunctions<DIM> source_term_name)
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars, auxvars, source_term_name) {
  this->get_parameters(params);
}

/**
 * @brief Destroy the Phase Field Operator Melting< T,  DIM,  NLFI>:: Phase Field Operator
 * object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
DiffusionOperator<T, DIM, NLFI>::~DiffusionOperator() {}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Set the NonLinearFormIntegrator dedicated to standard AllenCahn calculation
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM, class NLFI>
NLFI *DiffusionOperator<T, DIM, NLFI>::set_nlfi_ptr(const double dt, const mfem::Vector &u) {
  Catch_Time_Section("DiffusionOperator::set_nlfi_ptr");

  mfem::GridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  NLFI *nlfi_ptr = new NLFI(un, this->alpha_, this->kappa_);
  return nlfi_ptr;
}

/**
 * @brief Get parameters values for the given Operator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param params
 */
template <class T, int DIM, class NLFI>
void DiffusionOperator<T, DIM, NLFI>::get_parameters(const Parameters &params) {
  this->alpha_ = params.get_param_value<double>("alpha");
  this->kappa_ = params.get_param_value<double>("kappa");
}

/**
 * @brief Compute Energy contribution at the end of time step
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param iter
 * @param u
 */
template <class T, int DIM, class NLFI>
void DiffusionOperator<T, DIM, NLFI>::ComputeEnergies(const int &it, const double &dt,
                                                      const double &t, const mfem::Vector &u) {
  Catch_Time_Section("DiffusionOperator::ComputeEnergies");
  mfem::GridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::GridFunction gf(this->fespace_);
  EnergyCoefficient g(&un_gf, 0., 1.);
  gf.ProjectCoefficient(g);

  // Calcul de l'intégrale de l'objet FunctionCoefficient sur le domaine
  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeL1Error(zero);
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", energy));
}
