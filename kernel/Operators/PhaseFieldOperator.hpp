/*
 * Copyright © CEA 2022
 *
 * \brief PhaseField time-dependent operator built on the basis of ex16.cpp from MFEM repository
 *
 *   After spatial discretization, the phasefield model can be written as:
 *      dphi/dt = M^{-1}(-Kphi)
 *   where phi denotes the phase-field variable and K the phase-field operator (that does not
 *   depends on time)
 *
 * \file PhaseFieldOperator.hpp
 * \author ci230846
 * \date 12/01/2022
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
#include "Profiling/Profiling.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"

#pragma once

/**
 * @brief PhaseFieldOperator class
 *
 */
template <class T, int DIM, class NLFI>
class PhaseFieldOperator final : public PhaseFieldOperatorBase<T, DIM, NLFI> {
 private:
  double mobility_coeff_, omega_, lambda_;

 public:
  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars);
  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars, Variables<T, DIM> &auxvars);

  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars, AnalyticalFunctions<DIM> source_term_name);

  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars, Variables<T, DIM> &auxvars,
                     AnalyticalFunctions<DIM> source_term_name);

  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;
  void get_parameters(const Parameters &vectr_param) override;
  void ComputeEnergies(const int &it, const double &dt, const double &t,
                       const mfem::Vector &u) override;

  ~PhaseFieldOperator();
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
PhaseFieldOperator<T, DIM, NLFI>::PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial,
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
PhaseFieldOperator<T, DIM, NLFI>::PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial,
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
PhaseFieldOperator<T, DIM, NLFI>::PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial,
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
PhaseFieldOperator<T, DIM, NLFI>::PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial,
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
PhaseFieldOperator<T, DIM, NLFI>::~PhaseFieldOperator() {}

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
NLFI *PhaseFieldOperator<T, DIM, NLFI>::set_nlfi_ptr(const double dt, const mfem::Vector &u) {
  Catch_Time_Section("PhaseFieldOperator::set_nlfi_ptr");

  mfem::GridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);
  NLFI *nlfi_ptr = new NLFI(un, this->omega_, this->lambda_, this->mobility_coeff_);

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
void PhaseFieldOperator<T, DIM, NLFI>::get_parameters(const Parameters &params) {
  this->omega_ = params.get_param_value<double>("omega");
  this->lambda_ = params.get_param_value<double>("lambda");
  this->mobility_coeff_ = params.get_param_value<double>("mobility");
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
void PhaseFieldOperator<T, DIM, NLFI>::ComputeEnergies(const int &it, const double &dt,
                                                       const double &t, const mfem::Vector &u) {
  Catch_Time_Section("PhaseFieldOperator::ComputeEnergies");

  mfem::GridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::GridFunction gf(this->fespace_);
  EnergyCoefficient g(&un_gf, 0.5 * this->lambda_, this->omega_);
  gf.ProjectCoefficient(g);
  mfem::GridFunction sigf(this->fespace_);
  EnergyCoefficient sig(&un_gf, this->lambda_, 0.);
  sigf.ProjectCoefficient(sig);

  // Calcul de l'intégrale de l'objet FunctionCoefficient sur le domaine
  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeL1Error(zero);
  const auto interfacial_energy = sigf.ComputeL1Error(zero);

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", energy));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Sigma[J.m-3]", interfacial_energy));
}
