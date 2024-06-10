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
#include "Integrators/DiffusionNLFIntegrator.hpp"
#include "Operators/PhaseFieldOperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/UtilsForSolvers.hpp"
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
class PhaseFieldOperator : public PhaseFieldOperatorBase<T, DIM, NLFI> {
 public:
  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars, Variables<T, DIM> &auxvars);

  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars);

  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars, AnalyticalFunctions<DIM> source_term_name);

  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;

  virtual ~PhaseFieldOperator();
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
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars, auxvars) {}

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
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars) {}

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
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars, source_term_name) {}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Set the NonLinearFormIntegrator dedicated to standard phasefield calculations
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
  //------Start profiling-------------------------
  Timers timer_set_nlfi_ptr("set_nlfi_ptr");
  timer_set_nlfi_ptr.start();
  //------------------------------------------------

  mfem::GridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);
  NLFI *nlfi_ptr = new NLFI(un, this->omega_, this->lambda_, this->mobility_coeff_);

  // save the results of profiling 
  timer_set_nlfi_ptr.stop();
  UtilsForOutput::getInstance().update_timer("set_nlfi_ptr",timer_set_nlfi_ptr );
  //------------------------------------------------
  return nlfi_ptr;
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
