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
#include "Coefficients/PhaseChangeCoefficient.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Integrators/DiffusionNLFIntegrator.hpp"
#include "Operators/PhaseFieldOperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Utils/UtilsForSolvers.hpp"
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
                     Variables<T, DIM> &vars);

  template <class... Args>
  PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                     Variables<T, DIM> &vars, const std::string &source_term_name, Args... args);

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
template <class... Args>
PhaseFieldOperator<T, DIM, NLFI>::PhaseFieldOperator(SpatialDiscretization<T, DIM> *spatial,
                                                     const Parameters &params,
                                                     Variables<T, DIM> &vars,
                                                     const std::string &source_term_name,
                                                     Args... args)
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars, source_term_name, args...) {}

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
  mfem::GridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);
  NLFI *nlfi_ptr = new NLFI(un, this->omega_, this->lambda_, this->mobility_coeff_);
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
