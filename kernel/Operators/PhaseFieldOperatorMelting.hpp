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
 * \file PhaseFieldOperatorMelting.hpp
 * \author ci230846
 * \date 12/01/2022
 */

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <vector>
#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Integrators/AllenCahnMeltingNLFormIntegrator.hpp"
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
 * @brief PhaseFieldOperatorMelting class
 *
 */
template <class T, int DIM, class NLFI>
class PhaseFieldOperatorMelting final : public PhaseFieldOperatorBase<T, DIM, NLFI> {
 public:
  PhaseFieldOperatorMelting(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                            Variables<T, DIM> &vars);

  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;

  virtual ~PhaseFieldOperatorMelting();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/**
 * @brief Construct a new Phase Field Operator:: Phase Field Operator object
 *
 * @param fespace Finite Element space
 * @param params list of Parameters
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorMelting<T, DIM, NLFI>::PhaseFieldOperatorMelting(
    SpatialDiscretization<T, DIM> *spatial, const Parameters &params, Variables<T, DIM> &vars)
    : PhaseFieldOperatorBase<T, DIM, NLFI>(spatial, params, vars) {}

/**
 * @brief Set the NonLinearFormIntegrator dedicated to melting coupling
 *
 * @tparam T
 * @tparam DIM
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM, class NLFI>
NLFI *PhaseFieldOperatorMelting<T, DIM, NLFI>::set_nlfi_ptr(const double dt,
                                                            const mfem::Vector &u) {
  mfem::GridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  NLFI *nlfi_ptr =
      new NLFI(un, this->omega_, this->lambda_, this->mobility_coeff_, this->phase_change_coeff_);
  return nlfi_ptr;
}

/**
 * @brief Destroy the Phase Field Operator Melting< T,  DIM,  NLFI>:: Phase Field Operator Melting
 * object
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorMelting<T, DIM, NLFI>::~PhaseFieldOperatorMelting() {}
