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
#include <utility>
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
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief PhaseFieldOperator class
 *
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
class PhaseFieldOperator final : public OPEBASE<T, DIM, NLFI> {
 private:
  double omega_, lambda_;
  const Parameters &params_;

 public:
  template <typename... Args>
  PhaseFieldOperator(SpatialDiscretization<T, DIM> const *spatial, const Parameters &params,
                     Args &&...args)
      : OPEBASE<T, DIM, NLFI>(spatial, params, std::forward<Args>(args)...), params_(params) {
    this->get_parameters(params);
  }

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
 * @brief Destroy the PhaseFieldOperator< T, DIM, NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
PhaseFieldOperator<T, DIM, NLFI, OPEBASE>::~PhaseFieldOperator() {}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Set the NonLinearFormIntegrator dedicated to AllenCahn
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
NLFI *PhaseFieldOperator<T, DIM, NLFI, OPEBASE>::set_nlfi_ptr(const double dt,
                                                              const mfem::Vector &u) {
  Catch_Time_Section("PhaseFieldOperator::set_nlfi_ptr");

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
void PhaseFieldOperator<T, DIM, NLFI, OPEBASE>::get_parameters(const Parameters &params) {
  this->omega_ = params.get_param_value<double>("omega");
  this->lambda_ = params.get_param_value<double>("lambda");
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
void PhaseFieldOperator<T, DIM, NLFI, OPEBASE>::ComputeEnergies(const int &it, const double &dt,
                                                                const double &t,
                                                                const mfem::Vector &u) {
  Catch_Time_Section("PhaseFieldOperator::ComputeEnergies");

  mfem::ParGridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::ParGridFunction gf(this->fespace_);
  EnergyCoefficient g(&un_gf, 0.5 * this->lambda_, this->omega_);
  gf.ProjectCoefficient(g);
  mfem::ParGridFunction sigf(this->fespace_);
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
