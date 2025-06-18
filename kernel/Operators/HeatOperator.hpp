/**
 * @file HeatOperator.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Heat operator (Base, Steady and TimeDependent)
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
#include "Operators/PhaseFieldOperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Operators/SteadyPhaseFieldOperatorBase.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Base class for Mass Heat
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
class HeatOperatorBase : public OPEBASE<T, DIM, NLFI, LHS_NLFI> {
 protected:
  Parameters conductivity_params_;

 public:
  template <typename... Args>
  explicit HeatOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials, Args &&...args)
      : OPEBASE<T, DIM, NLFI, LHS_NLFI>(spatials, std::forward<Args>(args)...) {
    this->get_parameters();
  }
  template <typename... Args>
  HeatOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params,
                   Args &&...args)
      : OPEBASE<T, DIM, NLFI, LHS_NLFI>(spatials, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;
  NLFI *set_nlfi_ptr(const double dt, const std::vector<mfem::Vector> &u) override;
  void get_parameters() override;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const std::vector<mfem::Vector> &u) override;

  ~HeatOperatorBase();
};

/**
 * @brief Destroy the HeatOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::~HeatOperatorBase() {}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief  Set the NonLinearFormIntegrator dedicated to Heat
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param dt
 * @param u
 * @return NLFI*
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
NLFI *HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::set_nlfi_ptr(
    const double dt, const std::vector<mfem::Vector> &u) {
  Catch_Time_Section("HeatOperatorBase::set_nlfi_ptr");

  std::vector<mfem::ParGridFunction> vun;
  for (int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }

  const Parameters &all_params = this->conductivity_params_ + this->params_ - this->default_p_;
  NLFI *nlfi_ptr = new NLFI(vun, all_params, this->auxvariables_);
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
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
void HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description", "Heat operator");
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
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
void HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::ComputeEnergies(
    const int &it, const double &t, const double &dt, const std::vector<mfem::Vector> &u) {
  Catch_Time_Section("HeatOperatorBase::ComputeEnergies");

  std::vector<mfem::ParGridFunction> vun;
  vun.reserve(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(std::move(un));
  }

  std::vector<mfem::ParGridFunction *> vun_ptr;
  for (auto &un : vun) {
    vun_ptr.push_back(&un);
  }

  auto gphi = this->nlfi_ptr_->get_energy(vun_ptr, 1.);
  auto energy_coeff = gphi.get();
  mfem::ParGridFunction gf(this->fes_[0]);

  gf.ProjectCoefficient(*energy_coeff);

  // Calcul de l'intégrale de l'objet FunctionCoefficient sur le domaine
  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeLpError(1, zero);
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", energy));
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class SteadyHeatOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI = mfem::BlockNonlinearFormIntegrator>
class SteadyHeatOperator final
    : public HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyPhaseFieldOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyHeatOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials, Args &&...args)
      : HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyPhaseFieldOperatorBase>(
            spatials, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyHeatOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                     const Parameters &params, Args &&...args)
      : HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyPhaseFieldOperatorBase>(
            spatials, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  void overload_conductivity(const Parameters &p_params);

  ~SteadyHeatOperator() {}
};

/**
 * @brief Set the default options for properties
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyHeatOperator<T, DIM, NLFI, LHS_NLFI>::set_default_properties() {
  this->conductivity_params_ = Parameters(Parameter("lambda", 1.0));
}

/**
 * @brief Overload the default options for conductivity
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param COND
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyHeatOperator<T, DIM, NLFI, LHS_NLFI>::overload_conductivity(const Parameters &p_params) {
  this->conductivity_params_ = p_params;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class HeatOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI, Density RHO, HeatCapacity CP>
class HeatOperator final : public HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, PhaseFieldOperatorBase> {
 protected:
  Parameters density_params_;
  Parameters heat_capacity_params_;
  void set_default_properties() override;
  void get_mass_coefficient(const mfem::Vector &u) override;

 public:
  template <typename... Args>
  HeatOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials, Args &&...args)
      : HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, PhaseFieldOperatorBase>(
            spatials, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  HeatOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params,
               Args &&...args)
      : HeatOperatorBase<T, DIM, NLFI, LHS_NLFI, PhaseFieldOperatorBase>(
            spatials, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }

  void overload_density(const Parameters &p_params);
  void overload_heat_capacity(const Parameters &p_params);
  void overload_conductivity(const Parameters &p_params);

  ~HeatOperator() {}
};

/**
 * @brief Overload the MassMatrix coefficient definition
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 * @param u
 */
template <class T, int DIM, class NLFI, class LHS_NLFI, Density RHO, HeatCapacity CP>
void HeatOperator<T, DIM, NLFI, LHS_NLFI, RHO, CP>::get_mass_coefficient(const mfem::Vector &u) {
  if (this->MassCoeff_ != nullptr) {
    delete this->MassCoeff_;
  }
  if (this->mass_gf_ != nullptr) {
    delete this->mass_gf_;
  }
  this->mass_gf_ = new mfem::ParGridFunction(this->fes_[0]);
  this->mass_gf_->SetFromTrueDofs(u);
  auto density_coefficient = new DensityCoefficient<0, RHO>(this->mass_gf_, this->density_params_);
  auto cp_coefficient =
      new HeatCapacityCoefficient<0, CP>(this->mass_gf_, this->heat_capacity_params_);
  this->MassCoeff_ = new mfem::ProductCoefficient(*density_coefficient, *cp_coefficient);
}

/**
 * @brief Set the default options for properties
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI, Density RHO, HeatCapacity CP>
void HeatOperator<T, DIM, NLFI, LHS_NLFI, RHO, CP>::set_default_properties() {
  this->density_params_ = Parameters(Parameter("rho", 1.0));

  this->heat_capacity_params_ = Parameters(Parameter("cp", 1.0));
}
/**
 * @brief Overload the default options for density
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI, Density RHO, HeatCapacity CP>
void HeatOperator<T, DIM, NLFI, LHS_NLFI, RHO, CP>::overload_density(const Parameters &p_params) {
  this->density_params_ = p_params;
}

/**
 * @brief Overload the default options for heat capacity
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI, Density RHO, HeatCapacity CP>
void HeatOperator<T, DIM, NLFI, LHS_NLFI, RHO, CP>::overload_heat_capacity(
    const Parameters &p_params) {
  this->heat_capacity_params_ = p_params;
}

/**
 * @brief Overload the default options for conductivity
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI, Density RHO, HeatCapacity CP>
void HeatOperator<T, DIM, NLFI, LHS_NLFI, RHO, CP>::overload_conductivity(
    const Parameters &p_params) {
  this->conductivity_params_ = p_params;
}
