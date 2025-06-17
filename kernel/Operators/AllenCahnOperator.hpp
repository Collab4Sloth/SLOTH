/**
 * @file AllenCahnOperator.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief AllenCahn operator (Base, Steady and TimeDependent)
 * @version 0.1
 * @date 2024-09-04
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
#include "Operators/SteadyReducedOperator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Base class for AllenCahn
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
class AllenCahnOperatorBase : public OPEBASE<T, DIM, NLFI, LHS_NLFI> {
 protected:
  Parameters mobility_params_;
  double omega_, lambda_;

 public:
  template <typename... Args>
  explicit AllenCahnOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                                 Args &&...args)
      : OPEBASE<T, DIM, NLFI, LHS_NLFI>(spatials, std::forward<Args>(args)...) {
    this->get_parameters();
  }
  template <typename... Args>
  AllenCahnOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                        const Parameters &params, Args &&...args)
      : OPEBASE<T, DIM, NLFI, LHS_NLFI>(spatials, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;
  NLFI *set_nlfi_ptr(const double dt, const std::vector<mfem::Vector> &u) override;
  void get_parameters() override;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const std::vector<mfem::Vector> &u) override;

  ~AllenCahnOperatorBase();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Destroy the AllenCahnOperatorBase< T, DIM, NLFI>:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::~AllenCahnOperatorBase() {}

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
template <class T, int DIM, class NLFI, class LHS_NLFI,
          template <class, int, class, class> class OPEBASE>
NLFI *AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::set_nlfi_ptr(
    const double dt, const std::vector<mfem::Vector> &u) {
  Catch_Time_Section("AllenCahnOperatorBase::set_nlfi_ptr");
  std::vector<mfem::ParGridFunction> vun;
  for (int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }

  const Parameters &all_params = this->mobility_params_ + this->params_ - this->default_p_;

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
void AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::get_parameters() {
  this->omega_ = this->params_.template get_param_value<double>("omega");
  this->lambda_ = this->params_.template get_param_value<double>("lambda");
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
void AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, OPEBASE>::ComputeEnergies(
    const int &it, const double &t, const double &dt, const std::vector<mfem::Vector> &u) {
  Catch_Time_Section("AllenCahnOperatorBase::ComputeEnergies");

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

  auto gphi = this->nlfi_ptr_->get_energy(vun_ptr, this->omega_);

  std::unique_ptr<GradientEnergyCoefficient> gradGphi =
      this->nlfi_ptr_->get_grad_energy(vun_ptr, 0.5 * this->lambda_);
  auto energy_coeff = gphi.get();

  mfem::ParGridFunction gf(this->fes_[0]);
  gf.ProjectCoefficient(*energy_coeff);

  mfem::ParGridFunction sigf(this->fes_[0]);
  auto grad_energy_coeff = gradGphi.get();
  sigf.ProjectCoefficient(*grad_energy_coeff);

  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeLpError(1, zero);
  const auto interfacial_energy = sigf.ComputeLpError(1, zero);

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Density[J.m-3]", energy + interfacial_energy));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Sigma[J.m-3]", 2. * interfacial_energy));
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class SteadyAllenCahnOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI = mfem::BlockNonlinearFormIntegrator>
class SteadyAllenCahnOperator final
    : public AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyPhaseFieldOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyAllenCahnOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials, Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyPhaseFieldOperatorBase>(
            spatials, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyAllenCahnOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                          const Parameters &params, Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, SteadyPhaseFieldOperatorBase>(
            spatials, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  void overload_mobility(const Parameters &p_params);

  ~SteadyAllenCahnOperator() {}
};

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyAllenCahnOperator<T, DIM, NLFI, LHS_NLFI>::set_default_properties() {
  this->mobility_params_ = Parameters(Parameter("mob", 1.0));
}
/**
 * @brief Overload the default options for mobility
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyAllenCahnOperator<T, DIM, NLFI, LHS_NLFI>::overload_mobility(
    const Parameters &p_params) {
  this->mobility_params_ = p_params;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class AllenCahnOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
class AllenCahnOperator final
    : public AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, PhaseFieldOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  AllenCahnOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials, Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, PhaseFieldOperatorBase>(
            spatials, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  AllenCahnOperator(std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params,
                    Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, LHS_NLFI, PhaseFieldOperatorBase>(
            spatials, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }

  void overload_mobility(const Parameters &p_params);
  void get_mass_coefficient(const mfem::Vector &u) override;

  ~AllenCahnOperator() {}
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void AllenCahnOperator<T, DIM, NLFI, LHS_NLFI>::get_mass_coefficient(const mfem::Vector &u) {
  PhaseFieldOperatorBase<T, DIM, NLFI, LHS_NLFI>::get_mass_coefficient(u);
}

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void AllenCahnOperator<T, DIM, NLFI, LHS_NLFI>::set_default_properties() {
  this->mobility_params_ = Parameters(Parameter("mob", 1.0));
}
/**
 * @brief Overload the default options for mobility
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void AllenCahnOperator<T, DIM, NLFI, LHS_NLFI>::overload_mobility(const Parameters &p_params) {
  this->mobility_params_ = p_params;
}
