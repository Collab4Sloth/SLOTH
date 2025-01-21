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
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
class AllenCahnOperatorBase : public OPEBASE<T, DIM, NLFI> {
 protected:
  Parameters mobility_params_;
  double omega_, lambda_;

 public:
  template <typename... Args>
  explicit AllenCahnOperatorBase(SpatialDiscretization<T, DIM> const *spatial, Args &&...args)
      : OPEBASE<T, DIM, NLFI>(spatial, std::forward<Args>(args)...) {
    this->get_parameters();
  }
  template <typename... Args>
  AllenCahnOperatorBase(SpatialDiscretization<T, DIM> const *spatial, const Parameters &params,
                        Args &&...args)
      : OPEBASE<T, DIM, NLFI>(spatial, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;
  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;
  void get_parameters() override;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const mfem::Vector &u) override;

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
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
AllenCahnOperatorBase<T, DIM, NLFI, OPEBASE>::~AllenCahnOperatorBase() {}

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
NLFI *AllenCahnOperatorBase<T, DIM, NLFI, OPEBASE>::set_nlfi_ptr(const double dt,
                                                                 const mfem::Vector &u) {
  Catch_Time_Section("AllenCahnOperatorBase::set_nlfi_ptr");

  mfem::ParGridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  std::vector<mfem::ParGridFunction> aux_gf = this->get_auxiliary_gf();

  const Parameters &all_params = this->mobility_params_ + this->params_ - this->default_p_;

  NLFI *nlfi_ptr = new NLFI(un, all_params, aux_gf);

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
void AllenCahnOperatorBase<T, DIM, NLFI, OPEBASE>::get_parameters() {
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
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
void AllenCahnOperatorBase<T, DIM, NLFI, OPEBASE>::ComputeEnergies(const int &it, const double &t,
                                                                   const double &dt,
                                                                   const mfem::Vector &u) {
  Catch_Time_Section("AllenCahnOperatorBase::ComputeEnergies");

  mfem::ParGridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::ParGridFunction gf(this->fespace_);
  auto gphi = this->nlfi_ptr_->get_energy(&un_gf, this->omega_);
  auto gradGphi = this->nlfi_ptr_->get_grad_energy(&un_gf, 0.5 * this->lambda_);
  auto energy_coeff = gphi.get();
  gf.ProjectCoefficient(*energy_coeff);

  mfem::ParGridFunction sigf(this->fespace_);
  auto grad_energy_coeff = gradGphi.get();
  sigf.ProjectCoefficient(*grad_energy_coeff);

  // Calcul de l'intégrale de l'objet FunctionCoefficient sur le domaine
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
template <class T, int DIM, class NLFI>
class SteadyAllenCahnOperator final
    : public AllenCahnOperatorBase<T, DIM, NLFI, SteadyPhaseFieldOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyAllenCahnOperator(SpatialDiscretization<T, DIM> *spatial, Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, SteadyPhaseFieldOperatorBase>(
            spatial, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyAllenCahnOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                          Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, SteadyPhaseFieldOperatorBase>(
            spatial, params, std::forward<Args>(args)...) {
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
template <class T, int DIM, class NLFI>
void SteadyAllenCahnOperator<T, DIM, NLFI>::set_default_properties() {
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
template <class T, int DIM, class NLFI>
void SteadyAllenCahnOperator<T, DIM, NLFI>::overload_mobility(const Parameters &p_params) {
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
template <class T, int DIM, class NLFI>
class AllenCahnOperator final : public AllenCahnOperatorBase<T, DIM, NLFI, PhaseFieldOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  AllenCahnOperator(SpatialDiscretization<T, DIM> *spatial, Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, PhaseFieldOperatorBase>(spatial,
                                                                    std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  AllenCahnOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Args &&...args)
      : AllenCahnOperatorBase<T, DIM, NLFI, PhaseFieldOperatorBase>(spatial, params,
                                                                    std::forward<Args>(args)...) {
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
template <class T, int DIM, class NLFI>
void AllenCahnOperator<T, DIM, NLFI>::get_mass_coefficient(const mfem::Vector &u) {
  PhaseFieldOperatorBase<T, DIM, NLFI>::get_mass_coefficient(u);
}

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
void AllenCahnOperator<T, DIM, NLFI>::set_default_properties() {
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
template <class T, int DIM, class NLFI>
void AllenCahnOperator<T, DIM, NLFI>::overload_mobility(const Parameters &p_params) {
  this->mobility_params_ = p_params;
}
