/**
 * @file DiffusionOperator.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Diffusion operator (Base, Steady and TimeDependent)
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
#include "Coefficients/DensityCoefficient.hpp"
#include "Coefficients/DiffusionCoefficient.hpp"
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
 * @brief Base class for Mass Diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
class DiffusionOperatorBase : public OPEBASE<T, DIM, NLFI> {
 protected:
  Parameters diffusion_params_;

 public:
  template <typename... Args>
  explicit DiffusionOperatorBase(SpatialDiscretization<T, DIM> *spatial, Args &&...args)
      : OPEBASE<T, DIM, NLFI>(spatial, std::forward<Args>(args)...) {}

  template <typename... Args>
  DiffusionOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                        Args &&...args)
      : OPEBASE<T, DIM, NLFI>(spatial, params, std::forward<Args>(args)...) {
    this->get_parameters();
  }

  void set_default_properties() override = 0;

  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override;
  void get_parameters() override;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const mfem::Vector &u) override;

  ~DiffusionOperatorBase();
};

/**
 * @brief Destroy the DiffusionOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam OPEBASE
 */
template <class T, int DIM, class NLFI, template <class, int, class> class OPEBASE>
DiffusionOperatorBase<T, DIM, NLFI, OPEBASE>::~DiffusionOperatorBase() {}

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
NLFI *DiffusionOperatorBase<T, DIM, NLFI, OPEBASE>::set_nlfi_ptr(const double dt,
                                                                 const mfem::Vector &u) {
  Catch_Time_Section("DiffusionOperatorBase::set_nlfi_ptr");

  mfem::ParGridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  std::vector<mfem::ParGridFunction> aux_gf = this->get_auxiliary_gf();
  const Parameters &all_params = this->diffusion_params_ + this->params_ - this->default_p_;
  NLFI *nlfi_ptr = new NLFI(un, all_params, this->auxvariables_);
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
void DiffusionOperatorBase<T, DIM, NLFI, OPEBASE>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description", "Diffusion operator");
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
void DiffusionOperatorBase<T, DIM, NLFI, OPEBASE>::ComputeEnergies(const int &it, const double &t,
                                                                   const double &dt,
                                                                   const mfem::Vector &u) {
  Catch_Time_Section("DiffusionOperatorBase::ComputeEnergies");
  mfem::ParGridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::ParGridFunction gf(this->fespace_);
  auto gphi = this->nlfi_ptr_->get_energy(&un_gf, 1.);
  auto energy_coeff = gphi.get();

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
 * @brief Class SteadyDiffusionOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
class SteadyDiffusionOperator final
    : public DiffusionOperatorBase<T, DIM, NLFI, SteadyPhaseFieldOperatorBase> {
 protected:
  void set_default_properties() override;

 public:
  template <typename... Args>
  SteadyDiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                          Args &&...args)
      : DiffusionOperatorBase<T, DIM, NLFI, SteadyPhaseFieldOperatorBase>(
            spatial, params, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  SteadyDiffusionOperator(SpatialDiscretization<T, DIM> *spatial, Args &&...args)
      : DiffusionOperatorBase<T, DIM, NLFI, SteadyPhaseFieldOperatorBase>(
            spatial, std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  void overload_diffusion(const Parameters &p_params);

  ~SteadyDiffusionOperator() {}
};
/**
 * @brief Set the default options for properties
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
void SteadyDiffusionOperator<T, DIM, NLFI>::set_default_properties() {
  this->diffusion_params_ = Parameters(Parameter("D", 1.0), Parameter("D_stab", 0.0));
}

/**
 * @brief Overload the default options for diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param DIFFU
 * @param p_params
 */
template <class T, int DIM, class NLFI>
void SteadyDiffusionOperator<T, DIM, NLFI>::overload_diffusion(const Parameters &p_params) {
  this->diffusion_params_ = p_params;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class DiffusionOperator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, Density DENS>
class DiffusionOperator final : public DiffusionOperatorBase<T, DIM, NLFI, PhaseFieldOperatorBase> {
 protected:
  Parameters density_params_;
  void set_default_properties() override;
  void get_mass_coefficient(const mfem::Vector &u) override;

 public:
  template <typename... Args>
  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                    Args &&...args)
      : DiffusionOperatorBase<T, DIM, NLFI, PhaseFieldOperatorBase>(spatial, params,
                                                                    std::forward<Args>(args)...) {
    this->set_default_properties();
  }
  template <typename... Args>
  DiffusionOperator(SpatialDiscretization<T, DIM> *spatial, Args &&...args)
      : DiffusionOperatorBase<T, DIM, NLFI, PhaseFieldOperatorBase>(spatial,
                                                                    std::forward<Args>(args)...) {
    this->set_default_properties();
  }

  void overload_density(const Parameters &p_params);
  void overload_diffusion(const Parameters &p_params);

  ~DiffusionOperator() {}
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
template <class T, int DIM, class NLFI, Density DENS>
void DiffusionOperator<T, DIM, NLFI, DENS>::get_mass_coefficient(const mfem::Vector &u) {
  if (this->MassCoeff_ != nullptr) {
    delete this->MassCoeff_;
  }
  if (this->mass_gf_ != nullptr) {
    delete this->mass_gf_;
  }
  this->mass_gf_ = new mfem::ParGridFunction(this->fespace_);
  this->mass_gf_->SetFromTrueDofs(u);
  auto density_coefficient = new DensityCoefficient<0, DENS>(this->mass_gf_, this->density_params_);
  auto one = new mfem::ConstantCoefficient(1.0);
  this->MassCoeff_ = new mfem::ProductCoefficient(*density_coefficient, *one);
}

/**
 * @brief Set the default options for properties

 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, Density DENS>
void DiffusionOperator<T, DIM, NLFI, DENS>::set_default_properties() {
  this->density_params_ = Parameters(Parameter("rho", 1.0));

  this->diffusion_params_ = Parameters(Parameter("D", 1.0));
}
/**
 * @brief Overload the default options for density
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, Density DENS>
void DiffusionOperator<T, DIM, NLFI, DENS>::overload_density(const Parameters &p_params) {
  this->density_params_ = p_params;
}

/**
 * @brief Overload the default options for diffusion
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param p_params
 */
template <class T, int DIM, class NLFI, Density DENS>
void DiffusionOperator<T, DIM, NLFI, DENS>::overload_diffusion(const Parameters &p_params) {
  this->diffusion_params_ = p_params;
}
