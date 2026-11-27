/**
 * @file MeltingCalphadNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief NonlinearFormIntegrator for ad-hoc melting term based on temperature
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#pragma once
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/AuxiliaryCoefficient.hpp"
#include "Integrators/MeltingBaseNLFormIntegrator.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Phase change of MeltingCalphadNLFormIntegrator: melting factor times the difference
 * of the driving forces of the secondary and primary phases, zero if one of them is zero
 * (|dgm| <= 1e-12).
 */
struct MeltingCalphadPhaseChangeEval {
  static constexpr int nb_variables = 1;
  static constexpr int nb_explicit_variables = 0;
  static constexpr int nb_auxiliary_variables = 0;  ///< read by index: see max_auxiliary_index
  static constexpr int nb_coordinates = 0;
  int primary_index = 0;    ///< index of the driving force of the primary phase
  int secondary_index = 0;  ///< index of the driving force of the secondary phase
  double melting_factor = 1.0;

  int max_auxiliary_index() const { return std::max(primary_index, secondary_index); }

  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,
                                   [[maybe_unused]] const double* explicit_values,
                                   const double* auxiliary_values,
                                   [[maybe_unused]] const double* coordinates,
                                   [[maybe_unused]] const double time, double* F, double* gradient,
                                   double* hessian) const {
    if (F != nullptr) {
      const double primary_dgm = auxiliary_values[primary_index];
      const double secondary_dgm = auxiliary_values[secondary_index];
      constexpr double epsilon = 1.e-12;
      *F = (fabs(primary_dgm) > epsilon && fabs(secondary_dgm) > epsilon)
               ? melting_factor * (secondary_dgm - primary_dgm)
               : 0.;
    }
    if (gradient != nullptr) gradient[0] = 0.0;
    if (hessian != nullptr) hessian[0] = 0.0;
  }
};

/**
 * @brief Class dedicated to the VF of Calphad melting contribution found in Allen-Cahn-type
 * equations.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class MeltingCalphadNLFormIntegrator : public MeltingBaseNLFormIntegrator<VARS> {
 private:
  double alpha_;
  std::string primary_phase_;
  std::string secondary_phase_;
  int dgm_primary_phase_index_;
  int dgm_secondary_phase_index_;
  int primary_nucleus_index_;
  int secondary_nucleus_index_;

  double melting_temperature_;
  double melting_enthalpy_;
  void get_parameters();

  void check_driving_forces();

  void check_nucleus();

  mfem::Vector Psi_;
  double get_seed_at_ip(const std::span<const double>& aux_values) const;

 public:
  MeltingCalphadNLFormIntegrator(Geometry geometry, const double time_step,
                                 const std::vector<mfem::ParGridFunction>& u_old,
                                 const std::vector<mfem::ParGridFunction>& aux_old,
                                 const Parameters& params, std::vector<VARS*> auxvars,
                                 const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvec) override;

  virtual ~MeltingCalphadNLFormIntegrator() = default;
};

#include "Integrators/MeltingCalphadNLFormIntegrator.tpp"
