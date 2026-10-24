/**
 * @file SlothNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class prodiving common methods to all NLFormIntegrators
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * @anchor SlothNLFormIntegrator
 *
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
#include <cstddef>
#include <list>
#include <optional>
#include <span>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/Coefficients.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Options/ProblemsOptions.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class prodiving common methods to all NLFormIntegrators
 *
 */
template <class VARS>
class SlothNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator {
 private:
  void manage_auxiliary_variables(std::vector<VARS*> auxvars);
  std::vector<mfem::ParGridFunction> vect_aux_gf_;
  std::vector<mfem::Vector> vect_aux_old_gf_;
  std::vector<std::vector<std::string>> vect_aux_infos_;

 protected:
  Geometry geometry_ = Geometry::Cartesian;
  double time_step_;
  std::string integrator_name_ = "";

  // Current integration point (used by file-based coefficients)
  mfem::ElementTransformation* current_tr_ = nullptr;  ///< Transformation of the current element
  mfem::Vector x_ip_;        ///< Physical coordinates of the current integration point
  bool x_ip_valid_ = false;  ///< True if x_ip_ matches the current integration point

  void set_int_point(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip);
  std::span<const double> get_physical_point();

  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvect) = 0;

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmat) = 0;
  virtual void get_coefficients() = 0;

  std::vector<mfem::ParGridFunction> u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  std::vector<mfem::ParGridFunction> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_infos_;

  std::vector<VARS*> auxvariables_;
  Parameters params_;
  std::vector<Coefficients> coefficients_;
  unsigned int nb_blk_;
  unsigned int nb_vaux_;

  std::vector<mfem::ParGridFunction> get_aux_gf();
  std::vector<std::vector<std::string>> get_aux_infos();

  void check_coefficient_types(std::list<GlossaryType> expected_type);
  std::optional<Coefficient> get_coefficient(const int blk, GlossaryType type, unsigned int id,
                                             std::optional<int> bdr_id = std::nullopt);

  // Evaluation of a coefficient and of its derivatives at the current integration point.
  // The coefficient is passed by reference (no copy at each integration point) and may be
  // modified (coordinates set before evaluation).
  // var_blk: block of the variable of a coefficient of one variable (W, H...: see
  // FunctionCoefficient::is_single_variable); such a coefficient is evaluated on u[var_blk]
  // only, and its derivatives with respect to the other blocks are zero. Ignored for the
  // other coefficients, which receive all the blocks.
  double eval_coefficient(Coefficient& coef, const std::span<const double>& values,
                          const std::span<const double>& aux_values, const int var_blk = 0);
  double eval_gradient_coefficient(Coefficient& coef, const int blk,
                                   const std::span<const double>& values,
                                   const std::span<const double>& aux_values,
                                   const int var_blk = 0);
  double eval_gradient_energy_coefficient(Coefficient& coef, const int blk,
                                          const std::span<const double>& values,
                                          const std::span<const double>& aux_values,
                                          const int var_blk = 0);
  double eval_hessian_coefficient(Coefficient& coef, const int iblk, const int jblk,
                                  const std::span<const double>& values,
                                  const std::span<const double>& aux_values,
                                  const int var_blk = 0);

 private:
  /// Values given to a coefficient: all the blocks {u, u_old}, or block var_blk only for a
  /// coefficient of one variable
  std::pair<std::span<const double>, std::span<const double>> coefficient_inputs(
      const Coefficient& coef, const std::span<const double>& values, const int var_blk) const;

 protected:

  bool isAxisymmetric() const noexcept;

  //////////////////////////////////////////////////////////////////////////////
  // Evaluation of {u, u_old} and of the auxiliary fields at integration points
  //////////////////////////////////////////////////////////////////////////////
  /// Values returned by values_at_ip()
  struct IpValues {
    std::span<const double> u;    ///< [u_0..u_{n-1}, u_old_0..u_old_{n-1}]
    std::span<const double> aux;  ///< auxiliary fields
  };

  static bool same_finite_element(const mfem::FiniteElement& a, const mfem::FiniteElement& b);
  static bool same_finite_elements(const mfem::Array<const mfem::FiniteElement*>& el);

  void prepare_element_values(const mfem::Array<const mfem::FiniteElement*>& el,
                              mfem::ElementTransformation& Tr,
                              const mfem::Array<const mfem::Vector*>& elfun,
                              const std::vector<mfem::ParGridFunction>& aux_gf,
                              const mfem::IntegrationRule& ir);

  IpValues values_at_ip(const mfem::IntegrationRule& ir, const int i, const mfem::Vector& Psi,
                        const mfem::Array<const mfem::Vector*>& elfun,
                        const std::vector<mfem::ParGridFunction>& aux_gf,
                        mfem::ElementTransformation& Tr);

  void element_gradient(const mfem::ParGridFunction& gf, const mfem::ElementTransformation& Tr,
                        const mfem::DenseMatrix& gradPsi, mfem::Vector& grad);

 private:
  mfem::Vector grad_dofs_;  ///< work vector of element_gradient()

  // Work data of prepare_element_values() / values_at_ip()
  const mfem::IntegrationRule* ipv_ir_ = nullptr;  ///< rule of the cached table (nullptr: none)
  int ipv_element_ = -1;                           ///< element of the cached table
  int ipv_nb_ = 0;                                 ///< number of blocks of the cached table
  std::size_t ipv_naux_ = 0;                       ///< number of auxiliary fields of the table
  std::vector<double> ipv_u_table_;                ///< {u, u_old} at each point (nip x 2 nb)
  std::vector<double> ipv_aux_table_;              ///< auxiliary fields at each point (nip x naux)
  std::vector<double> ipv_u_;                      ///< {u, u_old} at one point (no table)
  std::vector<double> ipv_aux_;                    ///< auxiliary fields at one point (no table)
  std::vector<mfem::Vector> ipv_uold_dofs_;        ///< element dofs of u_old (one per block)
  std::vector<mfem::Vector> ipv_aux_dofs_;         ///< element dofs of the auxiliary fields
  std::vector<char> ipv_aux_dot_;                  ///< aux k evaluated from its element dofs
  mfem::Vector ipv_shape_;                         ///< shape functions at one point

 public:
  virtual void init() = 0;
  SlothNLFormIntegrator(Geometry geometry, const double time_step,
                        const std::vector<mfem::ParGridFunction> u_old,
                        const std::vector<mfem::ParGridFunction> aux_old, const Parameters& params,
                        std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients);
  virtual ~SlothNLFormIntegrator() = default;
};

#include "Integrators/SlothNLFormIntegrator.tpp"
