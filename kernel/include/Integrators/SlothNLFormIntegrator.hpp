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
#include <list>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/Coefficients.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
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
  std::string integrator_name_ = "";
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

  std::vector<mfem::ParGridFunction> get_aux_gf();
  std::vector<std::vector<std::string>> get_aux_infos();

  void check_coefficient_types(std::list<GlossaryType> expected_type);
  std::optional<Coefficient> get_coefficient(const int blk, GlossaryType type, unsigned int id,
                                             std::optional<int> bdr_id = std::nullopt);

  double compute_coefficient(Coefficient coef, const std::span<const double>& values,
                             const std::span<const double>& aux_values);

  double compute_gradient_coefficient(Coefficient coef, const int blk,
                                      const std::span<const double>& values,
                                      const std::span<const double>& aux_values);
  double compute_hessian_coefficient(Coefficient coef, const int iblk, const int jblk,
                                     const std::span<const double>& values,
                                     const std::span<const double>& aux_values);

 public:
  virtual void init() = 0;
  SlothNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                        const std::vector<mfem::ParGridFunction> aux_old, const Parameters& params,
                        std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients);
  virtual ~SlothNLFormIntegrator() = default;
};

#include "Integrators/SlothNLFormIntegrator.tpp"
