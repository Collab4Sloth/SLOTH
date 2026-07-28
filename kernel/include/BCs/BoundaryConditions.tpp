/**
 * @file BoundaryConditions.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief BoundaryConditions class used to build and manage boundary conditions
 * @version 0.1
 * @date 2025-09-05
 *
 * @anchor bcs
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

#include <limits>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "BCs/Boundary.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Boundary Conditions:: Boundary Conditions object
 *
 * @tparam Args
 * @param fespace
 * @param mesh_max_bdr_attributes
 * @param boundaries
 */
template <class T, int DIM>
template <class... Args>
BoundaryConditions<T, DIM>::BoundaryConditions(SpatialDiscretization<T, DIM>* spatial,
                                               const Args&... boundaries) {
  bool must_be_periodic = spatial->is_periodic();

  this->fespace_ = spatial->get_finite_element_space();

  auto bdrs = std::vector<Boundary>{boundaries...};
  auto mesh_max_bdr_attributes = bdrs.size();
  // If periodic, the number of boundary conditions is used
  // Consistency is checked later
  if (!must_be_periodic) {
    mesh_max_bdr_attributes = spatial->get_max_bdr_attributes();
  }

  Neumann_bdr_.SetSize(mesh_max_bdr_attributes);
  Robin_bdr_.SetSize(mesh_max_bdr_attributes);
  Dirichlet_bdr_.SetSize(mesh_max_bdr_attributes);
  Dirichlet_value_.SetSize(mesh_max_bdr_attributes);
  bool exist_periodic_bdr = false;
  for (const auto& bdr : bdrs) {
    if (bdr.is_periodic_boundary()) {
      exist_periodic_bdr = true;
      break;
    }
  }

  if (!must_be_periodic && exist_periodic_bdr) {
    mfem::mfem_error(
        "BoundaryConditions::BoundaryConditions(): mesh is not defined as periodic but at least "
        "one boundary is flagged periodic. Please check your data");
  }

  if (must_be_periodic && !exist_periodic_bdr) {
    mfem::mfem_error(
        "BoundaryConditions::BoundaryConditions(): mesh is defined as periodic but no boundary is "
        "flagged periodic. Please check your data");
  }

  bool test_standard_bdr = mesh_max_bdr_attributes == bdrs.size();

  if (exist_periodic_bdr || test_standard_bdr) {
    std::unordered_set<int> id_seen;

    for (const auto& bdr : bdrs) {
      const auto& id = bdr.get_boundary_index();
      // Check index value
      if (id >= static_cast<int>(mesh_max_bdr_attributes)) {
        std::string msg = "BoundaryConditions::BoundaryConditions(): bad index " +
                          std::to_string(id) + ". Index should be lower than " +
                          std::to_string(mesh_max_bdr_attributes);
        mfem::mfem_error(msg.c_str());
      }

      // Check unicity
      if (!id_seen.insert(id).second) {
        std::string msg = "BoundaryConditions::BoundaryConditions(): duplicated index " +
                          std::to_string(id) + ". Index should be unique";
        mfem::mfem_error(msg.c_str());
      }

      // Create marker arrays
      const auto& boundary_type_ = bdr.get_boundary_type();
      Dirichlet_bdr_[id] = 0;
      Neumann_bdr_[id] = 0;
      Robin_bdr_[id] = 0;
      if (boundary_type_ == BoundaryConditionType::Dirichlet) {
        Dirichlet_bdr_[id] = 1;
        Dirichlet_value_[id] = bdr.get_boundary_value();
      } else if (boundary_type_ == BoundaryConditionType::Neumann) {
        Neumann_bdr_[id] = 1;
      } else if (boundary_type_ == BoundaryConditionType::Robin) {
        Robin_bdr_[id] = 1;
      }
    }
    this->fespace_->GetEssentialTrueDofs(this->Dirichlet_bdr_, this->ess_tdof_list_);
  } else {
    std::string msg =
        "BoundaryConditions::BoundaryConditions(): user-defined boundaries  " +
        std::to_string(bdrs.size()) +
        " are unconsistent with the total number of boundaries associated to the mesh " +
        std::to_string(mesh_max_bdr_attributes);
    mfem::mfem_error(msg.c_str());
  }
}

/**
 * @brief return the list of essential dofs
 *
 * @return mfem::Array<int> array of essential dofs
 */
template <class T, int DIM>
mfem::Array<int> BoundaryConditions<T, DIM>::GetEssentialDofs() {
  return this->ess_tdof_list_;
}

/**
 * @brief return the marker array for the boundary condition
 *
 * @return mfem::Array<int> marker array for the boundary condition
 */
template <class T, int DIM>
mfem::Array<int> BoundaryConditions<T, DIM>::get_marker_array(const std::string& boundary_type) {
  if (boundary_type == "Neumann") {
    return this->Neumann_bdr_;
  }
  if (boundary_type == "Robin") {
    return this->Robin_bdr_;
  }

  mfem::Array<int> null_bdr(this->Neumann_bdr_.Size());
  null_bdr = 0;
  return null_bdr;
}

/**
 * @brief Set dirichlet boundary conditions
 *
 * @param u unknown vector
 * @param coefficients coefficients list for the variable
 * @param auxvars_unk unknown vectors of the auxiliary variables
 * 
 */
template <class T, int DIM>
void BoundaryConditions<T, DIM>::SetDirichletBoundaryConditions(mfem::Vector& u, Coefficients& coefficients, std::vector<mfem::Vector> auxvars_unk) {

  const int nb_bdr = this->Dirichlet_bdr_.Size();
  mfem::Array<int> tmp_array_bdr(this->Dirichlet_bdr_.Size());
  const int coef_size = coefficients.size();
  Coefficient dirichlet_coef = Coefficient(Glossary::Default, 0.0);

  std::vector<double> u_values;
  std::vector<double> vaux_values;

  // Inline function to compute Dirichlet coefficient
  auto compute_dirichlet_coefficient = [&](Coefficient& coef, 
    const std::span<const double>& values,
    const std::span<const double>& aux_values) -> double {
    if (coef.is_scalar()) {
      return coef.compute();
    } else {
      return coef.compute(values, aux_values);
    }
  };

  // Loop over boundaries
  for (int i = 0; i < nb_bdr; i++) {

    // If the boundary is dirichlet
    if (this->Dirichlet_bdr_[i] > 0) {

      // Check if there is a coefficient given and store it
      bool has_dirichlet_coef = false;
      for (int l = 0; l < coef_size; l++) {
        const auto& coef = coefficients[l];
        if (coef.get_type() == GlossaryType::Dirichlet) {
          auto bdr_ids = coef.get_bdr_index_coef();
          if (std::find(bdr_ids.begin(), bdr_ids.end(), i) != bdr_ids.end()) {
            dirichlet_coef = coef;
            has_dirichlet_coef = true;
            break;
          }
        }
      }

      // Get the list of essential true dofs
      tmp_array_bdr = 0;
      tmp_array_bdr[i] = 1;
      mfem::Array<int> dof_list;
      this->fespace_->GetEssentialTrueDofs(tmp_array_bdr, dof_list);

      if (has_dirichlet_coef) {
        mfem::Vector dirichlet_at_dofs(dof_list.Size());
        // Loop over essential dofs and compute the coefficient at the dofs
        for (int j = 0; j < dof_list.Size(); j++) {
          u_values.clear();
          vaux_values.clear();

          int dof = dof_list[j];
          u_values.push_back(u(dof));
          for (auto aux : auxvars_unk)
              vaux_values.emplace_back(aux(dof));

          dirichlet_at_dofs[j] = compute_dirichlet_coefficient(dirichlet_coef, std::span<const double>(u_values),
          std::span<const double>(vaux_values));
        }

        // SetSubVector with the calculated dirichlet values
        u.SetSubVector(dof_list, dirichlet_at_dofs);
      } else {
        // If no dirichlet coefficient is given, use the constant dirichlet value
        u.SetSubVector(dof_list, this->Dirichlet_value_[i]);
      }
    }
  }
}  // end of SetDirichletBoundaryConditions

/**
 * @brief Destroy the Boundary Conditions:: Boundary Conditions object
 *
 */
template <class T, int DIM>
BoundaryConditions<T, DIM>::~BoundaryConditions() {}
