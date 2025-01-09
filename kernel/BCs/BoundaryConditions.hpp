/*
 * Copyright © CEA 2023
 *
 * \brief BoundaryConditions class used to build and manage boundary conditions
 *
 * \file BoundaryConditions.hpp
 * \author ci230846
 * \date 23/03/2023
 */
#include <limits>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "BCs/Boundary.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

/**
 * @brief Class used to manage boundary conditions
 *
 */
template <class T, int DIM>
class BoundaryConditions {
 private:
  T *fecollection_;
  mfem::ParFiniteElementSpace *fespace_;
  mfem::Array<int> Dirichlet_bdr_;
  mfem::Array<double> Dirichlet_value_;
  mfem::Array<int> ess_tdof_list_;
  std::initializer_list<Boundary> boundaries_;

 public:
  template <class... Args>
  BoundaryConditions(SpatialDiscretization<T, DIM> *spatial, const Args &...boundaries);
  void SetBoundaryConditions(mfem::Vector &u);
  mfem::Array<int> GetEssentialDofs();
  ~BoundaryConditions();
};

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
BoundaryConditions<T, DIM>::BoundaryConditions(SpatialDiscretization<T, DIM> *spatial,
                                               const Args &...boundaries) {
  bool must_be_periodic = spatial->is_periodic();

  this->fespace_ = spatial->get_finite_element_space();

  auto bdrs = std::vector<Boundary>{boundaries...};
  auto mesh_max_bdr_attributes = bdrs.size();
  // If periodic, the number of boundary conditions is used
  // Consistency is checked later
  if (!must_be_periodic) {
    mesh_max_bdr_attributes = spatial->get_max_bdr_attributes();
  }

  Dirichlet_bdr_.SetSize(mesh_max_bdr_attributes);
  Dirichlet_value_.SetSize(mesh_max_bdr_attributes);
  bool exist_periodic_bdr = false;
  for (const auto &bdr : bdrs) {
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

    for (const auto &bdr : bdrs) {
      const auto &id = bdr.get_boundary_index();
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

      if (bdr.is_essential_boundary()) {
        Dirichlet_bdr_[id] = 1;
      } else {
        Dirichlet_bdr_[id] = 0;
      }
      Dirichlet_value_[id] = bdr.get_boundary_value();
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
 * @brief Set boundary conditions
 *
 * @param u unknown vector
 */
template <class T, int DIM>
void BoundaryConditions<T, DIM>::SetBoundaryConditions(mfem::Vector &u) {
  mfem::Array<int> tmp_array_bdr(this->Dirichlet_bdr_.Size());
  for (auto i = 0; i < this->Dirichlet_bdr_.Size(); i++) {
    tmp_array_bdr = 0;
    mfem::Array<int> dof;
    if (this->Dirichlet_bdr_[i] > 0) {
      tmp_array_bdr[i] = 1;
      this->fespace_->GetEssentialTrueDofs(tmp_array_bdr, dof);
      u.SetSubVector(dof, this->Dirichlet_value_[i]);
    }
  }
}  // end of SetBoundaryConditions

/**
 * @brief Destroy the Boundary Conditions:: Boundary Conditions object
 *
 */
template <class T, int DIM>
BoundaryConditions<T, DIM>::~BoundaryConditions() {}
