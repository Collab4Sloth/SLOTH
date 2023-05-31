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
#include <vector>
#include "BCs/Boundary.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"
#pragma once

/**
 * @brief Class used to manage boundary conditions
 *
 */
template <class T, int DIM>
class BoundaryConditions {
 private:
  T *fecollection_;
  mfem::FiniteElementSpace *fespace_;
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
  this->fespace_ = spatial->get_finite_element_space();
  const auto &mesh_max_bdr_attributes = spatial->get_max_bdr_attributes();

  auto bdrs = std::vector<Boundary>{boundaries...};

  Dirichlet_bdr_.SetSize(mesh_max_bdr_attributes);
  Dirichlet_value_.SetSize(mesh_max_bdr_attributes);
  bool exist_periodic_bdr = {false};
  [bdrs, &exist_periodic_bdr]() {
    bool is_periodic = false;
    for (const auto &bdr : bdrs) {
      if (bdr.is_essential_boundary()) {
        is_periodic = true;
        break;
      }
    }
    return is_periodic;
  };
  bool test_periodic_bdr = (mesh_max_bdr_attributes != bdrs.size()) && exist_periodic_bdr;
  bool test_standard_bdr = mesh_max_bdr_attributes == bdrs.size();
  if (test_periodic_bdr || test_standard_bdr) {
    for (const auto &bdr : bdrs) {
      const auto &id = bdr.get_boundary_index();
      if (bdr.is_essential_boundary()) {
        Dirichlet_bdr_[id] = 1;
      } else {
        Dirichlet_bdr_[id] = 0;
      }
      Dirichlet_value_[id] = bdr.get_boundary_value();
    }
    this->fespace_->GetEssentialTrueDofs(this->Dirichlet_bdr_, this->ess_tdof_list_);
  } else {
    throw std::runtime_error(
        "BoundaryConditions::BoundaryConditions(): user-defined boundaries  " +
        std::to_string(bdrs.size()) +
        " are unconsistent with the total number of boundaries associated to the mesh " +
        std::to_string(mesh_max_bdr_attributes));
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
