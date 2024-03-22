/*
 * Copyright © CEA 2023
 *
 * \brief Spatial discretization used by phase-field models
 *
 * \file SpatialDiscretization.hpp
 * \author ci230846
 * \date 24/03/2023
 */

/*!
 * \page __spatial SpatialDiscretization
 * \section _spatial_sec1 Description
 */

#pragma once
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"

/**
 * @brief specialized_spatial_constructor
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
struct specialized_spatial_constructor {};

/**
 * @brief Class for defining the spatial discretization with a mesh coming from a GMSH file or built
 * by using dedicated  MFEM methods
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
class SpatialDiscretization {
 private:
  int size_;
  std::string existing_mesh_name_;
  mfem::FiniteElementSpace *fespace_;
  T *fecollection_;

 public:
  SpatialDiscretization(const std::string &mesh_type, const int &fe_order,
                        const std::string &mesh_file, bool periodic_mesh = false) {
    specialized_spatial_constructor<T, DIM> init;
    init(*this, mesh_type, fe_order, mesh_file, periodic_mesh);
  }

  template <class... Args>
  explicit SpatialDiscretization(const std::string &mesh_type, const int &fe_order,
                                 std::tuple<Args...> tup_args) {
    specialized_spatial_constructor<T, DIM> init;
    init(*this, mesh_type, fe_order, tup_args);
  }

  template <class... Args>
  explicit SpatialDiscretization(const std::string &mesh_type, const int &fe_order,
                                 std::tuple<Args...> tup_args,
                                 std::vector<mfem::Vector> translations) {
    specialized_spatial_constructor<T, DIM> init;
    init(*this, mesh_type, fe_order, tup_args, translations);
  }

  int fe_order_;
  int dimension_;
  mfem::Mesh mesh_;
  int mesh_max_bdr_attributes_;
  bool is_periodic_mesh_ = {false};

  void set_finite_element_space();

  mfem::Mesh &get_mesh();
  mfem::FiniteElementSpace *get_finite_element_space();

  std::size_t getSize();
  std::size_t get_max_bdr_attributes();
  int get_dimension();

  void apply_uniform_refinement(const int &level);

  void make_periodic_mesh(std::vector<mfem::Vector>);
  void make_periodic_line();

  bool is_periodic();

  ~SpatialDiscretization();
};

///////////////////////
// DIM = 1
///////////////////////

/**
 * @brief Specialization for dimension one
 *
 * @tparam T
 */
template <typename T>
struct specialized_spatial_constructor<T, 1> {
  /**
   * @brief operator() specialized for a mesh built from GMSH file
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param file
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 1> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, const std::string &file, bool periodic_mesh) {
    a_my_class.fe_order_ = fe_order;
    a_my_class.dimension_ = 1;

    switch (Meshes::from(mesh_type)) {
      case Meshes::GMSH: {
        if (std::filesystem::exists(file)) {
          const char *mesh_file = file.c_str();
          a_my_class.mesh_ = mfem::Mesh::LoadFromFile(mesh_file, 1, 1);
          break;
        } else {
          throw std::runtime_error("SpatialDiscretization::SpatialDiscretization: " + file +
                                   " doesn't exist. Please check your data.");
        }
      }
      default:
        throw std::runtime_error(
            "SpatialDiscretization::SpatialDiscretization: here, only GMSH mesh type is "
            "allowed");
        break;
    }

    a_my_class.is_periodic_mesh_ = periodic_mesh;
    if (!a_my_class.is_periodic_mesh_) {
      a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_.bdr_attributes.Max();
    } else {
      a_my_class.mesh_max_bdr_attributes_ = -1;
    }

    a_my_class.set_finite_element_space();
  }

  /**
   * @brief operator() specialized for a mesh built from MFEM methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 1> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);
    a_my_class.set_finite_element_space();
  }

  /**
   * @brief operator() specialized for a periodic mesh built from MFEM methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 1> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args,
                  std::vector<mfem::Vector> translations) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.make_periodic_line();

    a_my_class.is_periodic_mesh_ = true;
    a_my_class.set_finite_element_space();
  }

  /**
   * @brief Build a one-dimensional mesh by using MFEM inline methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void build_mesh(SpatialDiscretization<T, 1> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args) {
    a_my_class.fe_order_ = fe_order;
    auto tup_size = std::tuple_size<decltype(tup_args)>::value;
    a_my_class.dimension_ = 1;

    switch (Meshes::from(mesh_type)) {
      case Meshes::InlineLineWithSegments: {
        if (tup_size == 2) {
          const auto nx = std::get<0>(tup_args);
          const auto sx = std::get<1>(tup_args);
          a_my_class.mesh_ = mfem::Mesh::MakeCartesian1D(nx, sx);
        } else {
          throw std::runtime_error(
              "SpatialDiscretization::SpatialDiscretization: InlineLineWithSegments "
              "requires "
              "two "
              "argument, the number of segments");
        }
        break;
      }
      default:
        throw std::runtime_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineLineWithSegments, "
            "InlineSquareWithTriangles, InlineSquareWithQuadrangles mesh types are allowed");
        break;
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_.bdr_attributes.Max();
  }
};

///////////////////////
// DIM = 2
///////////////////////

/**
 * @brief Specialization for dimension two
 *
 * @tparam T
 */
template <typename T>
struct specialized_spatial_constructor<T, 2> {
  /**
   * @brief operator() specialized for a mesh built from GMSH file
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param file
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 2> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, const std::string &file, bool periodic_mesh) {
    a_my_class.fe_order_ = fe_order;
    a_my_class.dimension_ = 2;

    switch (Meshes::from(mesh_type)) {
      case Meshes::GMSH: {
        if (std::filesystem::exists(file)) {
          const char *mesh_file = file.c_str();
          a_my_class.mesh_ = mfem::Mesh::LoadFromFile(mesh_file, 1, 1);
          break;
        } else {
          throw std::runtime_error("SpatialDiscretization::SpatialDiscretization: " + file +
                                   " doesn't exist. Please check your data.");
        }
      }
      default:
        throw std::runtime_error(
            "SpatialDiscretization::SpatialDiscretization: here, only GMSH mesh type is "
            "allowed");
        break;
    }
    a_my_class.is_periodic_mesh_ = periodic_mesh;
    if (!a_my_class.is_periodic_mesh_) {
      a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_.bdr_attributes.Max();
    } else {
      a_my_class.mesh_max_bdr_attributes_ = -1;
    }
    a_my_class.set_finite_element_space();
  }

  /**
   * @brief operator() specialized for a mesh built from MFEM methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 2> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);
    a_my_class.set_finite_element_space();
  }

  /**
   * @brief operator() specialized for a periodic mesh built from MFEM methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   * @param translations
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 2> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args,
                  std::vector<mfem::Vector> translations) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.make_periodic_mesh(translations);
    a_my_class.is_periodic_mesh_ = true;
    a_my_class.set_finite_element_space();
  }

  /**
   * @brief Build a two-dimensional mesh by using MFEM inline methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void build_mesh(SpatialDiscretization<T, 2> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args) {
    a_my_class.fe_order_ = fe_order;
    auto tup_size = std::tuple_size<decltype(tup_args)>::value;
    a_my_class.dimension_ = 2;

    auto element = mfem::Element::QUADRILATERAL;
    switch (Meshes::from(mesh_type)) {
      case Meshes::InlineSquareWithQuadrangles: {
        element = mfem::Element::QUADRILATERAL;
        break;
      }
      case Meshes::InlineSquareWithTriangles: {
        element = mfem::Element::TRIANGLE;
        break;
      }
      default:
        throw std::runtime_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineLineWithSegments, "
            "InlineSquareWithTriangles, InlineSquareWithQuadrangles mesh types are allowed");
        break;
    }
    if (tup_size == 4) {
      const auto nx = std::get<0>(tup_args);
      const auto ny = std::get<1>(tup_args);
      const auto sx = std::get<2>(tup_args);
      const auto sy = std::get<3>(tup_args);
      a_my_class.mesh_ = mfem::Mesh::MakeCartesian2D(nx, ny, element, false, sx, sy, false);

    } else {
      throw std::runtime_error(
          "SpatialDiscretization::SpatialDiscretization: InlineSquareWithQuadrangles "
          "requires "
          "four arguments, the number of nodes and the length along each direction");
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_.bdr_attributes.Max();
  }
};

///////////////////////
// DIM = 3
///////////////////////

/**
 * @brief Specialization for dimension three
 *
 * @tparam T
 */
template <typename T>
struct specialized_spatial_constructor<T, 3> {
  /**
   * @brief operator() specialized for a mesh built from GMSH file
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param file
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 3> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, const std::string &file, bool periodic_mesh) {
    a_my_class.fe_order_ = fe_order;
    a_my_class.dimension_ = 3;

    switch (Meshes::from(mesh_type)) {
      case Meshes::GMSH: {
        if (std::filesystem::exists(file)) {
          const char *mesh_file = file.c_str();
          a_my_class.mesh_ = mfem::Mesh::LoadFromFile(mesh_file, 1, 1);
          break;
        } else {
          throw std::runtime_error("SpatialDiscretization::SpatialDiscretization: " + file +
                                   " doesn't exist. Please check your data.");
        }
      }
      default:
        throw std::runtime_error(
            "SpatialDiscretization::SpatialDiscretization: here, only GMSH mesh type is "
            "allowed");
        break;
    }
    a_my_class.is_periodic_mesh_ = periodic_mesh;
    if (!a_my_class.is_periodic_mesh_) {
      a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_.bdr_attributes.Max();
    } else {
      a_my_class.mesh_max_bdr_attributes_ = -1;
    }

    a_my_class.set_finite_element_space();
  }

  /**
   * @brief operator() specialized for a mesh built from MFEM methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 3> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.set_finite_element_space();
  }

  /**
   * @brief operator() specialized for a periodic mesh built from MFEM methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void operator()(SpatialDiscretization<T, 3> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args,
                  std::vector<mfem::Vector> translations) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.make_periodic_mesh(translations);
    a_my_class.is_periodic_mesh_ = true;
    a_my_class.set_finite_element_space();
  }

  /**
   * @brief Build a three-dimensional mesh by using MFEM inline methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void build_mesh(SpatialDiscretization<T, 3> &a_my_class, const std::string &mesh_type,
                  const int &fe_order, std::tuple<Args...> tup_args) {
    a_my_class.fe_order_ = fe_order;
    auto tup_size = std::tuple_size<decltype(tup_args)>::value;
    a_my_class.dimension_ = 3;

    auto element = mfem::Element::TETRAHEDRON;
    switch (Meshes::from(mesh_type)) {
      case Meshes::InlineSquareWithTetraedres: {
        element = mfem::Element::TETRAHEDRON;
        break;
      }
      case Meshes::InlineSquareWithHexaedres: {
        element = mfem::Element::HEXAHEDRON;
        break;
      }
      default:
        throw std::runtime_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineLineWithSegments, "
            "InlineSquareWithTriangles, InlineSquareWithQuadrangles mesh types are allowed");
        break;
    }

    if (tup_size == 6) {
      a_my_class.dimension_ = 3;
      const auto nx = std::get<0>(tup_args);
      const auto ny = std::get<1>(tup_args);
      const auto nz = std::get<2>(tup_args);
      const auto sx = std::get<3>(tup_args);
      const auto sy = std::get<4>(tup_args);
      const auto sz = std::get<5>(tup_args);
      a_my_class.mesh_ = mfem::Mesh::MakeCartesian3D(nx, ny, nz, element, sx, sy, sz);
    } else {
      throw std::runtime_error(
          "SpatialDiscretization::SpatialDiscretization: " + mesh_type +
          " requires six arguments, the number of nodes and the length along each direction");
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_.bdr_attributes.Max();
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief  return a pointer of Mesh
 *
 * @tparam T
 * @return mfem::Mesh&
 */
template <class T, int DIM>
mfem::Mesh &SpatialDiscretization<T, DIM>::get_mesh() {
  return this->mesh_;
}

/**
 * @brief Set the FE_Collection, the FE_Space and associated size
 *
 * @tparam T
 * @return mfem::FiniteElementSpace*
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::set_finite_element_space() {
  this->fecollection_ = new T(this->fe_order_, this->dimension_);
  this->fespace_ = new mfem::FiniteElementSpace(&this->mesh_, this->fecollection_);
  this->size_ = this->fespace_->GetTrueVSize();
}

/**
 * @brief return a pointer toward the finite element space
 *
 * @tparam T
 * @return mfem::FiniteElementSpace*
 */
template <class T, int DIM>
mfem::FiniteElementSpace *SpatialDiscretization<T, DIM>::get_finite_element_space() {
  return this->fespace_;
}

/**
 * @brief get the size of the Finite Element Space
 *
 * @tparam T
 * @return int
 */
template <class T, int DIM>
std::size_t SpatialDiscretization<T, DIM>::getSize() {
  return this->size_;
}

/**
 * @brief get the maximum number of boundaries
 *
 * @tparam T
 * @return int
 */
template <class T, int DIM>
std::size_t SpatialDiscretization<T, DIM>::get_max_bdr_attributes() {
  return this->mesh_max_bdr_attributes_;
}

/**
 * @brief get the dimension of the problem
 *
 * @tparam T
 * @return int
 */
template <class T, int DIM>
int SpatialDiscretization<T, DIM>::get_dimension() {
  return this->dimension_;
}

/**
 * @brief Apply nb_ref uniform refinement
 *
 * @tparam T
 * @param nb_ref
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::apply_uniform_refinement(const int &nb_ref) {
  for (auto l = 0; l < nb_ref; l++) {
    this->mesh_.UniformRefinement();
  }
}

/**
 * @brief Create the periodic mesh using the vertex mapping defined by the translations vector
 *
 * @tparam T
 * @tparam DIM
 * @param translations
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::make_periodic_mesh(std::vector<mfem::Vector> translations) {
  const auto tol = 1.e-6;
  std::vector<int> periodicMap = this->mesh_.CreatePeriodicVertexMapping(translations, tol);
  auto periodic_mesh = mfem::Mesh::MakePeriodic(this->mesh_, periodicMap);
  this->mesh_ = mfem::Mesh(periodic_mesh, true);  // replace the input mesh with the periodic one
}

/**
 * @brief Create a periodic line segment
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::make_periodic_line() {
  // Based on mfem.org example
  // Create the vertex mapping. To begin, create the identity mapping.
  std::vector<int> periodicMap(this->mesh_.GetNV());
  for (int i = 0; i < periodicMap.size(); ++i) {
    periodicMap[i] = i;
  }
  // Modify the mapping so that the last vertex gets mapped to the first vertex.
  periodicMap.back() = 0;
  auto periodic_mesh = mfem::Mesh::MakePeriodic(this->mesh_, periodicMap);
  this->mesh_ = mfem::Mesh(periodic_mesh, true);  // replace the input mesh with the periodic one
}

/**
 * @brief Return the flag to know whether the mesh is periodic or not
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
bool SpatialDiscretization<T, DIM>::is_periodic() {
  return this->is_periodic_mesh_;
}

/**
 * @brief Destroy the Spatial Discretization< T>:: Spatial Discretization object
 *
 * @tparam T
 */
template <class T, int DIM>
SpatialDiscretization<T, DIM>::~SpatialDiscretization() {}
