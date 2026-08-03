/**
 * @file Spatial.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Spatial discretization
 * @version 0.1
 * @date 2025-09-05
 *
 * @anchor meshing
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
#include <filesystem>  // NOLINT [avoid  <filesystem> is an unapproved C++17 header.]
#include <functional>
#include <memory>
#include <regex>  // NOLINT [avoid  <filesystem> is an unapproved C++11 header.]
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Options/Options.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Helper structure for handling mesh file patterns and counts.
 *
 * This structure provides functionality to count the number of mesh files
 * in the current directory that match a given prefix pattern, and to compare
 * the count to an expected number.
 */
struct split_mesh_helper {
  /**
   * @brief Counts the number of files in the current directory matching the given prefix.
   *
   * @param mesh_prefix
   * @return int
   */
  int count_mesh_files(const std::string& mesh_prefix) {
    namespace fs = std::filesystem;
    std::regex pattern("^" + mesh_prefix + ".*$");

    int count = 0;

    for (const auto& entry : fs::directory_iterator(".")) {
      if (entry.is_regular_file()) {
        const std::string filename = entry.path().filename().string();
        if (std::regex_match(filename, pattern)) {
          ++count;
        }
      }
    }
    return count;
  }

  /**
   * @brief Checks whether the number of mesh files matches the expected count.
   *
   * @param n_files
   * @param mesh_pattern
   * @return true
   * @return false
   */
  bool operator()(int n_files, std::string mesh_pattern) {
    return (count_mesh_files(mesh_pattern) == n_files);
  }
};

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
  HYPRE_BigInt size_;
  std::string existing_mesh_name_;
  mfem::ParFiniteElementSpace* fespace_;
  T* fecollection_;
  bool owns_mesh_{true};
  bool is_nc_simplices_ = {false};

 public:
  /**
   * @brief Construct a new SpatialDiscretization object from a mesh file
   *        (e.g. GMSH), or using a dedicated MFEM mesh generation method
   *        identified by `mesh_type`.
   *
   * @details Delegates the actual mesh construction to
   *          `specialized_spatial_constructor<T, DIM>`, which reads/builds
   *          the underlying serial mesh according to `mesh_type` and
   *          `mesh_file`, calls `EnsureNCMesh(is_nc_simplices_)` before
   *          finalizing it, and converts it to a distributed `mfem::ParMesh`.
   *          Mesh attributes (element/boundary attribute sets) found in the
   *          file are then loaded via `set_mesh_attributes_from_file()`.
   *
   * @tparam T Finite element collection type.
   * @tparam DIM Spatial dimension.
   *
   * @param mesh_type Identifier of the mesh source/generation strategy (e.g.
   *                  a GMSH file type, or a dedicated MFEM inline mesh
   *                  builder).
   * @param fe_order Finite element order used to build the finite element
   *                collection and space.
   * @param ref_level Number of uniform refinement levels applied at
   *                  construction time, before any AMR.
   * @param mesh_file Path (or naming pattern) of the mesh file to read.
   * @param periodic_mesh Whether the mesh should be treated as periodic.
   * @param allow_nc_simplices Whether nonconforming refinement is allowed on
   *                          simplex elements (see `is_nc_simplices()`). Must
   *                          be set correctly here since the underlying
   *                          `ParMesh` cannot be converted to nonconforming
   *                          after construction.
   */
  SpatialDiscretization(const std::string& mesh_type, const int& fe_order, const int& ref_level,
                        const std::string& mesh_file, bool periodic_mesh = false,
                        bool allow_nc_simplices = false) {
    specialized_spatial_constructor<T, DIM> init;
    this->is_nc_simplices_ = allow_nc_simplices;
    init(*this, mesh_type, fe_order, ref_level, mesh_file, periodic_mesh);
    this->set_mesh_attributes_from_file();
  }

  /**
   * @brief Construct a new SpatialDiscretization object
   *
   * @details Delegates the actual mesh construction to
   *          `specialized_spatial_constructor<T, DIM>`, which dispatches on
   *          `mesh_type` to build the underlying serial mesh, call
   *          `EnsureNCMesh(is_nc_simplices_)`
   *          before finalizing it, and convert it to a distributed
   *          `mfem::ParMesh`.
   *
   * @tparam Args Types of the extra arguments forwarded to the specialized
   *             constructor (mesh dimensions, lengths, etc.), packed in
   *             `tup_args`.
   * @tparam T Finite element collection type.
   * @tparam DIM Spatial dimension.
   *
   * @param mesh_type Identifier of the mesh generation strategy to use (e.g.
   *                  "InlineSquareWithQuadrangles").
   * @param fe_order Finite element order used to build the finite element
   *                collection and space.
   * @param ref_level Number of uniform refinement levels applied at
   *                  construction time, before any AMR.
   * @param tup_args Tuple of extra parameters specific to `mesh_type` (e.g.
   *                number of elements and physical extents per direction).
   * @param allow_nc_simplices Whether nonconforming refinement is allowed on
   *                          simplex elements (see `is_nc_simplices()`). Must
   *                          be set correctly here since the underlying
   *                          `ParMesh` cannot be converted to nonconforming
   *                          after construction.
   */
  template <class... Args>
  explicit SpatialDiscretization(const std::string& mesh_type, const int& fe_order,
                                 const int& ref_level, std::tuple<Args...> tup_args,
                                 bool allow_nc_simplices = false) {
    specialized_spatial_constructor<T, DIM> init;
    this->is_nc_simplices_ = allow_nc_simplices;
    init(*this, mesh_type, fe_order, ref_level, tup_args);
  }

  /**
   * @brief Construct a new SpatialDiscretization object on a periodic mesh
   *        built with explicit translation vectors.
   *
   * @details Delegates the actual mesh construction to
   *          `specialized_spatial_constructor<T, DIM>`, which dispatches on
   *          `mesh_type` to build the underlying serial mesh, apply the given
   *          periodic `translations`, call `EnsureNCMesh(is_nc_simplices_)`
   *          before finalizing it, and convert it to a distributed
   *          `mfem::ParMesh`. See the non-periodic overload for the
   *          equivalent construction without translations.
   *
   * @tparam Args Types of the extra arguments forwarded to the specialized
   *             constructor (mesh dimensions, lengths, etc.), packed in
   *             `tup_args`.
   * @tparam T Finite element collection type.
   * @tparam DIM Spatial dimension.
   *
   * @param mesh_type Identifier of the mesh generation strategy to use (e.g.
   *                  "InlineSquareWithQuadrangles").
   * @param fe_order Finite element order used to build the finite element
   *                collection and space.
   * @param ref_level Number of uniform refinement levels applied at
   *                  construction time, before any AMR.
   * @param tup_args Tuple of extra parameters specific to `mesh_type` (e.g.
   *                number of elements and physical extents per direction).
   * @param translations Periodic translation vectors defining the
   *                     periodicity of the mesh.
   * @param allow_nc_simplices Whether nonconforming refinement is allowed on
   *                          simplex elements (see `is_nc_simplices()`). Must
   *                          be set correctly here since the underlying
   *                          `ParMesh` cannot be converted to nonconforming
   *                          after construction.
   */
  template <class... Args>
  explicit SpatialDiscretization(const std::string& mesh_type, const int& fe_order,
                                 const int& ref_level, std::tuple<Args...> tup_args,
                                 std::vector<mfem::Vector> translations,
                                 bool allow_nc_simplices = false) {
    specialized_spatial_constructor<T, DIM> init;
    this->is_nc_simplices_ = allow_nc_simplices;
    init(*this, mesh_type, fe_order, ref_level, tup_args, translations);
  }

  /**
   * @brief Construct a new SpatialDiscretization object that shares an
   *        already-existing parallel mesh, but builds its own independent
   *        finite element space on top of it.
   *
   * @details Use this constructor when several independent fields (e.g. `c`
   *          and `mu` in a Cahn-Hilliard system) must be discretized on the
   *          *same* mesh (so that mesh operations such as AMR stay consistent
   *          across all fields) while still owning distinct
   *          `mfem::ParFiniteElementSpace` objects. This instance does NOT
   *          take ownership of `existing_mesh`: it will never delete it, since
   *          that responsibility remains with the SpatialDiscretization that
   *          originally created the mesh. The caller must ensure that the
   *          owning instance outlives this one.
   *
   * @tparam T Finite element collection type.
   * @tparam DIM Spatial dimension.
   *
   * @param existing_mesh Pointer to an already-built parallel mesh, owned by
   *                      another SpatialDiscretization instance. Not deleted
   *                      by this instance's destructor.
   * @param fe_order Finite element order used to build the finite element
   *                collection and space attached to `existing_mesh`.
   * @param is_periodic_mesh Whether the shared mesh is periodic. Should match
   *                        the value used when `existing_mesh` was originally
   *                        built.
   */
  SpatialDiscretization(mfem::ParMesh* existing_mesh, const int& fe_order,
                        bool is_periodic_mesh = false) {
    this->mesh_ = existing_mesh;
    this->is_periodic_mesh_ = is_periodic_mesh;
    this->owns_mesh_ = false;  // ce mesh appartient à une autre instance -> ne pas le delete
    this->fe_order_ = fe_order;
    this->dimension_ = existing_mesh->Dimension();
    this->mesh_max_bdr_attributes_ =
        existing_mesh->bdr_attributes.Size() ? existing_mesh->bdr_attributes.Max() : 0;
    this->set_finite_element_space();
  }

  /**
   * @brief Build the parallel mesh by reading pre-split per-rank GMSH files.
   *
   * @details Expects `mesh_file` to be a naming prefix ending with a dot
   *          (e.g. "mesh."), with one file per MPI rank following MFEM's
   *          parallel file naming convention (via `mfem::MakeParFilename`,
   *          e.g. "mesh.000000", "mesh.000001", ...). The number of such
   *          files must exactly match the number of MPI processes used to
   *          run the simulation. Each rank reads only its own file and
   *          builds its local piece of the distributed `mfem::ParMesh`
   *          directly, without going through a serial `mfem::Mesh` first.
   *
   * @warning This path does not call `EnsureNCMesh()` before constructing
   *          the `ParMesh`, unlike the other mesh-construction paths in this
   *          class. As established, a conforming `ParMesh` cannot be
   *          converted to nonconforming after construction (in parallel), so
   *          nonconforming AMR will not be usable on a mesh built through
   *          this function unless the split files themselves already encode
   *          NC information, or this function is updated to match the other
   *          constructors.
   *
   * @tparam T Finite element collection type.
   * @tparam DIM Spatial dimension.
   *
   * @param mesh_file Naming prefix for the per-rank split mesh files; must
   *                  end with "." (e.g. "mesh.").
   * @return true  The mesh was successfully read and `mesh_` was built.
   * @return false `mesh_file` does not end with "." (naming convention not
   *               respected); no mesh is built and `mesh_` is left untouched.
   */
  bool GMSHReaderSplitFiles(const std::string mesh_file, bool allow_nc_simplices = false) {
    if (!mesh_file.ends_with(".")) return false;

    split_mesh_helper checker;

    int myid = mfem::Mpi::WorldRank();
    int mpi_size = mfem::Mpi::WorldSize();

    if (!checker(mpi_size, mesh_file)) {
      std::string msg = "SpatialDiscretization::SpatialDiscretization: " + mesh_file +
                        "* files is not correctly used. The number of MPI processes is different "
                        "of the number of files.";
      mfem::mfem_error(msg.c_str());
    }
    std::string fname(mfem::MakeParFilename(mesh_file, myid));
    std::ifstream ifs(fname);
    if (!ifs.good()) {
      std::string msg = "SpatialDiscretization::SpatialDiscretization: " + fname +
                        " doesn't exist. Please check your data.";
      mfem::mfem_error(msg.c_str());
    }

    this->is_nc_simplices_ = allow_nc_simplices;

    this->mesh_ = new mfem::ParMesh(MPI_COMM_WORLD, ifs);

    if (allow_nc_simplices) {
      MFEM_VERIFY(this->mesh_->Nonconforming(),
                  "GMSHReaderSplitFiles: the loaded mesh is conforming, but nonconforming "
                  "AMR was requested. The mesh files must be regenerated from a mesh that "
                  "had EnsureNCMesh() called before being split, since a parallel ParMesh "
                  "cannot be converted to nonconforming after construction.");
    }

    return true;
  }

  int fe_order_;
  int dimension_;
  mfem::ParMesh* mesh_;
  int mesh_max_bdr_attributes_;
  bool is_periodic_mesh_ = {false};

  void set_finite_element_space();

  void set_mesh_attributes_from_file();
  std::shared_ptr<mfem::AttributeSets> get_elem_attributes();
  std::shared_ptr<mfem::AttributeSets> get_bdr_attributes();

  std::shared_ptr<mfem::AttributeSets> elem_attr_sets_;
  std::shared_ptr<mfem::AttributeSets> bdr_attr_sets_;

  mfem::ParMesh* get_mesh();
  mfem::ParFiniteElementSpace* get_finite_element_space() const;
  mfem::ParFiniteElementSpace& get_ref_finite_element_space();

  std::size_t getSize() const;
  std::size_t get_max_bdr_attributes() const;
  int get_dimension() const;

  void apply_uniform_refinement(const int& level);

  bool is_periodic();
  bool is_nc_simplices() const;

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
  void operator()(SpatialDiscretization<T, 1>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, const std::string& file,
                  bool periodic_mesh) {
    a_my_class.fe_order_ = fe_order;
    a_my_class.dimension_ = 1;

    switch (Meshes::from(mesh_type)) {
      case Meshes::GMSH: {
        if (std::filesystem::exists(file)) {
          const char* mesh_file = file.c_str();
          mfem::Mesh tmp_mesh = mfem::Mesh::LoadFromFile(mesh_file, 1, 1);
          tmp_mesh.EnsureNCMesh(a_my_class.is_nc_simplices_);
          tmp_mesh.Finalize(true);

          a_my_class.mesh_ =
              mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh);  // definition of the parallel mesh
          tmp_mesh.Clear();
          break;
        } else if (a_my_class.GMSHReaderSplitFiles(file, a_my_class.is_nc_simplices_)) {
          SlothInfo::verbose(
              "SpatialDiscretization: enable GMSH reader from split files based on the pattern ",
              file);
          break;
        } else {
          std::string msg = "SpatialDiscretization::SpatialDiscretization: " + file +
                            " doesn't exist. Please check your data.";
          mfem::mfem_error(msg.c_str());
        }
      }
      default:
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only GMSH mesh type is "
            "allowed");
        break;
    }
    a_my_class.apply_uniform_refinement(ref_level);
    a_my_class.is_periodic_mesh_ = periodic_mesh;
    if (!a_my_class.is_periodic_mesh_) {
      a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
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
  void operator()(SpatialDiscretization<T, 1>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, std::tuple<Args...> tup_args) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);
    a_my_class.apply_uniform_refinement(ref_level);
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
  void operator()(SpatialDiscretization<T, 1>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, std::tuple<Args...> tup_args,
                  [[maybe_unused]] std::vector<mfem::Vector> translations) {
    this->build_periodic_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.apply_uniform_refinement(ref_level);

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
  void build_mesh(SpatialDiscretization<T, 1>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, std::tuple<Args...> tup_args) {
    a_my_class.fe_order_ = fe_order;
    auto tup_size = std::tuple_size<decltype(tup_args)>::value;
    a_my_class.dimension_ = 1;

    switch (Meshes::from(mesh_type)) {
      case Meshes::InlineLineWithSegments: {
        if (tup_size == 2) {
          const auto nx = std::get<0>(tup_args);
          const auto sx = std::get<1>(tup_args);
          mfem::Mesh tmp_mesh = mfem::Mesh::MakeCartesian1D(nx, sx);
          tmp_mesh.EnsureNCMesh(a_my_class.is_nc_simplices_);
          tmp_mesh.Finalize(true);
          a_my_class.mesh_ =
              new mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh);  // definition of the parallel mesh
          tmp_mesh.Clear();
        } else {
          mfem::mfem_error(
              "SpatialDiscretization::SpatialDiscretization: InlineLineWithSegments "
              "requires "
              "two "
              "argument, the number of segments");
        }
        break;
      }
      default:
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineLineWithSegments, "
            "InlineSquareWithTriangles, InlineSquareWithQuadrangles mesh types are allowed");
        break;
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
  }

  /**
   * @brief Build a periodic one-dimensional mesh by using MFEM inline methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void build_periodic_mesh(SpatialDiscretization<T, 1>& a_my_class, const std::string& mesh_type,
                           const int& fe_order, std::tuple<Args...> tup_args) {
    a_my_class.fe_order_ = fe_order;
    auto tup_size = std::tuple_size<decltype(tup_args)>::value;
    a_my_class.dimension_ = 1;

    switch (Meshes::from(mesh_type)) {
      case Meshes::InlineLineWithSegments: {
        if (tup_size == 2) {
          const auto nx = std::get<0>(tup_args);
          const auto sx = std::get<1>(tup_args);
          mfem::Mesh tmp_mesh = mfem::Mesh::MakeCartesian1D(nx, sx);

          // Based on mfem.org example
          // Create the vertex mapping. To begin, create the identity mapping.
          std::vector<int> periodicMap(tmp_mesh.GetNV());
          std::iota(periodicMap.begin(), periodicMap.end(), 0);
          // for (std::size_t i = 0; i < periodicMap.size(); ++i) {
          //   periodicMap[i] = i;
          // }
          // Modify the mapping so that the last vertex gets mapped to the first vertex.
          periodicMap.back() = 0;
          auto periodic_mesh = mfem::Mesh::MakePeriodic(tmp_mesh, periodicMap);
          tmp_mesh.Clear();
          mfem::Mesh tmp_mesh_periodic =
              mfem::Mesh(periodic_mesh, true);  // replace the input mesh with the periodic one
          tmp_mesh_periodic.EnsureNCMesh(a_my_class.is_nc_simplices_);
          tmp_mesh.Finalize(true);
          a_my_class.mesh_ = new mfem::ParMesh(
              MPI_COMM_WORLD, tmp_mesh_periodic);  // definition of the parallel mesh
          tmp_mesh_periodic.Clear();
        } else {
          mfem::mfem_error(
              "SpatialDiscretization::SpatialDiscretization: InlineLineWithSegments "
              "requires "
              "two "
              "argument, the number of segments");
        }
        break;
      }
      default:
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineLineWithSegments, "
            "InlineSquareWithTriangles, InlineSquareWithQuadrangles mesh types are allowed");
        break;
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
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
  void operator()(SpatialDiscretization<T, 2>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, const std::string& file,
                  bool periodic_mesh) {
    a_my_class.fe_order_ = fe_order;
    a_my_class.dimension_ = 2;

    switch (Meshes::from(mesh_type)) {
      case Meshes::GMSH: {
        if (std::filesystem::exists(file)) {
          const char* mesh_file = file.c_str();
          mfem::Mesh tmp_mesh(mesh_file, 1, 1);
          tmp_mesh.EnsureNCMesh(a_my_class.is_nc_simplices_);
          tmp_mesh.Finalize(true);

          a_my_class.mesh_ = new mfem::ParMesh(MPI_COMM_WORLD,
                                               tmp_mesh);  // definition of the parallel mesh
          tmp_mesh.Clear();
          break;
        } else if (a_my_class.GMSHReaderSplitFiles(file, a_my_class.is_nc_simplices_)) {
          SlothInfo::verbose(
              "SpatialDiscretization: enable GMSH reader from split files based on the pattern ",
              file);
          break;
        } else {
          std::string msg = "SpatialDiscretization::SpatialDiscretization: " + file +
                            " doesn't exist. Please check your data.";
          mfem::mfem_error(msg.c_str());
        }
      }
      default:
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only GMSH mesh type is "
            "allowed");
        break;
    }
    a_my_class.apply_uniform_refinement(ref_level);
    a_my_class.is_periodic_mesh_ = periodic_mesh;
    if (!a_my_class.is_periodic_mesh_) {
      a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
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
  void operator()(SpatialDiscretization<T, 2>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, std::tuple<Args...> tup_args) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.apply_uniform_refinement(ref_level);
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
  void operator()(SpatialDiscretization<T, 2>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, std::tuple<Args...> tup_args,
                  std::vector<mfem::Vector> translations) {
    this->build_periodic_mesh(a_my_class, mesh_type, fe_order, tup_args, translations);

    a_my_class.apply_uniform_refinement(ref_level);

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
  void build_mesh(SpatialDiscretization<T, 2>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, std::tuple<Args...> tup_args) {
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
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineSquareWithQuadrangles, InlineSquareWithTriangles mesh types are allowed");
        break;
    }
    if (tup_size == 4) {
      const auto nx = std::get<0>(tup_args);
      const auto ny = std::get<1>(tup_args);
      const auto sx = std::get<2>(tup_args);
      const auto sy = std::get<3>(tup_args);
      mfem::Mesh tmp_mesh = mfem::Mesh::MakeCartesian2D(nx, ny, element, false, sx, sy, false);
      tmp_mesh.EnsureNCMesh(a_my_class.is_nc_simplices_);
      tmp_mesh.Finalize(true);
      a_my_class.mesh_ = new mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh);  // definition of the parallel
                                                                       // mesh
      tmp_mesh.Clear();
    } else {
      std::string msg =
          "SpatialDiscretization::SpatialDiscretization: " + mesh_type +
          " requires 4 arguments, the number of nodes and the length along each direction";
      mfem::mfem_error(msg.c_str());
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
  }

  /**
   * @brief Build a periodic two-dimensional mesh by using MFEM inline methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void build_periodic_mesh(SpatialDiscretization<T, 2>& a_my_class, const std::string& mesh_type,
                           const int& fe_order, std::tuple<Args...> tup_args,
                           std::vector<mfem::Vector> translations) {
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

        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineSquareWithQuadrangles, InlineSquareWithTriangles mesh types are allowed");
        break;
    }
    if (tup_size == 4) {
      const auto nx = std::get<0>(tup_args);
      const auto ny = std::get<1>(tup_args);
      const auto sx = std::get<2>(tup_args);
      const auto sy = std::get<3>(tup_args);

      mfem::Mesh tmp_mesh = mfem::Mesh::MakeCartesian2D(nx, ny, element, false, sx, sy, false);

      const auto tol = 1.e-6;
      std::vector<int> periodicMap = tmp_mesh.CreatePeriodicVertexMapping(translations, tol);
      auto periodic_mesh = mfem::Mesh::MakePeriodic(tmp_mesh, periodicMap);
      tmp_mesh.Clear();
      mfem::Mesh tmp_mesh_periodic =
          mfem::Mesh(periodic_mesh, true);  // replace the input mesh with the periodic one
      tmp_mesh_periodic.EnsureNCMesh(a_my_class.is_nc_simplices_);
      tmp_mesh.Finalize(true);
      a_my_class.mesh_ =
          new mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh_periodic);  // definition of the parallel mesh
      tmp_mesh_periodic.Clear();
    } else {
      std::string msg =
          "SpatialDiscretization::SpatialDiscretization: " + mesh_type +
          " requires 4 arguments, the number of nodes and the length along each direction";
      mfem::mfem_error(msg.c_str());
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
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
  void operator()(SpatialDiscretization<T, 3>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, const std::string& file,
                  bool periodic_mesh) {
    a_my_class.fe_order_ = fe_order;
    a_my_class.dimension_ = 3;

    switch (Meshes::from(mesh_type)) {
      case Meshes::GMSH: {
        if (std::filesystem::exists(file)) {
          const char* mesh_file = file.c_str();
          mfem::Mesh tmp_mesh = mfem::Mesh::LoadFromFile(mesh_file, 1, 1);
          tmp_mesh.EnsureNCMesh(a_my_class.is_nc_simplices_);
          tmp_mesh.Finalize(true);
          a_my_class.mesh_ =
              new mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh);  // definition of the parallel mesh
          tmp_mesh.Clear();
          break;
        } else if (a_my_class.GMSHReaderSplitFiles(file, a_my_class.is_nc_simplices_)) {
          SlothInfo::verbose(
              "SpatialDiscretization: enable GMSH reader from split files based on the pattern ",
              file);
          break;
        } else {
          std::string msg = "SpatialDiscretization::SpatialDiscretization: " + file +
                            " doesn't exist. Please check your data.";
          mfem::mfem_error(msg.c_str());
        }
      }
      default:
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only GMSH mesh type is "
            "allowed");
        break;
    }
    a_my_class.apply_uniform_refinement(ref_level);
    a_my_class.is_periodic_mesh_ = periodic_mesh;
    if (!a_my_class.is_periodic_mesh_) {
      a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
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
  void operator()(SpatialDiscretization<T, 3>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, std::tuple<Args...> tup_args) {
    this->build_mesh(a_my_class, mesh_type, fe_order, tup_args);

    a_my_class.apply_uniform_refinement(ref_level);
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
  void operator()(SpatialDiscretization<T, 3>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, const int& ref_level, std::tuple<Args...> tup_args,
                  std::vector<mfem::Vector> translations) {
    this->build_periodic_mesh(a_my_class, mesh_type, fe_order, tup_args, translations);

    a_my_class.apply_uniform_refinement(ref_level);
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
  void build_mesh(SpatialDiscretization<T, 3>& a_my_class, const std::string& mesh_type,
                  const int& fe_order, std::tuple<Args...> tup_args) {
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
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineSquareWithTetraedres, InlineSquareWithHexaedres mesh types are allowed");
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
      mfem::Mesh tmp_mesh = mfem::Mesh::MakeCartesian3D(nx, ny, nz, element, sx, sy, sz);
      tmp_mesh.EnsureNCMesh(a_my_class.is_nc_simplices_);
      tmp_mesh.Finalize(true);
      a_my_class.mesh_ = new mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh);
      tmp_mesh.Clear();
    } else {
      std::string msg =
          "SpatialDiscretization::SpatialDiscretization: " + mesh_type +
          " requires six arguments, the number of nodes and the length along each direction";
      mfem::mfem_error(msg.c_str());
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
  }

  /**
   * @brief Build a periodic two-dimensional mesh by using MFEM inline methods
   *
   * @tparam Args
   * @param a_my_class
   * @param mesh_type
   * @param fe_order
   * @param tup_args
   */
  template <typename... Args>
  void build_periodic_mesh(SpatialDiscretization<T, 3>& a_my_class, const std::string& mesh_type,
                           const int& fe_order, std::tuple<Args...> tup_args,
                           std::vector<mfem::Vector> translations) {
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
        mfem::mfem_error(
            "SpatialDiscretization::SpatialDiscretization: here, only "
            "InlineSquareWithTetraedres, InlineSquareWithHexaedres mesh types are allowed");
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

      mfem::Mesh tmp_mesh = mfem::Mesh::MakeCartesian3D(nx, ny, nz, element, sx, sy, sz);

      const auto tol = 1.e-6;
      std::vector<int> periodicMap = tmp_mesh.CreatePeriodicVertexMapping(translations, tol);
      auto periodic_mesh = mfem::Mesh::MakePeriodic(tmp_mesh, periodicMap);
      tmp_mesh.Clear();
      mfem::Mesh tmp_mesh_periodic =
          mfem::Mesh(periodic_mesh, true);  // replace the input mesh with the periodic one
      tmp_mesh_periodic.EnsureNCMesh(a_my_class.is_nc_simplices_);
      tmp_mesh_periodic.Finalize(true);
      a_my_class.mesh_ =
          new mfem::ParMesh(MPI_COMM_WORLD, tmp_mesh_periodic);  // definition of the parallel mesh
      tmp_mesh_periodic.Clear();

    } else {
      std::string msg =
          "SpatialDiscretization::SpatialDiscretization: " + mesh_type +
          " requires six arguments, the number of nodes and the length along each direction";
      mfem::mfem_error(msg.c_str());
    }
    a_my_class.mesh_max_bdr_attributes_ = a_my_class.mesh_->bdr_attributes.Max();
  }
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief  return a pointer of Mesh
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::Mesh&
 * @note This method returns a mfem::Mesh instead of a mfem::ParMesh because it is called in
 *       the constructor of the PostProcessing objet.
 */
template <class T, int DIM>
mfem::ParMesh* SpatialDiscretization<T, DIM>::get_mesh() {
  return this->mesh_;
}

/**
 * @brief Set the FE_Collection, the FE_Space and associated size
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::ParFiniteElementSpace*
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::set_finite_element_space() {
  this->fecollection_ = new T(this->fe_order_, this->dimension_);
  this->fespace_ = new mfem::ParFiniteElementSpace(this->mesh_, this->fecollection_);
  // CCI
  this->size_ = this->fespace_->GetTrueVSize();
  int rank = mfem::Mpi::WorldRank();
  int taille = this->fespace_->GlobalTrueVSize();
  SlothInfo::debug("My Id = ", rank, " TrueVSize = ", size_, " and GlobalTrueVSize = ", taille);
  // CCI
}

/**
 * @brief Set the mesh attributes find in a mesh file
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param imesh
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::set_mesh_attributes_from_file() {
  this->elem_attr_sets_ = std::make_shared<mfem::AttributeSets>(this->mesh_->attribute_sets);
  this->bdr_attr_sets_ = std::make_shared<mfem::AttributeSets>(this->mesh_->bdr_attribute_sets);
}

/**
 * @brief Return the set of attributes associated with the elements
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::shared_ptr<mfem::AttributeSets>
 */
template <class T, int DIM>
std::shared_ptr<mfem::AttributeSets> SpatialDiscretization<T, DIM>::get_elem_attributes() {
  return this->elem_attr_sets_;
}

/**
 * @brief Return the set of attributes associated with the boundaries
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return std::shared_ptr<mfem::AttributeSets>
 */
template <class T, int DIM>
std::shared_ptr<mfem::AttributeSets> SpatialDiscretization<T, DIM>::get_bdr_attributes() {
  return this->bdr_attr_sets_;
}

/**
 * @brief return a pointer toward the finite element space
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::ParFiniteElementSpace*
 */
template <class T, int DIM>
mfem::ParFiniteElementSpace* SpatialDiscretization<T, DIM>::get_finite_element_space() const {
  return this->fespace_;
}

/**
 * @brief return a reference of the finite element space
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return mfem::ParFiniteElementSpace&
 */
template <class T, int DIM>
mfem::ParFiniteElementSpace& SpatialDiscretization<T, DIM>::get_ref_finite_element_space() {
  return *this->fespace_;
}

/**
 * @brief get the size of the Finite Element Space
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return int
 */
template <class T, int DIM>
std::size_t SpatialDiscretization<T, DIM>::getSize() const {
  return this->fespace_->GetTrueVSize();
}

/**
 * @brief get the maximum number of boundaries
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return int
 */
template <class T, int DIM>
std::size_t SpatialDiscretization<T, DIM>::get_max_bdr_attributes() const {
  return this->mesh_max_bdr_attributes_;
}

/**
 * @brief get the dimension of the problem
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return int
 */
template <class T, int DIM>
int SpatialDiscretization<T, DIM>::get_dimension() const {
  return this->dimension_;
}

/**
 * @brief Apply nb_ref uniform refinement
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @param nb_ref
 */
template <class T, int DIM>
void SpatialDiscretization<T, DIM>::apply_uniform_refinement(const int& nb_ref) {
  for (auto l = 0; l < nb_ref; l++) {
    this->mesh_->UniformRefinement();
  }
}

/**
 * @brief Return the flag to know whether the mesh is periodic or not
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 */
template <class T, int DIM>
bool SpatialDiscretization<T, DIM>::is_periodic() {
  return this->is_periodic_mesh_;
}

/**
 * @brief Return whether nonconforming refinement is allowed on simplex
 *        (triangle/tetrahedron) elements of this mesh.
 *
 * @details This flag only affects simplex meshes: it controls the
 *          `simplices_nonconforming` argument passed to `mfem::Mesh::EnsureNCMesh()`
 *          when the mesh is built. It has no effect for quadrilateral/hexahedral
 *          meshes, which natively support nonconforming refinement regardless
 *          of this flag.
 *
 * @tparam T Finite element collection type.
 * @tparam DIM Spatial dimension.
 * @return true  Nonconforming (hanging-node) refinement is allowed on simplices.
 * @return false Simplices are kept conforming (default MFEM behavior).
 */
template <class T, int DIM>
bool SpatialDiscretization<T, DIM>::is_nc_simplices() const {
  return this->is_nc_simplices_;
}

/**
 * @brief Destroy the Spatial Discretization< T>:: Spatial Discretization object
 *
 * @tparam T Finite element collection type.
 */
template <class T, int DIM>
SpatialDiscretization<T, DIM>::~SpatialDiscretization() {
  delete fespace_;
  if (this->owns_mesh_) {
    delete mesh_;
  }
  delete fecollection_;
}
