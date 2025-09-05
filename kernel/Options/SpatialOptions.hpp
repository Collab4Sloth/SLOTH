/**
 * @file SpatialOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for Spatial discretization
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <string>

#include "Utils/Utils.hpp"

#pragma once

///////////////////////////////////////////////////
//////// MESHES
///////////////////////////////////////////////////
struct Meshes {
  enum value {
    InlineLineWithSegments,
    InlineSquareWithTriangles,
    InlineSquareWithQuadrangles,
    InlineSquareWithTetraedres,
    InlineSquareWithHexaedres,
    GMSH
  };
  static value from(const std::string&);
};

Meshes::value Meshes::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<Meshes::value> m{
      {"InlineLineWithSegments", Meshes::InlineLineWithSegments},
      {"InlineSquareWithTriangles", Meshes::InlineSquareWithTriangles},
      {"InlineSquareWithQuadrangles", Meshes::InlineSquareWithQuadrangles},
      {"InlineSquareWithTetraedres", Meshes::InlineSquareWithTetraedres},
      {"InlineSquareWithHexaedres", Meshes::InlineSquareWithHexaedres},
      {"GMSH", Meshes::GMSH}};
  return m.find("Meshes", v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// Boundary conditions
///////////////////////////////////////////////////
struct BoundaryConditionType {
  enum value { Dirichlet, Neumann, Periodic, Robin };
  static value from(const std::string&);
};
BoundaryConditionType::value BoundaryConditionType::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<BoundaryConditionType::value> m{
      {"Dirichlet", BoundaryConditionType::Dirichlet},
      {"Neumann", BoundaryConditionType::Neumann},
      {"Robin", BoundaryConditionType::Robin},
      {"Periodic", BoundaryConditionType::Periodic}};
  return m.find("BoundaryConditionType", v);
}
