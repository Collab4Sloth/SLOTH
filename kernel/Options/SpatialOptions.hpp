/**
 * @file SpatialOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for Spatial discretization
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
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
