/**
 * @file PhaseFieldOptions.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-05-01
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "Utils/UtilsForDebug.hpp"

#pragma once

/*!
 * \brief Lambda expression used to manage exception
 * \param[in] b: boolean variable used in a conditional test
 * \param[in] method: string variable specifying the method concerned by the
 * exception
 * \param[in] msg: string variable specifying the error message written on the
 * screen
 */
static auto throw_if = [](const bool b, const std::string& method, const std::string& msg) {
  if (b) {
    SlothInfo::error(method, ": ", msg);
    mfem::mfem_error("Error message");
  }
};

static auto throw_iff = [](const bool b, const std::string& msg, int& errorLevel) {
  try {
    if (b) {
      throw msg;
    }
  } catch (std::string const& error) {
    std::cerr << error << std::endl;
    errorLevel++;
  }
};

static auto stringfindInVectorOfString = [](const std::vector<std::string> v, const std::string w) {
  return find(v.begin(), v.end(), w) != v.end();
};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Custom IterationKey for specialized map
 *
 */
struct IterationKey {
  std::pair<std::string, int> iter_;
  std::pair<std::string, double> time_step_;
  std::pair<std::string, double> time_;

  /**
   * @brief Construct a new Iteration Key object
   *
   * @param iter
   * @param time_step
   * @param time
   * @param s_iter
   * @param s_time_step
   * @param s_time
   */
  IterationKey(int iter, double time_step, double time, std::string s_iter = "Iter[-]",
               std::string s_time_step = "Dt[s]", std::string s_time = "Time[s]")
      : iter_(s_iter, iter), time_step_(s_time_step, time_step), time_(s_time, time) {}

  /**
   * @brief Comparison operator : mandatory for using IterationKey as a key in map
   *
   * @param user_key
   * @return true
   * @return false
   */
  bool operator<(const IterationKey& user_key) const {
    return std::tie(iter_, time_step_, time_) <
           std::tie(user_key.iter_, user_key.time_step_, user_key.time_);
  }
};
using SpecializedValue = std::pair<std::string, double>;

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
namespace PhaseFieldPrivate {
template <typename EType>
struct mmap : private std::vector<std::pair<const char* const, EType>> {
  using mpair = std::pair<const char* const, EType>;
  mmap(const std::initializer_list<mpair>&);
  EType find(const char* const, const std::string&);
};

/**
 * @brief Construct a new mmap<E Type>::mmap object
 *
 * @tparam EType
 * @param values
 */
template <typename EType>
mmap<EType>::mmap(const std::initializer_list<typename mmap::mpair>& values)
    : std::vector<typename mmap::mpair>(values) {}

/**
 * @brief
 *
 * @tparam EType
 * @param n
 * @param v
 * @return EType
 */
template <typename EType>
EType mmap<EType>::find(const char* const n, const std::string& v) {
  const auto pe = this->end();
  const auto p = std::find_if(this->begin(), pe, [&v](const mpair& e) {
    // Convert in uppercase for sake of generality
    std::string s1 = v;
    std::string s2 = e.first;
    transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
    transform(s2.begin(), s2.end(), s2.begin(), ::toupper);

    return s1 == s2;
  });
  if (p == pe) {
    std::string msg =
        "EnumNotFound::EnumNotFound : "
        "string '" +
        std::string(n) +
        "' "
        "is not a valid value for '" +
        v + "' keyword. See documentation";
    mfem::mfem_error(msg.c_str());
  }
  return p->second;
}
}  // namespace PhaseFieldPrivate

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct SourceTerm {
  enum value { Null, Sinusoide2D };
  static value from(const std::string&);
};

SourceTerm::value SourceTerm::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<SourceTerm::value> m{{"Null", SourceTerm::Null},
                                                      {"Sinusoide2D", SourceTerm::Sinusoide2D}};
  return m.find("SourceTerm", v);
}

enum class PhaseChange { Null, Constant, Calphad };
// struct ThermodynamicsPotentials {
//   enum value { W, F, H, X };
//   static value from(const std::string&);
// };

enum class ThermodynamicsPotentials { W, F, H, X, LOG };
enum class ThermodynamicsPotentialDiscretization { Implicit, Explicit, SemiImplicit };

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct AnalyticalFunctionsType {
  enum value { Heaviside, Sinusoide, Sinusoide2, HyperbolicTangent, Parabolic, Uniform };
  static value from(const std::string&);
};
AnalyticalFunctionsType::value AnalyticalFunctionsType::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<AnalyticalFunctionsType::value> m{
      {"Heaviside", AnalyticalFunctionsType::Heaviside},
      {"Sinusoide", AnalyticalFunctionsType::Sinusoide},
      {"Sinusoide2", AnalyticalFunctionsType::Sinusoide2},
      {"HyperbolicTangent", AnalyticalFunctionsType::HyperbolicTangent},
      {"Parabolic", AnalyticalFunctionsType::Parabolic},
      {"Uniform", AnalyticalFunctionsType::Uniform}};
  return m.find("AnalyticalFunctionsType", v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// Diffusion
///////////////////////////////////////////////////
enum class DiffusionCoefficients { Linear };
enum class CoefficientDiscretization {
  Implicit,
  Explicit,
};

///////////////////////////////////////////////////
//////// CONVERGENCE
///////////////////////////////////////////////////
struct ConvergenceType {
  enum value { RELATIVE_MAX, ABSOLUTE_MAX };
  static value from(const std::string&);
};
///////////////////////////////////////////////////
//////// SOLVERS
///////////////////////////////////////////////////
enum class NLSolverType { NEWTON };
enum class IterativeSolverType { BICGSTAB, GMRES, CG, MINRES };
enum class DirectSolverType { UMFPACK };
enum class HypreSolverType { HYPRE_PCG, HYPRE_GMRES, HYPRE_FGMRES };

///////////////////////////////////////////////////
//////// PRECONDITIONERS
///////////////////////////////////////////////////
enum class PreconditionerType { SMOOTHER, NO };
enum class HyprePreconditionerType {
  HYPRE_ILU,
  HYPRE_BOOMER_AMG,
  HYPRE_DIAG_SCALE,
  HYPRE_SMOOTHER,
  NO
};

///////////////////////////////////////////////////
//////// ODE SOLVER
///////////////////////////////////////////////////
struct TimeScheme {
  enum value { EulerImplicit, EulerExplicit, RungeKutta4 };
  static value from(const std::string&);
};
TimeScheme::value TimeScheme::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<TimeScheme::value> m{{"EulerImplicit", TimeScheme::EulerImplicit},
                                                      {"EulerExplicit", TimeScheme::EulerExplicit},
                                                      {"RungeKutta4", TimeScheme::RungeKutta4}};
  return m.find("TimeScheme", v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// Variables
///////////////////////////////////////////////////

struct VariableType {
  enum value { Primary, Auxiliary, Internal };
  static value from(const std::string&);
};
VariableType::value VariableType::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<VariableType::value> m{{"Primary", VariableType::Primary},
                                                        {"Auxiliary", VariableType::Auxiliary},
                                                        {"Inetrnal", VariableType::Internal}};
  return m.find("VariableType", v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// Problems
///////////////////////////////////////////////////

struct Problems {
  enum value { Diffusion, AllenCahn, Calphad };
  static value from(const std::string&);
};
Problems::value Problems::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<Problems::value> m{{"Diffusion", Problems::Diffusion},
                                                    {"AllenCahn", Problems::AllenCahn},
                                                    {"Calphad", Problems::Calphad}};
  return m.find("Problems", v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// Properties
///////////////////////////////////////////////////
enum class Property { Conductivity, Density, HeatCapacity, Diffusion, Mobility };

////////////////////////
//// Conductivity
///////////////////////
enum class Conductivity { Constant, Linear };

////////////////////////
//// Density
///////////////////////
enum class Density { Constant, Linear };

////////////////////////
//// HeatCapacity
///////////////////////
enum class HeatCapacity { Constant, Linear };

////////////////////////
//// Diffusion
///////////////////////
enum class Diffusion { Constant, Linear };

////////////////////////
//// Mobility
///////////////////////
enum class Mobility { Constant, Degenerated };
