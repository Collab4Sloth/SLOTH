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

#pragma once

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

template <typename EType>
mmap<EType>::mmap(const std::initializer_list<typename mmap::mpair>& values)
    : std::vector<typename mmap::mpair>(values) {}

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
    std::runtime_error(
        "EnumNotFound::EnumNotFound : "
        "string '" +
        v +
        "' "
        "is not a valid value for '" +
        std::string(n) + "' keyword. See documentation");
  }
  return p->second;
}
}  // namespace PhaseFieldPrivate

///////////////////////////////////////////////////
//////// Thermodynamics
///////////////////////////////////////////////////

struct SourceTerm {
  enum value { Null, Sinusoide2D };
  static value from(const std::string&);
};

enum class Mobility { Constant, Degenerated };
enum class PhaseChange { Null, Constant, Calphad };
// struct ThermodynamicsPotentials {
//   enum value { W, F, H, X };
//   static value from(const std::string&);
// };

enum class ThermodynamicsPotentials { W, F, H, X };
enum class ThermodynamicsPotentialDiscretization { Implicit, Explicit, SemiImplicit };

struct AnalyticalFunctionsType {
  enum value { Heaviside, Sinusoide, Sinusoide2, HyperbolicTangent, Parabolic, Uniform };
  static value from(const std::string&);
};

///////////////////////////////////////////////////
//////// Diffusion
///////////////////////////////////////////////////
enum class DiffusionCoefficients { Linear };
enum class DiffusionCoefficientDiscretization { Implicit, Explicit, SemiImplicit };

///////////////////////////////////////////////////
//////// SOLVER
///////////////////////////////////////////////////
struct ConvergenceType {
  enum value { RELATIVE_MAX, ABSOLUTE_MAX };
  static value from(const std::string&);
};
///////////////////////////////////////////////////
//////// SOLVER
///////////////////////////////////////////////////
enum class NLSolverType { NEWTON };
enum class SolverType { BICGSTAB, GMRES, CG, MINRES, UMFPACK };
enum class IterativeSolverType { BICGSTAB, GMRES, CG, MINRES };
enum class DirectSolverType { UMFPACK };
enum class HypreSolverType { HYPRE_PCG, HYPRE_GMRES, HYPRE_FGMRES };

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

///////////////////////////////////////////////////
//////// Boundary conditions
///////////////////////////////////////////////////
struct BoundaryConditionType {
  enum value { Dirichlet, Neumann, Periodic, Robin };
  static value from(const std::string&);
};

///////////////////////////////////////////////////
//////// PROBLEMS
///////////////////////////////////////////////////
struct Problems {
  enum value { Diffusion, AllenCahn, CahnHilliard };
  static value from(const std::string&);
};

struct VariableType {
  enum value { Primary, Auxiliary };
  static value from(const std::string&);
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

//////////////////////////
// Conversion de l'option SourceTerm
SourceTerm::value SourceTerm::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<SourceTerm::value> m{{"Null", SourceTerm::Null},
                                                      {"Sinusoide2D", SourceTerm::Sinusoide2D}};
  return m.find("SourceTerm", v);
}

//////////////////////////
// Conversion de l'option BoundaryConditionType
BoundaryConditionType::value BoundaryConditionType::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<BoundaryConditionType::value> m{
      {"Dirichlet", BoundaryConditionType::Dirichlet},
      {"Neumann", BoundaryConditionType::Neumann},
      {"Robin", BoundaryConditionType::Robin},
      {"Periodic", BoundaryConditionType::Periodic}};
  return m.find("BoundaryConditionType", v);
}

//////////////////////////
// Conversion de l'option Problems
Problems::value Problems::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<Problems::value> m{{"Diffusion", Problems::Diffusion},
                                                    {"AllenCahn", Problems::AllenCahn},
                                                    {"CahnHilliard", Problems::CahnHilliard}};
  return m.find("Problems", v);
}

//////////////////////////
// Conversion de l'option  TimeScheme
TimeScheme::value TimeScheme::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<TimeScheme::value> m{{"EulerImplicit", TimeScheme::EulerImplicit},
                                                      {"EulerExplicit", TimeScheme::EulerExplicit},
                                                      {"RungeKutta4", TimeScheme::RungeKutta4}};
  return m.find("TimeScheme", v);
}

//////////////////////////
// Conversion de l'option  Meshes
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
//////////////////////////
// Conversion de l'option  AnalyticalFunctionsType
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

///////////////////////////////////////////////////
//////// Static methods
///////////////////////////////////////////////////
static auto throw_if = [](const bool b, const std::string& msg, int& errorLevel) {
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

///////////////////////////////////////////////////
////////
///////////////////////////////////////////////////
namespace utils {

void getInfo() {
  std::cout << " ################################\n";
  std::cout << " ## PF-MFEM small application    \n";
  std::cout << " ################################\n";
  std::cout << " - Available problems :\n";
  std::cout << "    * AnalyticalPoisson \n";
  std::cout << "    * Diffusion1D \n";
  std::cout << "    * AllenCahn1D \n";
  std::cout << "    * Diffusion2D \n";
  std::cout << "    * AllenCahn2D \n";
  std::cout << " - Available meshes :\n";
  std::cout << "    * InlineSquareWithTriangles \n";
  std::cout << "    * InlineSquareWithQuadrangles \n";
  std::cout << "    * InlineLineWithSegments \n";
  std::cout << "    * Camembert \n";
}

void checkOptions(const std::string& meshType, const std::string& targetedProblem) {
  int errorLevel = 0;
  // Meshes
  std::vector<std::string> availableMeshes{"InlineSquareWithTriangles", "InlineLineWithSegments",
                                           "InlineSquareWithQuadrangles", "Camembert"};
  throw_if(!stringfindInVectorOfString(availableMeshes, meshType),
           "\n>>> Error:  available meshes are InlineSquareWithTriangles, InlineLineWithSegments,  "
           "InlineSquareWithQuadrangles, Camembert",
           errorLevel);
  // Problems
  std::vector<std::string> availableProblems{"AnalyticalPoisson", "Diffusion1D", "Diffusion2D",
                                             "AllenCahn1D", "AllenCahn2D"};
  throw_if(!stringfindInVectorOfString(availableProblems, targetedProblem),
           "\n>>> Error:  available problems are AnalyticalPoisson, "
           "Diffusion1D, Diffusion2D, "
           "AllenCahn1D, AllenCahn2D",
           errorLevel);
  ////////////////////////////////
  // Check global level of errors
  ////////////////////////////////
  try {
    if (errorLevel > 0) {
      throw std::string("\n>>> Please check arguments: ") + std::to_string(errorLevel) +
          std::string(" errors detected");
    }
  } catch (std::string const& error) {
    std::cerr << error << std::endl;
    exit(1);
  }
}

std::tuple<Meshes::value, Problems> transformOptions(const std::string& meshType,
                                                     const std::string& targetedProblem) {
  Meshes::value m;
  Problems p;

  // if (targetedProblem == "AnalyticalPoisson") {
  //   p = Problems::AnalyticalPoisson;
  // } else if (targetedProblem == "Diffusion2D") {
  //   p = Problems::Diffusion2D;
  // } else if (targetedProblem == "Diffusion1D") {
  //   p = Problems::Diffusion1D;
  // } else if (targetedProblem == "AllenCahn2D") {
  //   p = Problems::AllenCahn2D;
  // } else if (targetedProblem == "AllenCahn1D") {
  //   p = Problems::AllenCahn1D;
  // } else {
  //   //
  // }
  if (meshType == "InlineLineWithSegments") {
    m = Meshes::InlineLineWithSegments;
  } else if (meshType == "InlineSquareWithTriangles") {
    m = Meshes::InlineSquareWithTriangles;
  } else if (meshType == "InlineSquareWithQuadrangles") {
    m = Meshes::InlineSquareWithQuadrangles;
  } else if (meshType == "Camembert") {
    m = Meshes::GMSH;
  } else {
    //
  }

  return std::make_tuple(m, p);
}
}  // namespace utils
