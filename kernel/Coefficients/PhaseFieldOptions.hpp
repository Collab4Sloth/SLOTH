/*
 * PhaseFieldOptions.hpp
 *
 *  Created on: 20 march 2023
 *      Author: ci230846
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

#ifndef THERMODYNAMICSOPTIONS_HH_
#define THERMODYNAMICSOPTIONS_HH_

// CALPHAD
struct ThermodynamicsInitialization {
  enum value { UserDefined, UserDefinedControlledAtmosphere, Cesar, Prodhel };
  static value from(const std::string&);
};

struct ThermodynamicSubModels {
  enum value {
    No,
    OnlyAllenCahn,
    ThermalDiffusion_U_O,
    ThermalDiffusion_U_PU_O,
    PhaseField_U_O,
    PhaseField_U_PU_O,
    PhaseField_Welland
  };
  static value from(const std::string&);
};

struct ThermodynamicsFGR {
  enum value { No, REP, RNR };
  static value from(const std::string&);
};

struct ThermodynamicsModel {
  enum value {
    No,
    OPENCALPHAD,
    OPENCALPHADwithOXITRAN,
    OPENCALPHADwithOXIRED,
    OPENCALPHADPHASEFIELD,
    OPENCALPHADPHASEFIELD_SPECIFICTIME
  };
  static value from(const std::string&);
};

struct ThermodynamicsUpdateQuantitiesMethod {
  enum value { MoleFractions, MoleNumbers, UnconservedMoleNumbers };
  static value from(const std::string&);
};

struct ThermodynamicsUnitForSettingConditions {
  enum value {
    MoleNumbers,
    MoleFractions,
    MoleNumbers2MolesFractions,
    OneMoleNumbers2MolesFractions
  };
  static value from(const std::string&);
};

// Global PhaseField system
enum class ThermodynamicsPhaseFieldModel {
  Calphad,
  AllenCahn,
  ThermalDiffusion,
  WellandThermalDiffusion
};
struct ThermodynamicsStabilizationScheme {
  // enum value { No, SU, Laplacian1, Laplacian2 };
  enum value { No, Laplacian1, Laplacian2 };
  static value from(const std::string&);
};
enum class ThermodynamicsTargetedProblem { Transient, Permanent };

enum class ThermodynamicsCoordinates { Cartesian, Cylindrical };
enum class ThermodynamicsBoundaryConditions {
  NeumannNeumann,
  NeumannDirichlet,
  DirichletDirichlet
};
enum class ThermodynamicsInitialConditions {
  Heaviside,
  HyperbolicProfile,
  Linear,
  Uniform,
  FromLoading
};

// HeatTransfer equation
enum class ThermodynamicsHeatCapacityType { Constant, FromLawWelland, FromThermodynamics };
enum class ThermodynamicsConductivityType { Constant, FromLawWelland };
enum class ThermodynamicsDensityType { Constant, FromLawWelland };
struct ThermodynamicsLinearPowerSourceTerm {
  enum value { Constant, TimeDependent, ModifiedBesselFunction };
  static value from(const std::string&);
};
// AllenCahn equation

struct ThermodynamicsMolarVolumeType {
  enum value { One, Constant };
  static value from(const std::string&);
};

enum class ThermodynamicsDoubleWellDiscretization { Implicit, Explicit, SemiImplicit };

struct ThermodynamicsPotentials {
  enum value { W, H, X };
  static value from(const std::string&);
};

struct ThermodynamicsAllenCahnMobility {
  enum value { Given, Logarithmic, LogarithmicMean };
  static value from(const std::string&);
};

// Oxygen ThermalDiffusion equation
enum class ThermodynamicsDiffusionType { Constant, FromMobilities, FromLawWelland, FromLawOXIRED };
enum class ThermodynamicsHeatTransportType {
  Constant,
  FromLawWelland,
  FromLawRamirez,
  FromLawKonarski
};

// ThermoPhysicalProperties
struct ThermodynamicsProperties {
  enum value {
    WPC,
    WTM,
    WME,
    WSD,
    WLD,
    OSD,
    WSHP,
    RSHP,
    KSHP,
    CMV,
    WLC,
    WSC,
    WSR,
    WLR,
    WLCP,
    WSCP
  };
  static value from(const std::string&);
};

// Root finding algorithm
struct ThermodynamicRootFindingAlgorithm {
  enum value { MU, DeltaMU };
  static value from(const std::string&);
};

// Thermodynamic system considered to calculate mobilities
struct ThermodynamicSystemsForMobilities {
  enum value { UO2, UPUO2 };
  static value from(const std::string&);
};

#endif /* __THERMODYNAMICSOPTIONS_HH_ */
