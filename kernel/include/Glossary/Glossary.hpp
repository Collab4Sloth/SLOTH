/**
 * @file Glossary.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Glossary used to defined variables and properties in SLOTH
 * @version 0.1
 * @date 2025-11-07
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
#include <map>
#include <optional>
#include <string>

/**
 * @brief Enumeration of units used in SLOTH
 *
 */
enum class GlossaryUnit {
  Kelvin,
  Pascal,
  Watt,
  Second,
  Meter,
  Joules,
  Mole,
  SquareMeterPerSecond,
  MolePerCubicMeter,
  JoulesPerMole,
  JoulesPerMeter,
  JoulesPerSquareMeter,
  JoulesPerCubicMeter,
  CubicMeterPerJoulesPerSecond,
  MolesSquareMeterPerJoulesPerSecond,
  JoulesPerMeterPerSecondPerKelvin,
  JoulesPerMolePerKelvin,
  JoulesPerKelvin,
  None
};

/**
 * @brief Map used to print litteral unit from GlossaryUnit
 *
 */
static std::map<GlossaryUnit, std::string> GlossaryUnitMap = {
    {GlossaryUnit::Kelvin, "K"},
    {GlossaryUnit::Pascal, "Pa"},
    {GlossaryUnit::Second, "s"},
    {GlossaryUnit::Meter, "m"},
    {GlossaryUnit::SquareMeterPerSecond, "m2.s-1"},
    {GlossaryUnit::Mole, "mol"},
    {GlossaryUnit::Joules, "J"},
    {GlossaryUnit::JoulesPerMole, "J.mol-1"},
    {GlossaryUnit::JoulesPerMeter, "J.m-1"},
    {GlossaryUnit::JoulesPerSquareMeter, "J.m-2"},
    {GlossaryUnit::JoulesPerCubicMeter, "J.m-3"},
    {GlossaryUnit::MolePerCubicMeter, "mol.m-3"},
    {GlossaryUnit::JoulesPerMeterPerSecondPerKelvin, "J.s-1.m-1.K-1"},
    {GlossaryUnit::JoulesPerKelvin, "J.K-1"},
    {GlossaryUnit::JoulesPerMolePerKelvin, "J.mol-1.K-1"},
    {GlossaryUnit::CubicMeterPerJoulesPerSecond, "m3.J-1.s-1"},
    {GlossaryUnit::MolesSquareMeterPerJoulesPerSecond, "mol.m2.J-1.s-1"},
    {GlossaryUnit::None, "-"}};

/**
 * @brief Print litteral unit from GlossaryUnit
 *
 * @param unit GlossaryUnit
 * @return std::string Litteral version of the GlossaryUnit
 */
static std::string toString(GlossaryUnit unit) {
  auto it = GlossaryUnitMap.find(unit);
  if (it != GlossaryUnitMap.end()) {
    return it->second;
  }
  return "Unknown";
}

/**
 * @brief Enumeration of type of variables and properties used in SLOTH
 *
 */
enum class GlossaryType {
  Conductivity,
  HeatCapacity,
  Diffusivity,
  Mobility,
  SurfaceTension,
  Capillary,
  Temperature,
  Pressure,
  Concentration,
  Mole,
  MolarFraction,
  SiteFraction,
  ChemicalPotential,
  ThermodynamicPotential,
  FreeEnergy,
  GradientEnergy,
  PhaseFieldPotential,
  PhaseField,
  Thickness,
  System,
  Neumann,
  RobinA,
  RobinB,
  ExplicitTime
};

struct GlossaryQuantity {
  /**
   * @brief Construct a new GlossaryQuantity object
   *
   * @param t GlossaryType of the quantity
   * @param u GloassaryUnit of the quantity
   * @param d Description of the quantity
   */
  GlossaryQuantity(GlossaryType t, GlossaryUnit u, std::string d) : GlossaryQuantity(t, u, d, 0) {}
  GlossaryQuantity(GlossaryType t, GlossaryUnit u, std::string d, unsigned int id_qty)
      : type(t), unit(u), id(id_qty), description(d) {}

  GlossaryType type;
  GlossaryUnit unit;
  unsigned int id;
  std::string description;

  inline void setUnit(GlossaryUnit newUnit) { unit = newUnit; }
  inline void setId(unsigned int newId) { id = newId; }
};

/**
 * @brief Namespace for pre-defined GlossaryQuantity
 *
 */
namespace Glossary {
/**
 * @brief Quantity associated with phase-field variables
 *
 */
static const GlossaryQuantity PhaseField = GlossaryQuantity(
    GlossaryType::PhaseField, GlossaryUnit::None, "PhaseField variable (dimensionless)");
/**
 * @brief Quantity associated with a number of moles
 *
 */
static const GlossaryQuantity MoleNumber =
    GlossaryQuantity(GlossaryType::Mole, GlossaryUnit::Mole, "Mole number");
/**
 * @brief Quantity associated with molar fraction variables
 *
 */
static const GlossaryQuantity MoleFraction = GlossaryQuantity(
    GlossaryType::MolarFraction, GlossaryUnit::None, "Molar fraction variable (dimensionless)");
/**
 * @brief Quantity associated with site fraction variables
 *
 */
static const GlossaryQuantity SiteFraction = GlossaryQuantity(
    GlossaryType::SiteFraction, GlossaryUnit::None, "Site fraction variable (dimensionless)");

/**
 * @brief Quantity associated with phase-field mobility coefficients
 *
 */
static GlossaryQuantity Thickness =
    GlossaryQuantity(GlossaryType::Thickness, GlossaryUnit::Meter,
                     "Thickness of interface in " + toString(GlossaryUnit::Meter));

/**
 * @brief Quantity associated with phase-field mobility coefficients
 *
 */
static GlossaryQuantity Mobility = GlossaryQuantity(
    GlossaryType::Mobility, GlossaryUnit::CubicMeterPerJoulesPerSecond,
    "Mobility coefficient in " + toString(GlossaryUnit::CubicMeterPerJoulesPerSecond));

/**
 * @brief Quantity associated with inter-diffusion mobility coefficients
 *
 */
static GlossaryQuantity InterDiffusionMobility =
    GlossaryQuantity(GlossaryType::Mobility, GlossaryUnit::MolesSquareMeterPerJoulesPerSecond,
                     "Inter-diffusion mobility coefficient in " +
                         toString(GlossaryUnit::MolesSquareMeterPerJoulesPerSecond));

/**
 * @brief Quantity associated with mass diffusion coefficient
 *
 */
static GlossaryQuantity Diffusivity = GlossaryQuantity(
    GlossaryType::Diffusivity, GlossaryUnit::SquareMeterPerSecond,
    "Mass diffusion coefficient in " + toString(GlossaryUnit::SquareMeterPerSecond));

/**
 * @brief Quantity associated with chemical potential variables
 *
 */
static const GlossaryQuantity ChemicalPotential =
    GlossaryQuantity(GlossaryType::ChemicalPotential, GlossaryUnit::JoulesPerMole,
                     "Chemical potential variable in " + toString(GlossaryUnit::JoulesPerMole));

/**
 * @brief Quantity associated with concetration variables
 *
 */
static const GlossaryQuantity Concentration =
    GlossaryQuantity(GlossaryType::Concentration, GlossaryUnit::MolePerCubicMeter,
                     "Concentration variable in " + toString(GlossaryUnit::MolePerCubicMeter));
/**
 * @brief Quantity associated with surface tension property
 *
 */
static const GlossaryQuantity SurfaceTension =
    GlossaryQuantity(GlossaryType::SurfaceTension, GlossaryUnit::JoulesPerSquareMeter,
                     "Surface tension in " + toString(GlossaryUnit::JoulesPerSquareMeter));

/**
 * @brief Quantity associated with the capillary terms in phase-field equations
 *
 */
static const GlossaryQuantity Capillary =
    GlossaryQuantity(GlossaryType::Capillary, GlossaryUnit::JoulesPerMeter,
                     "Capillary coefficient in " + toString(GlossaryUnit::JoulesPerMeter));

/**
 * @brief Quantity associated with temperature in Kelvin
 *
 */
static const GlossaryQuantity Temperature =
    GlossaryQuantity(GlossaryType::Temperature, GlossaryUnit::Kelvin,
                     "Temperature in " + toString(GlossaryUnit::Kelvin));

/**
 * @brief Quantity associated with pressure in Pascal
 *
 */
static const GlossaryQuantity Pressure = GlossaryQuantity(
    GlossaryType::Pressure, GlossaryUnit::Pascal, "Pressure in " + toString(GlossaryUnit::Pascal));

/**
 * @brief Quantity associated with the molar Gibbs energy
 *
 */
static const GlossaryQuantity MolarGibbsEnergy =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::JoulesPerMole,
                     "Molar Gibbs free energy  in " + toString(GlossaryUnit::JoulesPerMole));
/**
 * @brief Quantity associated with the Gibbs energy in Joules
 */
static const GlossaryQuantity GibbsEnergy =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::Joules,
                     "Gibbs free energy in " + toString(GlossaryUnit::Joules));
/**
 * @brief Quantity associated with the molar enthalpy
 *
 */
static const GlossaryQuantity MolarEnthalpy =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::JoulesPerMole,
                     "Molar Enthalpy  in " + toString(GlossaryUnit::JoulesPerMole));
/**
 * @brief Quantity associated with the enthalpy in Joules
 *
 */
static const GlossaryQuantity Enthalpy =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::Joules,
                     "Enthalpy  in " + toString(GlossaryUnit::Joules));

/**
 * @brief Quantity associated with the interpolation functions used in phase-field equations
 *
 */
static const GlossaryQuantity InterpolationFunction =
    GlossaryQuantity(GlossaryType::PhaseFieldPotential, GlossaryUnit::None,
                     "Interpolation function (dimensionless)");

/**
 * @brief Quantity associated with the nucleus seed
 *
 */
static const GlossaryQuantity Nucleus =
    GlossaryQuantity(GlossaryType::PhaseField, GlossaryUnit::None, "Nucleus seed (dimensionless)");

/**
 * @brief Quantity associated with the driving force
 *
 */
static const GlossaryQuantity DrivingForce =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::Joules,
                     "Driving force in " + toString(GlossaryUnit::Joules));

/**
 * @brief Quantity associated with the free energy used in phase-field equations
 *
 */
static const GlossaryQuantity FreeEnergy = GlossaryQuantity(
    GlossaryType::FreeEnergy, GlossaryUnit::None, "Free energy function (dimensionless)");

// CCI ameliorer les commentaires et dimension
static const GlossaryQuantity GradEnergy = GlossaryQuantity(
    GlossaryType::GradientEnergy, GlossaryUnit::None, "Free energy function (dimensionless)");

/**
 * @brief Quantity associated with the thermal conductivity
 *
 */
static const GlossaryQuantity Conductivity = GlossaryQuantity(
    GlossaryType::Conductivity, GlossaryUnit::JoulesPerMeterPerSecondPerKelvin,
    "Themal conductivity in " + toString(GlossaryUnit::JoulesPerMeterPerSecondPerKelvin));

/**
 * @brief Quantity associated with the heat capacity
 *
 */
static const GlossaryQuantity Cp =
    GlossaryQuantity(GlossaryType::HeatCapacity, GlossaryUnit::JoulesPerKelvin,
                     "Heat capacity in " + toString(GlossaryUnit::JoulesPerKelvin));
/**
 * @brief Quantity associated with the molar heat capacity
 *
 */
static const GlossaryQuantity Cpm =
    GlossaryQuantity(GlossaryType::HeatCapacity, GlossaryUnit::JoulesPerMolePerKelvin,
                     "Heat capacity in " + toString(GlossaryUnit::JoulesPerMolePerKelvin));

/**
 * @brief Quantity associated with the Neumann boundary condition
 *
 */
static const GlossaryQuantity Neumann =
    GlossaryQuantity(GlossaryType::Neumann, GlossaryUnit::None, "Neumann boundary condition");

/**
 * @brief Quantity associated with the Robin boundary condition
 *
 */
static const GlossaryQuantity Robin_a =
    GlossaryQuantity(GlossaryType::RobinA, GlossaryUnit::None, "A-Robin boundary condition");

/**
 * @brief Quantity associated with the Robin boundary condition
 *
 */
static const GlossaryQuantity Robin_b =
    GlossaryQuantity(GlossaryType::RobinB, GlossaryUnit::None, "B-Robin boundary condition");

/**
 * @brief Quantity associated with the MPI rank
 *
 */
static const GlossaryQuantity MPI =
    GlossaryQuantity(GlossaryType::System, GlossaryUnit::None, "MPI rank variable (dimensionless)");

/**
 * @brief Quantity associated with the spatial coordinate variable
 *
 */
static const GlossaryQuantity Coordinate = GlossaryQuantity(
    GlossaryType::System, GlossaryUnit::None, "Spatial coordinate variable (dimensionless)");

/**
 * @brief Default quantity usefull integrators
 *
 */
static const GlossaryQuantity Default =
    GlossaryQuantity(GlossaryType::System, GlossaryUnit::None, "Default quantity (dimensionless)");

/**
 * @brief Quantity used for coefficients associated with time derivative with explicit solvers
 *
 */
static const GlossaryQuantity ExplicitTime_A = GlossaryQuantity(
    GlossaryType::ExplicitTime, GlossaryUnit::None,
    "Quantity used for coefficients associated with time derivative with explicit solvers", 0);

/**
 * @brief Quantity used for coefficients associated with time derivative with explicit solvers
 *
 */
static const GlossaryQuantity ExplicitTime_B = GlossaryQuantity(
    GlossaryType::ExplicitTime, GlossaryUnit::None,
    "Quantity used for coefficients associated with time derivative with explicit solvers", 1);

}  // namespace Glossary
