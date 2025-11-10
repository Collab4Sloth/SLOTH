/**
 * @file Glossary.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Glossary used to defined variables and properties in SLOTH
 * @version 0.1
 * @date 2025-11-07
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
#include <unordered_map>
#pragma once

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
  MolePerCubicMeter,
  JoulesPerMole,
  JoulesPerMeter,
  JoulesPerSquareMeter,
  JoulesPerCubicMeter,
  CubicMeterPerJoulesPerSecond,
  MolesSquareMeterPerJoulesPerSecond,
  None
};

/**
 * @brief Map used to print litteral unit from GlossaryUnit
 *
 */
static std::unordered_map<GlossaryUnit, std::string> GlossaryUnitMap = {
    {GlossaryUnit::Kelvin, "K"},
    {GlossaryUnit::Pascal, "Pa"},
    {GlossaryUnit::Watt, "W"},
    {GlossaryUnit::Second, "s"},
    {GlossaryUnit::Meter, "m"},
    {GlossaryUnit::Mole, "mol"},
    {GlossaryUnit::Joules, "J"},
    {GlossaryUnit::JoulesPerMole, "J.mol-1"},
    {GlossaryUnit::JoulesPerMeter, "J.m-1"},
    {GlossaryUnit::JoulesPerSquareMeter, "J.m-2"},
    {GlossaryUnit::JoulesPerCubicMeter, "J.m-3"},
    {GlossaryUnit::MolePerCubicMeter, "mol.m-3"},
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
  Density,
  SurfaceTension,
  Capillary,
  Temperature,
  Pressure,
  Concentration,
  MolarFraction,
  SiteFraction,
  ChemicalPotential,
  ThermodynamicPotential,
  PhaseFieldPotential,
  PhaseField,
  System
};

/**
 * @brief Structure used to define a GlossaryQuantity
 *
 */
struct GlossaryQuantity {
  /**
   * @brief Construct a new GlossaryQuantity object
   *
   * @param t GlossaryType of the quantity
   * @param u GloassaryUnit of the quantity
   * @param d Description of the quantity
   */
  GlossaryQuantity(GlossaryType t, GlossaryUnit u, std::string d)
      : type(t), unit(u), description(d) {}

  GlossaryType type;
  GlossaryUnit unit;
  std::string description;
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
static const GlossaryQuantity Phi = GlossaryQuantity(GlossaryType::PhaseField, GlossaryUnit::None,
                                                     "PhaseField variable (dimensionless)");
/**
 * @brief Quantity associated with molar fraction variables
 *
 */
static const GlossaryQuantity X = GlossaryQuantity(GlossaryType::MolarFraction, GlossaryUnit::None,
                                                   "Molar fraction variable (dimensionless)");

/**
 * @brief Quantity associated with Allen-Cahn mobility coefficients
 *
 */
static const GlossaryQuantity MobPhi = GlossaryQuantity(
    GlossaryType::Mobility, GlossaryUnit::CubicMeterPerJoulesPerSecond,
    "Allen-Cahn Mobility coefficient in " + toString(GlossaryUnit::CubicMeterPerJoulesPerSecond));

/**
 * @brief Quantity associated with inter-diffusion mobility coefficients
 *
 */
static const GlossaryQuantity Mob =
    GlossaryQuantity(GlossaryType::Mobility, GlossaryUnit::MolesSquareMeterPerJoulesPerSecond,
                     "Inter-diffusion mobility coefficient in " +
                         toString(GlossaryUnit::MolesSquareMeterPerJoulesPerSecond));

/**
 * @brief Quantity associated with chemical potential variables
 *
 */
static const GlossaryQuantity Mu =
    GlossaryQuantity(GlossaryType::ChemicalPotential, GlossaryUnit::JoulesPerMole,
                     "Chemical potential variable in " + toString(GlossaryUnit::JoulesPerMole));

/**
 * @brief Quantity associated with concetration variables
 *
 */
static const GlossaryQuantity C =
    GlossaryQuantity(GlossaryType::Concentration, GlossaryUnit::MolePerCubicMeter,
                     "PhaseField variable in " + toString(GlossaryUnit::MolePerCubicMeter));
/**
 * @brief Quantity associated with surface tension property
 *
 */
static const GlossaryQuantity Sigma =
    GlossaryQuantity(GlossaryType::SurfaceTension, GlossaryUnit::JoulesPerSquareMeter,
                     "Surface tension in " + toString(GlossaryUnit::JoulesPerSquareMeter));

/**
 * @brief Quantity associated with the capillary terms in phase-field equations
 *
 */
static const GlossaryQuantity Kappa =
    GlossaryQuantity(GlossaryType::Capillary, GlossaryUnit::JoulesPerMeter,
                     "Capillary coefficient in " + toString(GlossaryUnit::JoulesPerMeter));

/**
 * @brief Quantity associated with temperature in Kelvin
 *
 */
static const GlossaryQuantity Tk =
    GlossaryQuantity(GlossaryType::Temperature, GlossaryUnit::Kelvin,
                     "Temperature in " + toString(GlossaryUnit::Kelvin));

/**
 * @brief Quantity associated with pressure in Pascal
 *
 */
static const GlossaryQuantity P = GlossaryQuantity(GlossaryType::Temperature, GlossaryUnit::Pascal,
                                                   "Pressure in " + toString(GlossaryUnit::Pascal));

/**
 * @brief Quantity associated with the molar Gibbs energy
 *
 */
static const GlossaryQuantity Gm =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::JoulesPerMole,
                     "Molar Gibbs free energy  in " + toString(GlossaryUnit::JoulesPerMole));
/**
 * @brief Quantity associated with the Gibbs energy in Joules
 */
static const GlossaryQuantity G =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::Joules,
                     "Gibbs free energy in " + toString(GlossaryUnit::Joules));
/**
 * @brief Quantity associated with the molar enthalpy
 *
 */
static const GlossaryQuantity Hm =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::JoulesPerMole,
                     "Molar Enthalpy  in " + toString(GlossaryUnit::JoulesPerMole));
/**
 * @brief Quantity associated with the enthalpy in Joules
 *
 */
static const GlossaryQuantity H =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::Joules,
                     "Enthalpy  in " + toString(GlossaryUnit::Joules));

/**
 * @brief Quantity associated with the interpolation functions used in phase-field equations
 *
 */
static const GlossaryQuantity Pint =
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
static const GlossaryQuantity Dgm =
    GlossaryQuantity(GlossaryType::ThermodynamicPotential, GlossaryUnit::Joules,
                     "Driving force in " + toString(GlossaryUnit::Joules));

/**
 * @brief Quantity associated with the free energy used in phase-field equations
 *
 */
static const GlossaryQuantity F = GlossaryQuantity(
    GlossaryType::PhaseFieldPotential, GlossaryUnit::None, "Free energy function (dimensionless)");

/**
 * @brief Quantity associated with the MPI rank
 *
 */
static const GlossaryQuantity MPI =
    GlossaryQuantity(GlossaryType::System, GlossaryUnit::None, "MPI rank variable (dimensionless)");

}  // namespace Glossary
