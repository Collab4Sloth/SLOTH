/**
 * @file sloth.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief List of SLOTH header files
 * @version 0.1
 * @date 2026-03-21
 *
 * @copyright CEA (C) 2026
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

#include "mfem.hpp"  // NOLINT

#ifdef SLOTH_USE_LIBTORCH
#include "Calphad/CalphadInformedNeuralNetwork.hpp"
#endif
#ifdef SLOTH_USE_EXPRTK
#include "Coefficients/ExprTkSloth.hpp"
#endif
#include "Calphad/AnalyticalIdealSolution.hpp"
#include "Calphad/CalphadBase.hpp"
#include "Calphad/KKS.hpp"
#include "Coefficients/ListCoefficients.hpp"
#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Couplings/Coupling.hpp"
#include "Operators/ListOperators.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/ListProblems.hpp"
#include "Property/InterDiffusionCoefficient.hpp"
#include "Property/PropertyBase.hpp"
#include "Spatial/Spatial.hpp"
#include "Time/Time.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
