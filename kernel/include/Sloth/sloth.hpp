/**
 * @file sloth.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief List of all header files contained in the kernel directory
 * @version 0.1
 * @date 2024-08-06
 *
 * Copyright CEA (c) 2024
 *
 */

#pragma once

#include "mfem.hpp"  // NOLINT

#ifdef SLOTH_USE_LIBTORCH
#include "include/Calphad/CalphadInformedNeuralNetwork.hpp"
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
