

/*
 * Copyright © CEA 2023
 *
 * pf_constants.hpp
 *
 *  Created on: 20 March 2023
 *      Author: ci230846
 */
#pragma once

/**
 * @brief Default constant used by Newton Solver
 *
 */
namespace NewtonDefaultConstant {
const auto iter_max = 100;
const auto abs_tol = 1.e-13;
const auto rel_tol = 1.e-13;
const bool iterative_mode = false;
const auto print_level = 1;
}  // namespace NewtonDefaultConstant

/**
 * @brief Default constant used by BICGSTAB Solver
 *
 */
namespace BICGSTABDefaultConstant {
const auto iter_max = 100;
const auto abs_tol = 1.e-13;
const auto rel_tol = 1.e-13;
const bool iterative_mode = false;
const auto print_level = 1;
}  // namespace BICGSTABDefaultConstant

/**
 * @brief Default constant used by CG Solver
 *
 */
namespace CGDefaultConstant {
const auto iter_max = 100;
const auto abs_tol = 1.e-13;
const auto rel_tol = 1.e-13;
const bool iterative_mode = false;
const auto print_level = 1;
}  // namespace CGDefaultConstant

/**
 * @brief Default constant used by GMRES Solver
 *
 */
namespace GMRESDefaultConstant {
const auto iter_max = 100;
const auto abs_tol = 1.e-13;
const auto rel_tol = 1.e-13;
const bool iterative_mode = false;
const auto print_level = 1;
}  // namespace GMRESDefaultConstant

/**
 * @brief Default constant used by Mass Solver
 *
 */
namespace MassDefaultConstant {
const auto iter_max = 30;
const auto abs_tol = 1.e-15;
const auto rel_tol = 1.e-15;
const bool iterative_mode = false;
const auto print_level = -1;
}  // namespace MassDefaultConstant
