/**
 * @file SlothInit.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Initialisation of a SLOTH program: MPI, HYPRE and device
 * @version 0.1
 * @date 2026-11-16
 *
 * Copyright CEA (C) 2026
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

#include <string>

namespace Sloth {

/**
 * @brief Initialises a SLOTH program: MPI (mfem::Mpi::Init), HYPRE (mfem::Hypre::Init) and the
 * device on which the mfem::forall loops run (mfem::Device).
 *
 * The device is given by the command-line option --device=<config> (or --device <config>), or
 * by the environment variable SLOTH_DEVICE, with the MFEM syntax: "cpu" (default), "debug"
 * (memory managed as with a GPU, loops run on the host), "cuda", "hip"... The command-line option
 * takes priority over the environment variable; it is removed from argc / argv, which can then be
 * parsed by the program.
 *
 * With "cpu", no mfem::Device is created: the program runs as without device configuration.
 *
 * To be called once, at the beginning of main, instead of mfem::Mpi::Init and mfem::Hypre::Init.
 *
 * @param argc Number of command-line arguments (option --device removed)
 * @param argv Command-line arguments (option --device removed)
 */
void Init(int& argc, char**& argv);

/// Device configuration given to Init ("cpu" if no device is configured)
const std::string& DeviceConfiguration();

}  // namespace Sloth
