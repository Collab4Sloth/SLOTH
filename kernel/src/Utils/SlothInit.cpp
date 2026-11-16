/**
 * @file SlothInit.cpp
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
#include "Utils/SlothInit.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

namespace Sloth {
namespace {
/// Device configuration given to Init
std::string& device_configuration() {
  static std::string configuration = "cpu";
  return configuration;
}
}  // namespace

void Init(int& argc, char**& argv) {
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();

  // Environment variable, then command-line option (removed from argc / argv)
  std::string configuration = "cpu";
  if (const char* env = std::getenv("SLOTH_DEVICE"); env != nullptr && *env != '\0') {
    configuration = env;
  }
  constexpr std::string_view option = "--device";
  int kept = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if (arg.starts_with(option) && arg.size() > option.size() && arg[option.size()] == '=') {
      configuration = arg.substr(option.size() + 1);
    } else if (arg == option && i + 1 < argc) {
      configuration = argv[++i];
    } else {
      argv[kept++] = argv[i];
    }
  }
  argc = kept;
  argv[argc] = nullptr;

  device_configuration() = configuration;
  if (configuration == "cpu") {
    return;
  }
  // Configured once, kept until the end of the program
  static const std::unique_ptr<mfem::Device> device = std::make_unique<mfem::Device>(configuration);
  if (mfem::Mpi::Root()) {
    device->Print();
  }
}

const std::string& DeviceConfiguration() { return device_configuration(); }

}  // namespace Sloth
