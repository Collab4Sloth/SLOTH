/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Covering test for FileCoefficient
 * @version 0.1
 * @date 2026-09-22
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
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

struct TestParameters {
  int tcase = 0;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.tcase, "-tc", "--test-case",
                 "0: FileCoefficient interpolation on a structured CSV grid (3D + time)");

  args.Parse();

  if (!args.Good()) {
    args.PrintUsage(mfem::out);
    std::exit(EXIT_FAILURE);
  }
  args.PrintOptions(mfem::out);
}

/**
 * @brief Affine reference field F(t, x, y, z) = t + 2x + 3y + 4z.
 */
double referenceF(double t, double x, double y, double z) {
  return t + 2.0 * x + 3.0 * y + 4.0 * z;
}

/**
 * @brief Write a CSV (t, x, y, z, F) on a regular structured grid with unit spacing,
 *        F being the affine reference field.
 *        Linear-in-time and trilinear-in-space interpolation is therefore
 *        exact everywhere in the domain, including between grid nodes.
 *
 * @param path Path of the CSV file to generate
 * @param nt,nx,ny,nz Number of points along each axis
 */
void generateAffineCsv(const std::filesystem::path& path, int nt, int nx, int ny, int nz) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Cannot create " + path.string());
  }

  out << "t,x,y,z,F\n";
  for (int it = 0; it < nt; ++it) {
    const double t = static_cast<double>(it);
    for (int ix = 0; ix < nx; ++ix) {
      const double x = static_cast<double>(ix);
      for (int iy = 0; iy < ny; ++iy) {
        const double y = static_cast<double>(iy);
        for (int iz = 0; iz < nz; ++iz) {
          const double z = static_cast<double>(iz);
          out << t << "," << x << "," << y << "," << z << "," << referenceF(t, x, y, z) << "\n";
        }
      }
    }
  }
}

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Verbose);

  constexpr double epsilon = 1.e-10;

  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);
  switch (p.tcase) {
    case 0: {
      // FileCoefficient interpolation on a structured CSV grid (3D + time)
      SlothInfo::print("Running test case 0: FileCoefficient CSV interpolation");

      const std::filesystem::path csvPath = "test_field_affine.csv";
      // 3x4x4x4 grid in (t, x, y, z), unit spacing
      generateAffineCsv(csvPath, /*nt=*/3, /*nx=*/4, /*ny=*/4, /*nz=*/4);

      try {
        Coefficient coeff{FileCoefficient(csvPath)};

        // Evaluate the coefficient at (t, x, y, z) and compare with the reference field
        auto checkAt = [&](const std::string& label, double t, double x, double y, double z) {
          coeff.set_time(t);
          const std::array<double, 3> pos{x, y, z};
          const std::span<const double> posSpan(pos.data(), pos.size());
          const double got =
              coeff.compute(posSpan, std::span<const double>{}, std::span<const double>{}, 3u);
          const double expected = referenceF(t, x, y, z);
          if (std::abs(got - expected) > epsilon) {
            throw std::runtime_error(label + ": expected " + std::to_string(expected) + ", got " +
                                     std::to_string(got));
          }
          SlothInfo::print(label, " OK: ", got, " expected ", expected);
        };

        // 1) Point on a grid node
        checkAt("Grid node", 1.0, 2.0, 1.0, 3.0);

        // 2) Point inside a cell, between time levels (linear + trilinear interpolation)
        checkAt("Interior point", 1.5, 1.3, 2.7, 0.4);

        // 3) Point on the domain boundary (last cell on each axis)
        checkAt("Boundary point", 2.0, 3.0, 3.0, 3.0);

        SlothInfo::print("Successful FileCoefficient interpolation test");
      } catch (const std::exception& e) {
        std::filesystem::remove(csvPath);
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
      }

      std::filesystem::remove(csvPath);
      break;
    }
    default:
      std::cerr << "Unknown test case: " << p.tcase << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
