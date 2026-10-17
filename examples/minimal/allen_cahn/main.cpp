/**
 * @file main.cpp
 * @brief 1D Allen-Cahn: equilibrium of a planar interface (hyperbolic tangent profile)
 *
 * The phase field starts from the analytical equilibrium profile
 *     phi(x) = 1/2 + 1/2 tanh(2 (x - x0) / epsilon)
 * associated with the double-well potential omega * phi^2 (1 - phi)^2 and the gradient
 * energy lambda/2 |grad phi|^2, with lambda = 3 sigma epsilon / 2 and omega = 12 sigma / epsilon.
 * At equilibrium the profile must not move and the interfacial energy must be sigma:
 * see verify.py.
 *
 * The double-well potential is not hard-coded: the class DoubleWell is generated from
 * coefficients.json (expression and derivatives) when the program is built.
 */
#include <cmath>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "AllenCahnCoefficients.hpp"  // generated from coefficients.json at build time
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"  // Sloth1D/2D/3D namespaces (SPA, VAR, TransientPB, ...)

int main(int argc, char* argv[]) {
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  Profiling::getInstance().enable();

  using namespace Sloth1D;

  // ##############################
  //   Mesh and boundary conditions
  // ##############################
  const double L = 1.e-3;  // length of the domain [m]
  const int nb_elements = 100;
  SPA spatial("InlineLineWithSegments", 1, 0, std::make_tuple(nb_elements, L));

  auto boundaries = {Boundary("left", 0, "Neumann"), Boundary("right", 1, "Neumann")};
  auto bcs = BCS(&spatial, boundaries);

  // ##############################
  //   Physical parameters
  // ##############################
  const double epsilon = 1.e-4;  // interface thickness [m]
  const double sigma = 6.e-2;    // interfacial energy [J.m-2]
  const double mob = 1.e-5;      // mobility
  const double lambda = 3. * sigma * epsilon / 2.;
  const double omega = 12. * sigma / epsilon;

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  // Double-well potential omega * phi^2 (1 - phi)^2, generated from coefficients.json
  Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, DoubleWell(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficients coefs(double_well, capillary, mobility, grad_energy);

  // ##############################
  //   Phase field: initial condition = analytical equilibrium profile
  // ##############################
  const double center_x = 0.;
  const double a_x = 1.;
  const double x0 = L / 2.;  // position of the interface
  auto equilibrium_profile = std::function<double(const mfem::Vector&, double)>(
      [x0, epsilon](const mfem::Vector& x, [[maybe_unused]] double time) {
        return 0.5 + 0.5 * std::tanh(2. * (x[0] - x0) / epsilon);
      });
  auto initial_condition = AnalyticalFunctions<DIM>(equilibrium_profile);
  // Used by SLOTH to compute the error with respect to the analytical solution
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::from("HyperbolicTangent"), center_x, a_x, epsilon, x0);
  auto vars = VARS(
      VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition, analytical_solution));

  // ##############################
  //   Problem and post-processing (Saves/Problem1)
  // ##############################
  const std::string main_folder_path = "Saves";
  const std::string calculation_path = "Problem1";
  auto pst = PST(&spatial, Parameters(Parameter("main_folder_path", main_folder_path),
                                      Parameter("calculation_path", calculation_path),
                                      Parameter("frequency", 1), Parameter("level_of_detail", 1)));
  TransientOPE oper({&spatial}, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  TransientPB problem(oper, vars, {coefs}, pst);
  auto cc = Coupling("Allen-Cahn", problem);

  // ##############################
  //   Time integration
  // ##############################
  auto time_params = Parameters(Parameter("initial_time", 0.), Parameter("final_time", 10.),
                                Parameter("time_step", 0.1));
  auto time = TimeDiscretization(time_params, cc);
  time.solve();

  Profiling::getInstance().print();
  MPI_Finalize();
  return 0;
}
