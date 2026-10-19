/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D Multigrains simulations from "Programming Phase-field modeling - Biner 2017"
 *        (monolithic scheme: all grains solved in one non-linear system)
 * @version 0.2
 * @date 2026-10-19
 *
 * @copyright Copyright (c) 2026
 *
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "./GrainsCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "voro++.hh"

using namespace voro;

/**
 * @brief Distance between two points of a periodic square domain [0, L)^2
 *        (minimum image convention: the closest periodic copy is used).
 */
double periodic_distance(double x1, double y1, double x2, double y2, double L) {
  double dx = std::abs(x2 - x1);
  double dy = std::abs(y2 - y1);
  dx = std::min(dx, L - dx);
  dy = std::min(dy, L - dy);
  return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief Map a coordinate into [0, L)
 */
double wrap(double x, double L) {
  double xw = std::fmod(x, L);
  if (xw < 0.) xw += L;
  return xw;
}

struct TestParameters {
  // Physics
  double control_mobility = 5.;
  double control_lambda = 0.1;

  // Mesh
  double control_length = 32.;
  int control_nb_elements = 32;
  int control_refinement_level = 0;

  // Initial polycrystal
  const char* control_seeds_file = "initial-grain.txt";

  // Non-linear solver
  double control_rel_tol = 1.e-10;
  double control_abs_tol = 1.e-14;
  int control_print_level = 1;

  // Time
  double control_initial_time = 0.0;
  double control_final_time = 0.2;
  double control_time_step = 1.e-1;

  // Post-processing
  int control_save_frequency = 2;
  int control_level_of_detail = 1;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.control_mobility, "-m", "--mobility", "AllenCahn mobility.");
  args.AddOption(&p.control_lambda, "-la", "--lambda", "Capillary coefficient.");
  args.AddOption(&p.control_length, "-L", "--length", "Size of the square domain.");
  args.AddOption(&p.control_nb_elements, "-n", "--nb_elements",
                 "Number of elements in each direction.");
  args.AddOption(&p.control_refinement_level, "-rl", "--refinement_level",
                 "Uniform refinement level of the mesh.");
  args.AddOption(&p.control_seeds_file, "-s", "--seeds", "Seed file: one grain per line (x y).");
  args.AddOption(&p.control_rel_tol, "-rtol", "--rel_tol", "Relative tolerance (LBFGS).");
  args.AddOption(&p.control_abs_tol, "-atol", "--abs_tol", "Absolute tolerance (LBFGS).");
  args.AddOption(&p.control_print_level, "-pl", "--print_level",
                 "Print level of the non-linear solver (-1: quiet).");
  args.AddOption(&p.control_initial_time, "-ti", "--initial_time",
                 "Initial time of the simulation.");
  args.AddOption(&p.control_final_time, "-t", "--final_time", "Final time of the simulation.");
  args.AddOption(&p.control_time_step, "-dt", "--time_step",
                 "Constant time-step of the simulation.");
  args.AddOption(&p.control_save_frequency, "-f", "--save_frequency", "Save frequency.");
  args.AddOption(&p.control_level_of_detail, "-lod", "--level_of_detail",
                 "Level of detail of the outputs.");

  args.Parse();

  if (!args.Good()) {
    if (mfem::Mpi::WorldRank() == 0) {
      args.PrintUsage(mfem::out);
    }
    std::exit(EXIT_FAILURE);
  }
  if (mfem::Mpi::WorldRank() == 0) args.PrintOptions(mfem::out);
}

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Quiet);
  //---------------------------------------
  // Initialize MPI
  //---------------------------------------
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------

  // ################ //
  //   Read options   //
  // ################ //
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);

  using namespace Sloth2D;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const auto refinement_level = p.control_refinement_level;
  const auto L = p.control_length;
  const auto NN = p.control_nb_elements;
  // Create translation vectors defining the periodicity
  mfem::Vector x_translation({L, 0.0});
  mfem::Vector y_translation({0.0, L});
  std::vector<mfem::Vector> translations = {x_translation, y_translation};
  SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithQuadrangles", 1,
                                                   refinement_level, std::make_tuple(NN, NN, L, L),
                                                   translations);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                     Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     Coefficients    //
  // ####################
  const double mob = p.control_mobility;
  const double lambda = p.control_lambda;

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GrainGradient());
  Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, GrainGw());
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficients coef_ac_grains(double_well_imp, capillary, mobility, grad_energy);

  // ####################
  //     variables     //
  // ####################
  const std::string filename = p.control_seeds_file;
  std::ifstream inputFile(filename);

  if (!inputFile) {
    std::cerr << "Error while opening " << filename << " for reading." << std::endl;
    return 1;
  }

  // Vector of tuples to store the data
  std::vector<std::tuple<double, double>> data;

  // Read and parse the file
  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    double col1, col2;
    if (!(iss >> col1 >> col2)) {
      std::cerr << "Error parsing line: " << line << std::endl;
      continue;
    }
    // Seeds are mapped into [0, L) so that the periodic distance is valid
    data.emplace_back(wrap(col1, L), wrap(col2, L));
  }

  // Close the file
  inputFile.close();

  // Number of grains given by the seed file
  const int ngrains = static_cast<int>(data.size());
  if (mfem::Mpi::WorldRank() == 0) std::cout << "Number of grains: " << ngrains << std::endl;

  container con(0., L,              // x bounds
                0., L,              // y bounds
                -0.5, 0.5,          // z bounds (set to zero for 2D case)
                NN, NN, 1,          // Number of grid subdivisions
                true, true, false,  // Periodic in x and y
                NN);

  for (int i = 0; i < ngrains; ++i) {
    con.put(i, std::get<0>(data[i]), std::get<1>(data[i]), 0.0);
  }

  std::vector<std::function<double(const mfem::Vector&, double)>> vect_user_func;
  vect_user_func.reserve(ngrains);
  for (int id_grain = 0; id_grain < ngrains; ++id_grain) {
    auto user_func = std::function<double(const mfem::Vector&, double)>(
        [id_grain, &con, L](const mfem::Vector& x, [[maybe_unused]] double time) {
          int closest_particle_id = -1;
          double min_distance = std::numeric_limits<double>::max();

          // Loop through all particles to find the one closest to the grid point
          c_loop_all vl(con);
          double xx, y, z;
          if (vl.start()) {
            do {
              // Get the current particle's position
              vl.pos(xx, y, z);

              // Periodic distance from the grid point to the current particle
              double dist = periodic_distance(x[0], x[1], xx, y, L);
              if (dist < min_distance) {
                min_distance = dist;
                closest_particle_id = vl.pid();  // Get the particle ID (Voronoi cell index)
              }

            } while (vl.inc());
          }

          const auto func = (closest_particle_id == id_grain) ? 1. : 0.;
          return func;
        });
    vect_user_func.emplace_back(user_func);
  }

  auto v0 = VAR(&spatial, bcs, "phi_0", Glossary::PhaseField, 2, vect_user_func[0]);
  auto var_grains = VARS(v0);
  for (int i = 1; i < ngrains; i++) {
    var_grains.add(Variable<FECollection, DIM>(&spatial, bcs, "phi_" + std::to_string(i),
                                               Glossary::PhaseField, 2, vect_user_func[i]));
  }

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string main_folder_path = "Saves";
  const auto level_of_detail = p.control_level_of_detail;
  const auto frequency = p.control_save_frequency;

  std::string calculation_path = "GrainsProblem";
  std::map<std::string, std::tuple<double, double>> map_threshold_integral;
  for (int i = 0; i < ngrains; ++i) {
    map_threshold_integral["phi_" + std::to_string(i)] = std::make_tuple(-1.1, 1.1);
  }
  bool enable_save_specialized_at_iter = true;
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail),
                 Parameter("integral_to_compute", map_threshold_integral),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto pst = PST(&spatial, p_pst);

  // ###########################################
  // ###########################################
  //      Problem                             //
  // ###########################################
  // ###########################################
  // Monolithic: all grains share the same spatial discretization and are solved together
  std::vector<SPA*> spatials(ngrains, &spatial);
  TransientOPE ope_ac_grains(spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");

  ope_ac_grains.overload_nl_solver(
      NLSolverType::LBFGS,
      Parameters(Parameter("description", "LBFGS solver"),
                 Parameter("print_level", p.control_print_level),
                 Parameter("rel_tol", p.control_rel_tol), Parameter("abs_tol", p.control_abs_tol)));

  std::vector<Coefficients> vcoeff(ngrains, coef_ac_grains);
  TransientPB problem_ac_grains(ope_ac_grains, var_grains, vcoeff, pst);
  // Coupling 1
  auto cc = Coupling("Multigrains ", problem_ac_grains);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto t_initial = p.control_initial_time;
  const auto t_final = p.control_final_time;
  const auto dt = p.control_time_step;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  time.solve();
  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();
  //---------------------------------------
  // MPI is finalized automatically by mfem::Mpi at exit, i.e. after the
  // destruction of the SLOTH/MFEM/hypre objects declared above.
  //---------------------------------------
  return 0;
}
