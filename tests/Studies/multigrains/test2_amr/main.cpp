/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 3D Multigrains (30) simulations from "Programming Phase-field modeling - Biner 2017" with
 * partitioned scheme (one problem per grain, the other grains being auxiliary variables) and
 * optional AMR
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
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "./GrainsCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "voro++.hh"

using namespace voro;

/**
 * @brief Distance between two points of a periodic cubic domain [0, L)^3
 *        (minimum image convention: the closest periodic copy is used).
 */
double periodic_distance(double x1, double y1, double z1, double x2, double y2, double z2,
                         double L) {
  double dx = std::abs(x2 - x1);
  double dy = std::abs(y2 - y1);
  double dz = std::abs(z2 - z1);
  dx = std::min(dx, L - dx);
  dy = std::min(dy, L - dy);
  dz = std::min(dz, L - dz);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
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
  int control_nb_elements = 8;  // coarse mesh when AMR is enabled
  int control_refinement_level = 0;

  // AMR
  bool control_amr = true;
  int control_amr_level = 1;
  double control_amr_max_elem_error = 1.e-3;

  // Initial polycrystal
  const char* control_seeds_file = "initial-grain.txt";

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
  args.AddOption(&p.control_length, "-L", "--length", "Size of the cubic domain.");
  args.AddOption(&p.control_nb_elements, "-n", "--nb_elements",
                 "Number of elements in each direction (coarse mesh if AMR).");
  args.AddOption(&p.control_refinement_level, "-rl", "--refinement_level",
                 "Uniform refinement level of the mesh.");
  args.AddOption(&p.control_amr, "-amr", "--amr", "-no-amr", "--no-amr", "Enables AMR.");
  args.AddOption(&p.control_amr_level, "-al", "--amr_level", "Maximum AMR level.");
  args.AddOption(&p.control_amr_max_elem_error, "-ae", "--amr_max_elem_error",
                 "Maximum element error for AMR.");
  args.AddOption(&p.control_seeds_file, "-s", "--seeds", "Seed file: one grain per line (x y z).");
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

  using namespace Sloth3D;

  // Number of grains (compile-time constant required by setCoupling<ngrains>)
  constexpr int ngrains = 30;

  // ###########################################
  // ###########################################
  //            Polycristaline                //
  // ###########################################
  // ###########################################
  const auto L = p.control_length;
  const auto NN = p.control_nb_elements;

  const std::string filename = p.control_seeds_file;
  std::ifstream inputFile(filename);

  if (!inputFile) {
    std::cerr << "Error while opening " << filename << " for reading." << std::endl;
    return 1;
  }

  // Vector of tuples to store the data
  std::vector<std::tuple<double, double, double>> data;

  // Read and parse the file
  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    double col1, col2, col3;
    if (!(iss >> col1 >> col2 >> col3)) {
      std::cerr << "Error parsing line: " << line << std::endl;
      continue;
    }
    // Seeds are mapped into [0, L) so that the periodic distance is valid
    data.emplace_back(wrap(col1, L), wrap(col2, L), wrap(col3, L));
  }

  // Close the file
  inputFile.close();

  if (static_cast<int>(data.size()) != ngrains) {
    if (mfem::Mpi::WorldRank() == 0) {
      std::cerr << "The seed file " << filename << " contains " << data.size()
                << " grains whereas this test is compiled for " << ngrains << " grains."
                << std::endl;
    }
    return 1;
  }

  container con(0., L,             // x bounds
                0., L,             // y bounds
                0., L,             // z bounds
                NN, NN, NN,        // Number of grid subdivisions
                true, true, true,  // Periodic in x, y and z
                NN);

  for (int i = 0; i < ngrains; ++i) {
    con.put(i, std::get<0>(data[i]), std::get<1>(data[i]), std::get<2>(data[i]));
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
              double dist = periodic_distance(x[0], x[1], x[2], xx, y, z, L);
              if (dist < min_distance) {
                min_distance = dist;
                // Get the particle ID (Voronoi cell index)
                closest_particle_id = vl.pid();
              }
            } while (vl.inc());
          }
          const auto func = (closest_particle_id == id_grain) ? 1. : 0.;
          return func;
        });
    vect_user_func.emplace_back(user_func);
  }

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  const auto refinement_level = p.control_refinement_level;
  // Create translation vectors defining the periodicity
  mfem::Vector x_translation({L, 0.0, 0.0});
  mfem::Vector y_translation({0.0, L, 0.0});
  mfem::Vector z_translation({0.0, 0.0, L});
  std::vector<mfem::Vector> translations = {x_translation, y_translation, z_translation};

  // One spatial discretization per grain (each grain has its own AMR mesh)
  SPAS spatials =
      setPeriodicSpatialDiscretization(ngrains, "InlineSquareWithHexaedres", 1, refinement_level,
                                       std::make_tuple(NN, NN, NN, L, L, L), translations, true);

  {  // scope: every object referring to the meshes is destroyed before deleting them
    // ##############################
    //     Boundary conditions     //
    // ##############################
    auto boundaries = {Boundary("rear", 0, "Periodic"),  Boundary("lower", 1, "Periodic"),
                       Boundary("right", 2, "Periodic"), Boundary("upper", 3, "Periodic"),
                       Boundary("left", 4, "Periodic"),  Boundary("front", 5, "Periodic")};
    auto bcs = setBoundaryConditions(ngrains, spatials, boundaries);

    // ###########################################
    // ###########################################
    //      Variables
    // ###########################################
    // ###########################################
    std::vector<VARS> vars_grains;
    vars_grains.reserve(ngrains);  // mandatory: auxiliary pointers below point into this vector
    for (int i = 0; i < ngrains; ++i) {
      auto vv = VAR(spatials[i], bcs[i], "phi_" + std::to_string(i), Glossary::PhaseField, 2,
                    AnalyticalFunctions<DIM>(vect_user_func[i]));
      vv.set_additional_information("phi");
      vars_grains.emplace_back(vv);
    }

    // ###########################################
    // ###########################################
    //      Post-processing
    // ###########################################
    // ###########################################
    const std::string main_folder_path = "Saves";
    const auto level_of_detail = p.control_level_of_detail;
    const auto frequency = p.control_save_frequency;
    const bool enable_save_specialized_at_iter = true;
    std::string calculation_path = "GrainsProblem";
    std::map<std::string, std::tuple<double, double>> map_threshold_integral;
    for (int i = 0; i < ngrains; ++i) {
      map_threshold_integral["phi_" + std::to_string(i)] = std::make_tuple(-1.1, 1.1);
    }
    auto p_pst =
        Parameters(Parameter("main_folder_path", main_folder_path),
                   Parameter("calculation_path", calculation_path),
                   Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
                   Parameter("integral_to_compute", map_threshold_integral),
                   Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
    auto pst_vec = PST(spatials[0], p_pst);

    // ###########################################
    // ###########################################
    //     Coefficients
    // ###########################################
    // ###########################################
    // Partitioned: the other grains are auxiliary variables -> explicit coupling
    const double mob = p.control_mobility;
    const double lambda = p.control_lambda;
    Coefficient grad_energy_i(Glossary::GradEnergy, Scheme::Implicit, GrainGradient());
    Coefficient double_well_exp_i(Glossary::FreeEnergy, Scheme::Implicit, ExplicitGrainGw());
    Coefficient capillary_i(Glossary::Capillary, lambda);
    Coefficient mobility_i(Glossary::Mobility, mob);
    Coefficient squares(Glossary::PhaseField, Scheme::Implicit, SumOfSquares());
    squares.set_name("Squares");
    std::vector<Coefficients> coeffs_i{
        Coefficients(double_well_exp_i, capillary_i, mobility_i, grad_energy_i)};

    // ###########################################
    // ###########################################
    //     Problems
    // ###########################################
    // ###########################################
    const int level_amr = p.control_amr_level;
    auto amr_params = Parameters(Parameter("max_elem_error", p.control_amr_max_elem_error),
                                 Parameter("amr_max_level", level_amr), Parameter("nc_limit", 0),
                                 Parameter("max_preref_cycles", level_amr));
    mfem::ConstantCoefficient amr_coef(1.0);  // must outlive the estimators

    // Operators are kept alive for the whole run: a TransientPB refers to its operator
    std::vector<std::unique_ptr<TransientOPE>> opes;
    opes.reserve(ngrains);
    std::vector<TransientPB> ac_pbs;
    ac_pbs.reserve(ngrains);

    for (int i = 0; i < ngrains; ++i) {
      // Operator
      opes.emplace_back(std::make_unique<TransientOPE>(
          std::vector<SPA*>{spatials[i]}, std::vector<std::string>{"AllenCahn"},
          TimeScheme::EulerImplicit, "TimeDerivative"));

      // Auxiliary variables: all the other grains
      std::vector<VARS*> aux_ptrs;
      aux_ptrs.reserve(ngrains - 1);
      for (int j = 0; j < ngrains; ++j) {
        if (j != i) aux_ptrs.push_back(&vars_grains[j]);
      }
      ac_pbs.emplace_back(*opes.back(), vars_grains[i], coeffs_i, pst_vec);
      ac_pbs.back().set_auxvariables(std::move(aux_ptrs));

      // AMR
      if (p.control_amr) {
        auto* amr_integ_i = new mfem::DiffusionIntegrator(amr_coef);
        auto* estimator_i = new SlothErrorEstimators(ErrorEstimatorType::KELLY, amr_integ_i);
        auto* amr_i =
            new MultiVariableMaxAMR<VARS>(*spatials[i]->get_mesh(), spatials[i]->is_nc_simplices());
        amr_i->SetCriteria(estimator_i, amr_params);
        ac_pbs.back().set_amr(amr_i);
      }
    }
    ac_pbs.back().set_vtk_coefficients({squares});

    // ###########################################
    // ###########################################
    //            Coupling                      //
    // ###########################################
    // ###########################################
    auto cc = setCoupling<ngrains>("Multigrains", ac_pbs);

    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const auto t_initial = p.control_initial_time;
    const auto t_final = p.control_final_time;
    const auto dt = p.control_time_step;
    auto time_params =
        Parameters(Parameter("initial_time", t_initial), Parameter("vtk_unified", true),
                   Parameter("final_time", t_final), Parameter("time_step", dt));
    auto time = TimeDiscretization(time_params, cc);

    time.solve();
  }

  deleteSpatialDiscretization(spatials);

  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();
  //---------------------------------------
  // MPI is finalized automatically by mfem::Mpi at exit, i.e. after the
  // destruction of the SLOTH/MFEM/hypre objects.
  //---------------------------------------
  return 0;
}