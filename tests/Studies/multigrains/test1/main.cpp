/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Multigrains simulations from "Programming Phase-field modeling - Biner 2017"
 * @version 0.1
 * @date 2026-05-08
 *
 * @copyright Copyright (c) 2026
 *
 */
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include "./GrainsCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "voro++.hh"

using namespace voro;

double distance(double x1, double y1, double x2, double y2) {
  return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
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
  //
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  /////////////////////////
  using OPE = TransientOperator<FECollection, DIM>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;  // 2;
  auto L = 32.;
  auto NN = 32;
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

  const double& mob(5.);
  const double& lambda(0.1);

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GrainGradient());
  Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, GrainGw());
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficients coef_ac_grains(double_well_imp, capillary, mobility, grad_energy);
  // ####################
  //     variables     //
  // ####################

  const std::string& filename = "initial-grain.txt";
  std::ifstream inputFile(filename);

  if (!inputFile) {
    std::cerr << "Error while opening " << filename << " for reading." << std::endl;
    return 1;
  }

  // Vector of tuples to store the data
  std::vector<std::tuple<double, double>> data;
  std::vector<double> dmin;

  // Read and parse the file
  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    double col1, col2;
    if (!(iss >> col1 >> col2)) {
      std::cerr << "Error parsing line: " << line << std::endl;
      continue;
    }
    data.emplace_back(col1, col2);
  }

  // Close the file
  inputFile.close();

  container con(0., L,              // x bounds
                0., L,              // y bounds
                -0.5, 0.5,          // z bounds (set to zero for 2D case)
                NN, NN, 1,          // Number of grid subdivisions
                true, true, false,  // No periodic boundaries
                NN);

  for (int i = 0; i < data.size(); ++i) {
    con.put(i, std::get<0>(data[i]), std::get<1>(data[i]), 0.0);
  }

  std::vector<std::function<double(const mfem::Vector&, double)>> vect_user_func;
  for (int id_grain = 0; id_grain < static_cast<int>(data.size()); ++id_grain) {
    auto user_func = std::function<double(const mfem::Vector&, double)>(
        [id_grain, &con](const mfem::Vector& x, [[maybe_unused]] double time) {
          int closest_particle_id = -1;
          double min_distance = std::numeric_limits<double>::max();

          // Loop through all particles to find the one closest to the grid point
          c_loop_all vl(con);
          double xx, y, z;
          if (vl.start()) {
            do {
              // Get the current particle's position
              vl.pos(xx, y, z);

              // Calculate the distance from the grid point to the current particle
              double dist = distance(x[0], x[1], xx, y);
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

  const int ngrains = 30;
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
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 2;

  std::string calculation_path = "GrainsProblem";
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
      {"phi_0", {-1.1, 1.1}},  {"phi_1", {-1.1, 1.1}},  {"phi_2", {-1.1, 1.1}},
      {"phi_3", {-1.1, 1.1}},  {"phi_4", {-1.1, 1.1}},  {"phi_5", {-1.1, 1.1}},
      {"phi_6", {-1.1, 1.1}},  {"phi_7", {-1.1, 1.1}},  {"phi_8", {-1.1, 1.1}},
      {"phi_9", {-1.1, 1.1}},  {"phi_10", {-1.1, 1.1}}, {"phi_11", {-1.1, 1.1}},
      {"phi_12", {-1.1, 1.1}}, {"phi_13", {-1.1, 1.1}}, {"phi_14", {-1.1, 1.1}},
      {"phi_15", {-1.1, 1.1}}, {"phi_16", {-1.1, 1.1}}, {"phi_17", {-1.1, 1.1}},
      {"phi_18", {-1.1, 1.1}}, {"phi_19", {-1.1, 1.1}}, {"phi_20", {-1.1, 1.1}},
      {"phi_21", {-1.1, 1.1}}, {"phi_22", {-1.1, 1.1}}, {"phi_23", {-1.1, 1.1}},
      {"phi_24", {-1.1, 1.1}}, {"phi_25", {-1.1, 1.1}}, {"phi_26", {-1.1, 1.1}},
      {"phi_27", {-1.1, 1.1}}, {"phi_28", {-1.1, 1.1}}, {"phi_29", {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail),
                 Parameter("integral_to_compute", map_threshold_integral),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto pst = PST(&spatial, p_pst);

  std::vector<SPA*> spatials{&spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial};
  OPE ope_ac_grains(spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");

  ope_ac_grains.overload_nl_solver(
      NLSolverType::LBFGS,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));

  std::vector<Coefficients> vcoeff{
      coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains,
      coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains,
      coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains,
      coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains,
      coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains,
      coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains, coef_ac_grains};
  PB problem_ac_grains(ope_ac_grains, var_grains, vcoeff, pst);
  // // Coupling 1
  auto cc = Coupling("Multigrains ", problem_ac_grains);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.2;  // 50.0;
  const auto& dt = 1.e-1;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  // time.get_tree();
  time.solve();
  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
