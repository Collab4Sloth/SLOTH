/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Multigrains simulations from "Programming Phase-field modeling - Biner 2017" with AMR and
 * partioned scheme
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
#include <utility>

#include "./GrainsCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "voro++.hh"

using namespace voro;

double distance(double x1, double y1, double x2, double y2, double L) {
  double dx = std::abs(x2 - x1);
  double dy = std::abs(y2 - y1);
  // periodicity
  dx = std::min(dx, L - dx);
  dy = std::min(dy, L - dy);

  return std::sqrt(dx * dx + dy * dy);
}

/////////////////////////

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  using namespace Sloth2D;
  setVerbosity(Verbosity::Normal);
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

  constexpr int ngrains = 30;

  // ###########################################
  // ###########################################
  //            Polycristaline                //
  // ###########################################
  // ###########################################
  const double L = 32.;
  const double SUBDIV = 8;

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
                SUBDIV, SUBDIV, 1,  // Number of grid subdivisions
                true, true, false,  // No periodic boundaries
                SUBDIV);

  for (int i = 0; i < data.size(); ++i) {
    con.put(i, std::get<0>(data[i]), std::get<1>(data[i]), 0.0);
  }

  std::vector<std::function<double(const mfem::Vector&, double)>> vect_user_func;
  for (int id_grain = 0; id_grain < static_cast<int>(data.size()); ++id_grain) {
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
              // Calculate the distance from the grid point to the current particle
              double dist = distance(x[0], x[1], xx, y, L);
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
  int refinement_level = 0;
  int NN = 8;
  // Create translation vectors defining the periodicity
  mfem::Vector x_translation({L, 0.0});
  mfem::Vector y_translation({0.0, L});
  std::vector<mfem::Vector> translations = {x_translation, y_translation};

  SPAS spatials =
      setPeriodicSpatialDiscretization(30, "InlineSquareWithQuadrangles", 1, refinement_level,
                                       std::make_tuple(NN, NN, L, L), translations, true);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                     Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
  auto bcs = setBoundaryConditions(ngrains, spatials, boundaries);

  // ###########################################
  // ###########################################
  //      Variables
  // ###########################################
  // ###########################################
  std::vector<VARS> vars_grains;
  vars_grains.reserve(ngrains);
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
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 2;
  const bool enable_save_specialized_at_iter = true;
  std::string calculation_path = "GrainsProblem";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto pst_vec = PST(spatials[0], p_pst);

  // ###########################################
  // ###########################################
  //     Coefficients
  // ###########################################
  // ###########################################
  const double& mob(5.);
  const double& lambda(0.1);
  Coefficient grad_energy_i(Glossary::GradEnergy, Scheme::Implicit, GrainGradient());
  Coefficient double_well_exp_i(Glossary::FreeEnergy, Scheme::Implicit, ExplicitGrainGw());
  Coefficient capillary_i(Glossary::Capillary, lambda);
  Coefficient mobility_i(Glossary::Mobility, mob);
  Coefficient swithchingFct(Glossary::PhaseField, Scheme::Implicit, SwitchingFunction());
  swithchingFct.set_name("SWF");
  std::vector<Coefficients> coeffs_i{
      Coefficients(double_well_exp_i, capillary_i, mobility_i, grad_energy_i)};

  // ###########################################
  // ###########################################
  //     Problems
  // ###########################################
  // ###########################################
  int level_amr = 3;  // 4;
  auto amr_params =
      Parameters(Parameter("max_elem_error", 1.e-3), Parameter("amr_max_level", level_amr),
                 Parameter("nc_limit", 0), Parameter("max_preref_cycles", level_amr));
  auto amr_coef = mfem::ConstantCoefficient(1.0);

  std::vector<TransientPB> ac_pbs;
  ac_pbs.reserve(ngrains);

  for (int i = 0; i < ngrains; ++i) {
    // Operator
    TransientOPE ope_i({spatials[i]}, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");

    // Auxiliary variables
    std::vector<VARS*> aux_ptrs;
    aux_ptrs.reserve(ngrains - 1);
    for (int j = 0; j < ngrains; ++j) {
      if (j != i) aux_ptrs.push_back(&vars_grains[j]);
    }
    ac_pbs.emplace_back(ope_i, vars_grains[i], coeffs_i, pst_vec);
    ac_pbs.back().set_auxvariables(std::move(aux_ptrs));

    // AMR
    auto* amr_integ_i = new mfem::DiffusionIntegrator(amr_coef);
    auto* estimator_i = new SlothErrorEstimators(ErrorEstimatorType::KELLY, amr_integ_i);
    auto* amr_i =
        new MultiVariableMaxAMR<VARS>(*spatials[i]->get_mesh(), spatials[i]->is_nc_simplices());
    amr_i->SetCriteria(estimator_i, amr_params);
    ac_pbs.back().set_amr(amr_i);
    //
  }
  ac_pbs.back().set_vtk_coefficients({swithchingFct});

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
  const double t_initial = 0.0;
  const auto& t_final = 0.2;  // 50.0;
  const auto& dt = 1.e-1;
  auto time_params =
      Parameters(Parameter("initial_time", t_initial), Parameter("vtk_unified", true),
                 Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  time.solve();

  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();

  deleteSpatialDiscretization(spatials);
  //---------------------------------------
  return 0;
}
