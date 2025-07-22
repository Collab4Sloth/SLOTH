/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Spinodal decomposition in a 3D domain
 * @version 0.1
 * @date 2024-09-3
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Quiet);
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 3;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////
  using NLFI = AllenCahnGrainNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                              ThermodynamicsPotentials::F, Mobility::Constant>;
  using LHS_NLFI = TimeNLFormIntegrator<VARS>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, LHS_NLFI>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################

  const int order_fe = 1;
  auto refinement_level = 0;
  SPA spatial("GMSH", order_fe, refinement_level, "mesh_poly.msh", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("tata1", 0, "Neumann", 0.),   Boundary("tata2", 1, "Neumann", 0.),
                     Boundary("tata3", 2, "Neumann", 0.),   Boundary("tata4", 3, "Neumann", 0.),
                     Boundary("tata5", 4, "Neumann", 0.),   Boundary("tata6", 5, "Neumann", 0.),
                     Boundary("tata7", 6, "Neumann", 0.),   Boundary("tata8", 7, "Neumann", 0.),
                     Boundary("tata9", 8, "Neumann", 0.),   Boundary("tata10", 9, "Neumann", 0.),
                     Boundary("tata11", 10, "Neumann", 0.), Boundary("tata12", 11, "Neumann", 0.),
                     Boundary("tata13", 12, "Neumann", 0.), Boundary("tata14", 13, "Neumann", 0.),
                     Boundary("tata15", 14, "Neumann", 0.), Boundary("tata16", 15, "Neumann", 0.),
                     Boundary("tata17", 16, "Neumann", 0.), Boundary("tata18", 17, "Neumann", 0.),
                     Boundary("tata19", 18, "Neumann", 0.), Boundary("tata20", 19, "Neumann", 0.),
                     Boundary("tata21", 20, "Neumann", 0.), Boundary("tata22", 21, "Neumann", 0.),
                     Boundary("tata23", 22, "Neumann", 0.), Boundary("tata24", 23, "Neumann", 0.),
                     Boundary("tata25", 24, "Neumann", 0.), Boundary("tata26", 25, "Neumann", 0.),
                     Boundary("tata27", 26, "Neumann", 0.), Boundary("tata28", 27, "Neumann", 0.),
                     Boundary("tata29", 28, "Neumann", 0.), Boundary("tata30", 29, "Neumann", 0.),
                     Boundary("tata31", 30, "Neumann", 0.), Boundary("tata32", 31, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  const auto& mob(5.);
  const auto& lambda = 0.1;  // 3. * sigma * epsilon / 2.;
  const auto& omega = 1.;    // 12. * sigma / epsilon; sig = eps/12 ;
  const auto& epsilon(std::sqrt(8. * lambda));
  const auto& sigma(epsilon / 12.);

  auto params =
      Parameters(Parameter("epsilon", epsilon), Parameter("mobility", mob),
                 Parameter("sigma", sigma), Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################

  auto v1 = VAR(&spatial, bcs, "phi1", 2, 1, {"toto1"});
  // auto v2 = VAR(&spatial, bcs, "phi2", 2, 1, {"toto2"});
  auto v3 = VAR(&spatial, bcs, "phi3", 2, 1, {"toto3"});
  auto v4 = VAR(&spatial, bcs, "phi4", 2, 1, {"toto4"});
  auto v5 = VAR(&spatial, bcs, "phi5", 2, 1, {"toto5"});
  auto v6 = VAR(&spatial, bcs, "phi6", 2, 1, {"toto6"});
  auto v7 = VAR(&spatial, bcs, "phi7", 2, 1, {"toto7"});
  auto v8 = VAR(&spatial, bcs, "phi8", 2, 1, {"toto8"});
  auto v9 = VAR(&spatial, bcs, "phi9", 2, 1, {"toto9"});
  auto v10 = VAR(&spatial, bcs, "phi10", 2, 1, {"toto10"});
  auto v11 = VAR(&spatial, bcs, "phi11", 2, 1, {"toto11"});
  auto v12 = VAR(&spatial, bcs, "phi12", 2, 1, {"toto12"});
  auto v13 = VAR(&spatial, bcs, "phi13", 2, 1, {"toto13"});
  auto v14 = VAR(&spatial, bcs, "phi14", 2, 1, {"toto14"});
  auto v15 = VAR(&spatial, bcs, "phi15", 2, 1, {"toto15"});
  auto v16 = VAR(&spatial, bcs, "phi16", 2, 1, {"toto16"});
  auto v17 = VAR(&spatial, bcs, "phi17", 2, 1, {"toto17"});
  auto v18 = VAR(&spatial, bcs, "phi18", 2, 1, {"toto18"});
  auto v19 = VAR(&spatial, bcs, "phi19", 2, 1, {"toto19"});
  auto v20 = VAR(&spatial, bcs, "phi20", 2, 1, {"toto20"});
  auto v21 = VAR(&spatial, bcs, "phi21", 2, 1, {"toto21"});
  auto v22 = VAR(&spatial, bcs, "phi22", 2, 1, {"toto22"});
  auto v23 = VAR(&spatial, bcs, "phi23", 2, 1, {"toto23"});
  auto v24 = VAR(&spatial, bcs, "phi24", 2, 1, {"toto24"});
  auto v25 = VAR(&spatial, bcs, "phi25", 2, 1, {"toto25"});
  auto v26 = VAR(&spatial, bcs, "phi26", 2, 1, {"toto26"});
  auto v27 = VAR(&spatial, bcs, "phi27", 2, 1, {"toto27"});
  auto v28 = VAR(&spatial, bcs, "phi28", 2, 1, {"toto28"});
  auto v29 = VAR(&spatial, bcs, "phi29", 2, 1, {"toto29"});
  auto v30 = VAR(&spatial, bcs, "phi30", 2, 1, {"toto30"});
  auto v31 = VAR(&spatial, bcs, "phi31", 2, 1, {"toto31"});

  auto vars = VARS(v1, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19,
                   v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 50;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail),
                 Parameter("enable_save_specialized_at_iter", true),
                 Parameter("force_clean_output_dir", true));
  auto pst = PST(&spatial, p_pst);

  // ####################
  //     operator     //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  std::vector<SPA*> spatials{&spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial,
                             &spatial, &spatial, &spatial, &spatial, &spatial, &spatial};
  OPE oper(spatials, params, TimeScheme::EulerImplicit);

  oper.overload_mobility(Parameters(Parameter("mob", mob)));

  setVerbosity(Verbosity::Debug);

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  // ####################
  //     Problem       //
  // ####################
  PB problem1(oper, vars, pst, convergence);

  // ####################
  //     Coupling      //
  // ####################
  auto cc = Coupling("Default Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 500.;
  const auto& dt = 1.e-1;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
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
  //---------------------------------------
  return 0;
}
