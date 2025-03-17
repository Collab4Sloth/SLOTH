/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Coupling Calphad(BinaryTabulation)/HeatTransfer : covering test
 * @version 0.1
 * @date 2025-01-14
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Quiet);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  constexpr int DIM = Test<1>::dim;
  using FECollection = Test<1>::FECollection;
  using VARS = Test<1>::VARS;
  using VAR = Test<1>::VAR;
  using PSTCollection = Test<1>::PSTCollection;
  using PST = Test<1>::PST;
  using SPA = Test<1>::SPA;
  using BCS = Test<1>::BCS;
  /////////////////////////
  using PB_TABULATED = Calphad_Problem<BinaryTabulation<mfem::Vector>, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  SPA spatial("InlineLineWithSegments", 1, refinement_level, std::make_tuple(1000, 1.));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto Calphadboundaries = boundaries;
  auto Tboundaries = boundaries;
  auto Calphadbcs = BCS(&spatial, Calphadboundaries);

  auto Tbcs = BCS(&spatial, Tboundaries);
  // ####################
  //     parameters    //
  // ####################
  // Heat
  const auto& rho(1.);
  const auto& cp(1.);
  const auto& cond(2.);

  // ############################
  //     variables IC + SRC    //
  // ###########################
  const auto& pellet_radius = 0.00465;
  // Heat

  auto temp = VAR(&spatial, Tbcs, "T", 2, 750.);
  temp.set_additional_information("Temperature", "K");

  auto heat_vars = VARS(temp);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  auto user_func_solution =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& v, double time) {
        const double x = v[0];
        double max_val = 0.1;
        const auto func = std::max(max_val, std::min(x, 0.9));
        return func;
      });

  auto muo = VAR(&spatial, Calphadbcs, "muO", 2, 0.);
  muo.set_additional_information("O", "mu");
  auto xso = VAR(&spatial, Calphadbcs, "xsO", 2, 0.);
  xso.set_additional_information("O", "SOLUTION", "x");
  auto gs = VAR(&spatial, Calphadbcs, "gs", 2, 0.);
  gs.set_additional_information("SOLUTION", "g");
  auto outputs = VARS(muo, xso, gs);

  auto pres = VAR(&spatial, Calphadbcs, "pressure", 2, 50.e5);
  pres.set_additional_information("Pressure", "Pa");
  auto p_vars = VARS(pres);

  auto x_o_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto xo = VAR(&spatial, Calphadbcs, "O", 2, x_o_initial_condition);
  xo.set_additional_information("O", "X");
  auto compo_vars = VARS(xo);

  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");
  auto calphad_parameters =
      Parameters(description_calphad, Parameter("data filename", "data.thermotab"));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Calphad_Tabulated";
  auto cpst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto Calphad_pst = PST(&spatial, cpst);

  // ####################
  //     Problems      //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  //---------------
  // Calphad
  //---------------
  PB_TABULATED Calphad_tabulated(calphad_parameters, outputs, Calphad_pst, convergence, heat_vars,
                                 p_vars, compo_vars);

  // ####################
  //     Coupling      //
  // ####################
  auto cc = Coupling("Calphad calculation", Calphad_tabulated);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.25;
  const auto& dt = 0.25;
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
