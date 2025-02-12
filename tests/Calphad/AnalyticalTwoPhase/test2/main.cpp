/**
 * @file main.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Coupling Calphad(AnalyticalIdealSolution)/HeatTransfer/Allen-Cahn : covering test
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
  const auto DIM = 1;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VARS = Variables<FECollection, DIM>;
  using VAR = Variable<FECollection, DIM>;
  using PB_CAL = Calphad_Problem<AnalyticalParaboloidForTwoPhase_2<mfem::Vector>, VARS, PST>;
  using PB_KKS = Interface_problem<AnalyticalParaboloidForTwoPhase_2<mfem::Vector>, VARS, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  double L = 2;  // 4.65e-3;
  int NN = 3000;
  SpatialDiscretization<FECollection, DIM> spatial("InlineLineWithSegments", 1, refinement_level,
                                                   std::make_tuple(NN, L));

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto Calphadboundaries = {Boundary("left", 0, "Neumann", 0.),
                            Boundary("right", 1, "Neumann", 0.)};
  auto Calphadbcs = BoundaryConditions<FECollection, DIM>(&spatial, Calphadboundaries);
  auto Tboundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};

  auto Tbcs = BoundaryConditions<FECollection, DIM>(&spatial, Tboundaries);
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

  const auto& center_x = 0.;
  // const double y = v[1];
  const double y = 0.;
  const auto& a_x = 1.;
  const auto& a_y = 1.;
  const auto& thickness = 5.e-4;
  const auto& radius = 0.45 * pellet_radius;
  double epsilon = 3 * L / NN;

  auto muO = std::function<double(const mfem::Vector&, double)>(
      [L, epsilon](const mfem::Vector& v, double time) {
        const double x = v[0];
        // const double y = v[1];
        const double y = 0.;
        const auto r = std::sqrt(x * x + y * y);
        double func;
        if (r < (0.5 * L)) {
          func = 0.4;  // 0.4
        } else {
          func = 0.3;  // 0.3
        }
        // func = (0.5 * (0.3 + 0.4) + 0.5 * (0.3 - 0.4) * std::tanh((r - 0.5 * L) / epsilon));
        return func;
      });

  auto muu_ = std::function<double(const mfem::Vector&, double)>(
      [L, epsilon](const mfem::Vector& v, double time) {
        const double x = v[0];
        // const double y = v[1];
        const double y = 0.;
        const auto r = std::sqrt(x * x + y * y);
        double func;
        if (r < (0.5 * L)) {
          func = -0.125 + 0.3;  // 0.175
        } else {
          func = 0.2 + 0.4;  // 0.6
        }
        // func = (0.5 * (0.6 + 0.175) + 0.5 * (0.6 - 0.175) * std::tanh((r - 0.5 * L) / epsilon));
        return func;
      });

        auto acIC = std::function<double(const mfem::Vector&, double)>(
      [L, epsilon](const mfem::Vector& v, double time) {
        const double x = v[0];
        const double y = 0.;
        const auto r = std::sqrt(x * x + y * y);
        double func;
        func = 1 - (0.5 + 0.5 * std::tanh((r - 0.5 * L) / epsilon));
        return func;
      });

  auto ac_ic = AnalyticalFunctions<DIM>(acIC);
  auto muO_ = AnalyticalFunctions<DIM>(muO);
  auto muU_ = AnalyticalFunctions<DIM>(muu_);
  auto mass_ic1 = AnalyticalFunctions<DIM>(muO_);
  auto mass_ic2 = AnalyticalFunctions<DIM>(muU_);

  auto vars_ac = VAR(&spatial, Calphadbcs, "phi", 2, ac_ic);
  vars_ac.set_additional_information("Phase indicator", "-");

  auto ac_vars = VARS(vars_ac);
  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################

  // mu
  auto muo = VAR(&spatial, Calphadbcs, "muO", 2, 0.);
  muo.set_additional_information("O", "mu");
  auto muu = VAR(&spatial, Calphadbcs, "muU", 2, 0.);
  muu.set_additional_information("U", "mu");

  // x
  auto xso = VAR(&spatial, Calphadbcs, "xsO", 2, 0.);
  xso.set_additional_information("O", "SOLID", "x");
  auto xsu = VAR(&spatial, Calphadbcs, "xsU", 2, 0.);
  xsu.set_additional_information("U", "SOLID", "x");

  auto xlo = VAR(&spatial, Calphadbcs, "xlO", 2, 0.);
  xlo.set_additional_information("O", "LIQUID", "x");
  auto xlu = VAR(&spatial, Calphadbcs, "xlU", 2, 0.);
  xlu.set_additional_information("U", "LIQUID", "x");

  // g
  // SOLID
  auto gs = VAR(&spatial, Calphadbcs, "gs", 2, 0.);
  gs.set_additional_information("SOLID", "g");
  // LIQUID
  auto gl = VAR(&spatial, Calphadbcs, "gl", 2, 0.);
  gl.set_additional_information("LIQUID", "g");

  auto pres = VAR(&spatial, Calphadbcs, "pressure", 2, 50.e5);
  pres.set_additional_information("Pressure", "Pa");
  auto p_vars = VARS(pres);

  auto outputs = VARS(gs, gl, muo, xso, xlo, muu, xsu, xlu);

  auto xo = VAR(Variable<FECollection, DIM>(&spatial, Calphadbcs, "O", 2, mass_ic1));
  xo.set_additional_information("O", "X");
  auto xu = VAR(Variable<FECollection, DIM>(&spatial, Calphadbcs, "U", 2, mass_ic2));
  xu.set_additional_information("U", "X");
  auto compo_vars1 = VARS(xo);
  auto compo_vars2 = VARS(xu);
  // auto compo_vars = VARS(xo,xu);
  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");

  std::map<std::tuple<std::string, std::string>, mfem::real_t> coeff_k;
  coeff_k[std::make_tuple("SOLID", "O")] = 1.0;
  coeff_k[std::make_tuple("SOLID", "U")] = 1.0;
  coeff_k[std::make_tuple("LIQUID", "O")] = 1.0;
  coeff_k[std::make_tuple("LIQUID", "U")] = 1.0;

  auto coeff_ks = Parameter("coefficient_k", coeff_k);

  std::map<std::tuple<std::string, std::string>, mfem::real_t> c_eq;
  c_eq[std::make_tuple("SOLID", "O")] = 0.3;
  c_eq[std::make_tuple("SOLID", "U")] = 0.3;
  c_eq[std::make_tuple("LIQUID", "O")] = 0.4;
  c_eq[std::make_tuple("LIQUID", "U")] = 0.4;
  auto c_eqs = Parameter("equilibrium_composition", c_eq);
  auto calphad_parameters = Parameters(description_calphad, coeff_ks, c_eqs);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Calphad";
  auto cpst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto Calphad_pst = PST(&spatial, cpst);
  calculation_path = "KKS";
  auto kkspst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto KKS_pst = PST(&spatial, kkspst);

  // ####################
  //     Problems      //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  //---------------
  // Calphad
  //---------------
  PB_CAL Calphad_pb(calphad_parameters, outputs, Calphad_pst, convergence, heat_vars, p_vars, ac_vars,
                compo_vars1, compo_vars2);
  //---------------
  // KKS
  //---------------               
  PB_KKS KKS_pb(calphad_parameters, outputs, KKS_pst, convergence, heat_vars, p_vars, ac_vars,
                compo_vars1, compo_vars2);
  // ####################
  //     Coupling      //
  // ####################
  auto cc = Coupling("Calphad calculation", Calphad_pb,KKS_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 1.;
  const auto& dt = 0.1;
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
