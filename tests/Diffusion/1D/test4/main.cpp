/**
 * @file main.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Comparaison analytical and numerical solution
 * @version 0.1
 * @date 2024-11-28
 *
 * Copyright CEA (c) 2024
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

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

  using NLFI = ThermoDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                               Diffusion::Constant>;
  using OPE = DiffusionOperator<FECollection, DIM, NLFI, Density::Constant>;
  using TD = Problem<OPE, VARS, PST>;

  using CC = Calphad_Problem<AnalyticalIdealSolution<mfem::Vector>, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  double L = 1e-3;
  int order = 1;
  int NN = 80;
  SPA spatial("InlineLineWithSegments", order, refinement_level, std::make_tuple(NN, L));
  // ##############################
  //     Boundary conditions     //
  // // ##############################
  auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  const auto& stabCoeff(1.e-8);
  const auto& diffusionCoeff(1.e-11);
  //  ####################
  //      variables     //
  //  ####################
  const auto& dt = 0.1;

  //  -----------------
  // THERMAL DIFFUSION
  //   -----------------
  auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
      [L, diffusionCoeff, dt](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        auto func = 0.5 * (1 + std::erf((xx - L / 2) / std::sqrt(4 * diffusionCoeff * 5. * dt)));

        return func;
      });

  auto user_func_analytical = std::function<double(const mfem::Vector&, double)>(
      [L, dt, diffusionCoeff](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        auto func =
            0.5 * (1 + std::erf((xx - L / 2) / std::sqrt(4 * diffusionCoeff * (time + 5. * dt))));

        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto analytical_solution = AnalyticalFunctions<DIM>(user_func_analytical);

  auto td_vars = VARS(VAR(&spatial, bcs, "c", 2, initial_condition));

  auto td_parameters = Parameters(Parameter("M", diffusionCoeff));
  //  -----------------
  // CALPHAD
  //  -----------------
  auto temp = VAR(&spatial, bcs, "T", 1, 1. / Physical::R);
  temp.set_additional_information("Temperature", "K");
  auto pres = VAR(&spatial, bcs, "pressure", 2, 1.);
  pres.set_additional_information("Pressure", "Pa");
  auto xo = VAR(&spatial, bcs, "O", 2, initial_condition, analytical_solution);
  xo.set_additional_information("O", "X");

  auto heat_vars = VARS(temp);
  auto p_vars = VARS(pres);
  auto compo_vars = VARS(xo);

  auto muo = VAR(&spatial, bcs, "muO", 2, 1.);
  muo.set_additional_information("O", "mu");
  auto mu_var = VARS(muo);

  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");
  auto calphad_parameters = Parameters(description_calphad);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 100;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst = PST(&spatial, p_pst);
  calculation_path = "Calphad";
  auto cpst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto Calphad_pst = PST(&spatial, cpst);
  // ####################
  //     PROBLEMS       //
  // ####################
  // THERMAL DIFFUSION
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, td_parameters, TimeScheme::EulerImplicit);
  oper.overload_diffusion(Parameters(Parameter("D", stabCoeff)));

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  // CALPHAD
  CC CALPHAD_PROBLEM(calphad_parameters, mu_var, Calphad_pst, convergence, heat_vars, p_vars,
                     compo_vars);
  TD THERMAL_DIFFUSION(oper, compo_vars, pst, convergence, mu_var);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", THERMAL_DIFFUSION, CALPHAD_PROBLEM);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 1.e3;
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
