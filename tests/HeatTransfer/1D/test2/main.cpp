/**
 * @file main.cpp
 * @author mh286406 (marine.harel@cea.fr)
 * @brief 1D heat transfer problem
 * @version 0.1
 * @date 2026-02-25
 *
 * @copyright Copyright (c) 2026
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Verbose);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM = 1;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;

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
  const std::string& mesh_type = "InlineLineWithSegments";  // type of mesh
  const int order_fe = 1;                                   // finite element order
  const int refinement_level = 0;  // number of levels of uniform refinement
  const int n_elements = 30;
  const double length = 1e-3;

  SPA spatial(mesh_type, order_fe, refinement_level, std::make_tuple(n_elements, length));

  // ##############################
  //     Boundary conditions     //
  // ##############################
  const double T_initial = 0.0;
  auto boundaries = {Boundary("left", 0, "Robin"), Boundary("right", 1, "Dirichlet", T_initial)};
  auto bcs = BCS(&spatial, boundaries);
  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  // Heat
  const auto& rho(10000.0);
  const auto& cp(10000.0);
  const auto& cond(1.0);
  const auto& diffusivity(cond / (rho * cp));
  const auto& convection(10.0);
  const auto& T_infinity(1.0);
  Coefficient density(Glossary::Concentration, rho);
  Coefficient heat_capacity(Glossary::Cp, cp);
  Coefficient conductivity(Glossary::Conductivity, cond);
  Coefficient robin_a(Glossary::Robin_a, convection);
  Coefficient robin_b(Glossary::Robin_b, convection * T_infinity);
  robin_a.set_bdr_index_coef(std::vector<int>{0});
  robin_b.set_bdr_index_coef(std::vector<int>{0});

  // ####################
  //     variables     //
  // ####################

  auto user_func_analytical = std::function<double(const mfem::Vector&, double)>(
      [diffusivity, T_infinity, T_initial, cond, convection](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        const auto L_c = std::sqrt(4 * diffusivity * time);
        const auto hk = convection / cond;
        const auto func =
            T_initial + (T_infinity - T_initial) *
                            (std::erfc(xx / L_c) - std::exp(hk * xx + 0.25 * hk * hk * L_c * L_c) *
                                                       std::erfc(xx / L_c + 0.5 * hk * L_c));
        return func;
      });
  auto T_analytical = AnalyticalFunctions<DIM>(user_func_analytical);
  auto heat_vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, T_initial, T_analytical));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  // Heat
  const std::string& calculation_path = "Heat";
  auto p_pst = Parameters(
      Parameter("main_folder_path", main_folder_path),
      Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
      Parameter("level_of_detail", level_of_detail), Parameter("enable_compute_energies", false));
  auto pst = PST(&spatial, p_pst);

  // ####################
  //     Problems      //
  // ####################

  // Heat
  Coefficients coef_pb(density, heat_capacity, conductivity, robin_a, robin_b);
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"Fourier"}, TimeScheme::EulerImplicit, "HeatTimeDerivative");

  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", "Newton solver "),
                                     Parameter("print_level", 1), Parameter("abs_tol", 1.e-10)));
  PB Heat_pb("Heat", oper, heat_vars, {coef_pb}, pst);

  // Coupling 1
  auto cc = Coupling("Heat transfer", Heat_pb);
  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 5.0;
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
