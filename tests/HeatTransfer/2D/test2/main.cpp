/**
 * @file main.cpp
 * @author mh286406 (marine.harel@cea.fr)
 * @brief 2D heat transfer problem
 * @version 0.1
 * @date 2026-03-02
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

#include "./RobinCoefficient.hpp"
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
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;

  using OPE = SteadyOperator<FECollection, DIM>;
  using PB = Problem<OPE, VARS, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string& mesh_type = "InlineSquareWithQuadrangles";  // type of mesh
  const int order_fe = 1;                                        // finite element order
  const int refinement_level = 0;  // number of levels of uniform refinement

  SPA spatial(mesh_type, order_fe, refinement_level, std::make_tuple(100, 50, 0.02, 0.01));

  // ##############################
  //     Boundary conditions     //
  // ##############################
  const double T_d = 1173.0;
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", T_d),
                     Boundary("upper", 2, "Robin"), Boundary("left", 3, "Dirichlet", T_d)};
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

  const double k(3.0);
  const double h(50.0);
  const double T_inf(323.0);
  const double epsilon(0.7);
  const double sigma(5.669e-8);
  Coefficient conductivity(Glossary::Conductivity, k);
  Coefficient robin_a(Glossary::Robin_a, Scheme::Implicit, RobinCoefficient());
  Coefficient robin_b(Glossary::Robin_b, h * T_inf + epsilon * sigma * std::pow(T_inf, 4));
  robin_a.set_bdr_index_coef(std::vector<int>{2});
  robin_b.set_bdr_index_coef(std::vector<int>{2});

  // ####################
  //     variables     //
  // ####################

  auto heat_vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, T_d));

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
  Coefficients coef_pb(conductivity, robin_a, robin_b);
  std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"Fourier"});

  oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
  PB Heat_pb("Heat", oper, heat_vars, {coef_pb}, pst);

  // Coupling 1
  auto cc = Coupling("Heat transfer", Heat_pb);
  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.1;
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
