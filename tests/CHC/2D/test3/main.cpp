/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Allen-Cahn problem solved in a square
 * @version 0.1
 * @date 2024-05-23
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>

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
  setVerbosity(Verbosity::Debug);

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
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  using NLFI = CahnHilliardNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                            ThermodynamicsPotentials::WW, Mobility::Constant>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
  using PB = Problem<OPE, VARS, PST>;
  using PB1 = MPI_Problem<VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string mesh_type =
      "InlineSquareWithHexahedres";  // type of mesh // "InlineSquareWithTriangles"
  const int order_fe = 1;            // finite element order
  const int refinement_level = 0;    // number of levels of uniform refinement
  SPA spatial("GMSH", 1, refinement_level, "slothLogo.msh", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Interface thickness
  // Interfacial energy
  const double sigma(1.);
  // Two-phase mobility
  const double mob(5.);
  const double lambda = 2.;
  const double omega = 5.;
  auto params =
      Parameters(Parameter("sigma", sigma), Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################

  auto user_func_solution =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double co = 0.5;
        double epsilon = 0.01;
        double xx = 200 * x[0];
        double yy = 200 * x[1];

        double sol =
            co + epsilon * (std::cos(0.105 * xx) * std::cos(0.11 * yy) +
                            (std::cos(0.13 * xx) * std::cos(0.087 * yy)) *
                                (std::cos(0.13 * xx) * std::cos(0.087 * yy)) +
                            (std::cos(0.025 * xx - 0.15 * yy) * std::cos(0.07 * xx - 0.02 * yy)));

        return sol;
      });

  auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto mu_initial_condition = 0.0;
  auto v1 = VAR(&spatial, bcs, "phi", 2, phi_initial_condition);
  auto v2 = VAR(&spatial, bcs, "mu", 2, mu_initial_condition);
  auto vars = VARS(v1, v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const int frequency = 10;
  std::string calculation_path = "CahnHilliard";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  std::vector<SPA*> spatials{&spatial, &spatial};
  OPE oper(spatials, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", "Newton solver "),
                                     Parameter("print_level", 1), Parameter("rel_tol", 1.e-8),
                                     Parameter("abs_tol", 1.e-11), Parameter("iter_max", 1000)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_BOOMER_AMG;
  oper.overload_solver(solver);
  oper.overload_preconditioner(precond);

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst);
  PB problem1(oper, vars, pst, convergence);

  // Coupling 1
  auto cc = Coupling("CahnHilliard Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 1.e4;
  const double dt = 1.;
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
