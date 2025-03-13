/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D heat transfer  problem
 * @version 0.1
 * @date 2024-09-3
 *
 * @copyright Copyright (c) 2024
 *
 */
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

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  const int DIM=2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  // Heat
  using NLFI2 =
      HeatNLFormIntegrator<VARS, CoefficientDiscretization::Explicit, Conductivity::Constant>;
  using OPE2 = HeatOperator<FECollection, DIM, NLFI2, Density::Constant, HeatCapacity::Constant>;
  using PB2 = Problem<OPE2, VARS, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  SPA spatial("GMSH", 1, refinement_level, "camembert2D.msh", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto Tboundaries = {Boundary("lower", 0, "Neumann", 0.),
                      Boundary("external", 2, "Dirichlet", 750.),
                      Boundary("upper", 1, "Neumann", 0.)};
  auto Tbcs = BCS(&spatial, Tboundaries);
  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
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
  auto heat_vars = VARS(VAR(&spatial, Tbcs, "T", 2, 750.));
  auto pl = 4.e4;
  auto src_func = std::function<double(const mfem::Vector&, double)>(
      [pl, pellet_radius](const mfem::Vector& vcoord, double time) {
        const auto func = pl / (M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });

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
  auto p_pst2 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst2 = PST(&spatial, p_pst2);

  // ####################
  //     Problems      //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  // Heat:
  auto source_term = AnalyticalFunctions<DIM>(src_func);
  OPE2 oper2(&spatial, TimeScheme::EulerImplicit, source_term);
  oper2.overload_density(Parameters(Parameter("rho", rho)));
  oper2.overload_heat_capacity(Parameters(Parameter("cp", cp)));
  oper2.overload_conductivity(Parameters(Parameter("lambda", cond)));

  oper2.overload_nl_solver(NLSolverType::NEWTON,
                           Parameters(Parameter("description", "Newton solver "),
                                      Parameter("print_level", 1), Parameter("abs_tol", 1.e-6)));
  PB2 Heat_pb("Heat", oper2, heat_vars, pst2, convergence);

  // MPI_Problem<VAR, PST> mpi(vars3, pst3, convergence);

  // Coupling 1
  auto cc = Coupling("Heat transfer", Heat_pb);
  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 1.;
  const auto& dt = 0.25;
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
