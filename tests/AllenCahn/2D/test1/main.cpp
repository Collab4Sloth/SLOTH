/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief  Steady solution of the Allen-Cahn equations
 * @version 0.1
 * @date 2025-07-11
 *
 * Copyright CEA (c) 2025
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
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

  using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
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

  // std::vector<std::string> vect_elem{"InlineSquareWithQuadrangles", "InlineSquareWithTriangles"};
  std::vector<std::string> vect_elem{"InlineSquareWithQuadrangles"};
  // std::vector<int> vect_order{1, 2};
  std::vector<int> vect_order{1};
  std::vector<int> vect_NN{30};  //, 60, 90, 120};
  for (const auto elem_type : vect_elem) {
    for (const auto order : vect_order) {
      for (const auto NN : vect_NN) {
        const int order_fe = order;      // finite element order
        const int refinement_level = 0;  // number of levels of uniform refinement
        const int nx = NN;
        const int ny = NN;
        const double lx = 1.;
        const double ly = 1.;
        const std::tuple<int, int, double, double>& tuple_of_dimensions = std::make_tuple(
            nx, ny, lx, ly);  // Number of elements and maximum length in each direction

        SPA spatial(elem_type, order_fe, refinement_level, tuple_of_dimensions);
        // ##############################
        //     Boundary conditions     //
        // ##############################
        auto boundaries_phi = {
            Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 1.),
            Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Dirichlet", 0.)};
        auto bcs_phi = BCS(&spatial, boundaries_phi);

        // ###########################################
        // ###########################################
        //            Physical models               //
        // ###########################################
        // ###########################################
        // ####################
        //     parameters    //
        // ####################
        //  Interface thickness
        const double epsilon(0.1);
        // Interfacial energy
        const double sigma(1.);
        // Two-phase mobility
        const double mob(1.);
        const double lambda = 1.5 * sigma * epsilon;
        const double omega = 12. * sigma / epsilon;
        auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                                 Parameter("lambda", lambda), Parameter("omega", omega));
        // ####################
        //     variables     //
        // ####################
        const double center_x = 0.;
        const double center_y = 0.;
        const double a_x = 1.;
        const double a_y = 0.;
        const double thickness = 5.e-3;
        const double radius = 5.e-1;

        auto initial_condition =
            AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y,
                                     a_x, a_y, thickness, radius);
        auto analytical_solution =
            AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y,
                                     a_x, a_y, epsilon, radius);

        const std::string& var_name_1 = "phi";
        auto v1 = VAR(&spatial, bcs_phi, var_name_1, 2, initial_condition, analytical_solution);
        auto vars = VARS(v1);

        // ###########################################
        // ###########################################
        //      Post-processing                     //
        // ###########################################
        // ###########################################

        const std::string& main_folder_path =
            "Saves_order_" + std::to_string(order) + "_Nx" + std::to_string(NN);
        const int level_of_detail = 1;
        const int frequency = 1;
        std::string calculation_path = "Problem1";
        bool enable_save_specialized_at_iter = true;
        auto p_pst = Parameters(
            Parameter("main_folder_path", main_folder_path),
            Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
            Parameter("level_of_detail", level_of_detail),
            Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
        // ####################
        //     operators     //
        // ####################

        // Problem 1:
        std::vector<SPA*> spatials{&spatial};
        OPE oper(spatials, params, TimeScheme::EulerImplicit);
        oper.overload_mobility(Parameters(Parameter("mob", mob)));
        oper.overload_nl_solver(
            NLSolverType::NEWTON,
            Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                       Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-12),
                       Parameter("iter_max", 1000)));
        const auto& solver = HypreSolverType::HYPRE_GMRES;
        const auto& precond = HyprePreconditionerType::HYPRE_ILU;
        oper.overload_solver(solver);
        oper.overload_preconditioner(precond);

        auto phi_cvg = PhysicalConvergence(ConvergenceType::ABSOLUTE_MAX, 1.e-12);
        auto CVG = Convergence(phi_cvg);
        auto pst = PST(&spatial, p_pst);

        PB problem1(oper, vars, pst, CVG);

        // Coupling 1
        auto cc = Coupling("AllenCahn Coupling", problem1);

        // ###########################################
        // ###########################################
        //            Time-integration              //
        // ###########################################
        // ###########################################
        const double t_initial = 0.0;
        const double t_final = 2.e-3;
        const double dt = 1.e-3;
        auto time_params = Parameters(Parameter("initial_time", t_initial),
                                      Parameter("final_time", t_final), Parameter("time_step", dt));
        auto time = TimeDiscretization(time_params, cc);

        time.solve();
        //---------------------------------------
        // Profiling stop
        //---------------------------------------
        Profiling::getInstance().print();
      }
    }
  }

  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
