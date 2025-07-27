/**
 * Copyright CEA (c) 2025
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief MMS convergence tests :
 * Liangzhe Zhang, Michael R. Tonks, Derek Gaston, John W. Peterson, David Andrs, Paul C. Millett,
 * Bulent S. Biner, A quantitative comparison between C0 and C1 elements for solving the
 * Cahn–Hilliard equation, Journal of Computational Physics, Volume 236, 2013, Pages 74-80, ISSN
 * 0021-9991, https://doi.org/10.1016/j.jcp.2012.12.001.
 *
 * @version 0.1
 * @date 2025-07-04
 *
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
  setVerbosity(Verbosity::Verbose);

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
                                            ThermodynamicsPotentials::F, Mobility::Constant>;

  using LHS_NLFI = TimeCHNLFormIntegrator<VARS>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, LHS_NLFI>;
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
  const int refinement_level = 0;  // number of levels of uniform refinement

  // std::vector<std::string> vect_elem{"InlineSquareWithQuadrangles", "InlineSquareWithTriangles"};
  std::vector<std::string> vect_elem{"InlineSquareWithQuadrangles"};
  std::vector<int> vect_order{2, 1};
  std::vector<int> vect_NN{160, 80, 40, 20};
  for (const auto elem_type : vect_elem) {
    for (const auto order_fe : vect_order) {
      for (const auto NN : vect_NN) {
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
        auto boundaries = {
            Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 0.),
            Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Dirichlet", 0.)};
        auto bcs_phi = BCS(&spatial, boundaries);
        auto boundaries_mu = {
            Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Dirichlet", 0.),
            Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Dirichlet", 0.)};
        auto bcs_mu = BCS(&spatial, boundaries_mu);

        // ###########################################
        // ###########################################
        //            Physical models               //
        // ###########################################
        // ###########################################
        // ####################
        //     parameters    //
        // ####################
        //  Interface thickness
        const double epsilon(1.0);
        // Interfacial energy
        const double sigma(1.);
        // Two-phase mobility
        const double mob(1.);
        const double lambda = 1.;
        const double omega = 1.;
        auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                                 Parameter("lambda", lambda), Parameter("omega", omega));
        // ####################
        //     variables     //
        // ####################

        // phi analytical solution (mms)
        auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
            [](const mfem::Vector& x, double time) {
              const double xx = x[0];
              const auto func = (time + 1) * std::sin(M_PI * xx);
              return func;
            });

        // mu analytical solution (mms)
        auto mu_user_func_solution = std::function<double(const mfem::Vector&, double)>(
            [&](const mfem::Vector& x, double time) {
              const double xx = x[0];
              const auto c = (time + 1) * std::sin(M_PI * xx);
              const double delta_c = -M_PI * M_PI * c;

              const auto sol = c * c * c - c - delta_c;
              return sol;
            });

        // 1st equation source term
        auto user_func_source_term = std::function<double(const mfem::Vector&, double)>(
            [&](const mfem::Vector& x, double time) {
              const double xx = x[0];
              const double yy = x[1];
              const auto c = (time + 1) * std::sin(M_PI * xx);
              const double delta_c = -M_PI * M_PI * c;
              const double dc_dx = (time + 1) * M_PI * std::cos(M_PI * xx);
              const double delta_mu = delta_c * (3. * c * c - 1.) + 6. * c * (dc_dx * dc_dx) -
                                      M_PI * M_PI * M_PI * M_PI * c;
              const double dc_dt = std::sin(M_PI * xx);
              const auto func = dc_dt - delta_mu;
              return func;
            });
        auto user_func_source_term2 = std::function<double(const mfem::Vector&, double)>(
            [](const mfem::Vector& v, double time) { return 0.; });

        auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
        auto phi_analytical_condition = AnalyticalFunctions<DIM>(user_func_solution);
        auto mu_initial_condition = AnalyticalFunctions<DIM>(mu_user_func_solution);
        auto mu_analytical_condition = AnalyticalFunctions<DIM>(mu_user_func_solution);
        const std::string& var_name_1 = "phi";
        const std::string& var_name_2 = "mu";
        auto v1 =
            VAR(&spatial, bcs_phi, var_name_1, 2, phi_initial_condition, phi_analytical_condition);
        auto v2 = VAR(&spatial, bcs_mu, var_name_2, 2, mu_initial_condition);
        auto vars = VARS(v1, v2);

        // ###########################################
        // ###########################################
        //      Post-processing                     //
        // ###########################################
        // ###########################################
        const std::string& main_folder_path =
            "Saves_order_" + std::to_string(order_fe) + "_Nx" + std::to_string(NN);

        const int level_of_detail = 1;
        const int frequency = 1;
        std::string calculation_path = "Problem1";
        const double threshold = 10.;
        std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
            {var_name_1, {-10.1, 10.1}}};
        bool enable_save_specialized_at_iter = true;
        auto p_pst = Parameters(
            Parameter("main_folder_path", main_folder_path),
            Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
            Parameter("level_of_detail", level_of_detail),
            Parameter("integral_to_compute", map_threshold_integral),
            Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
        // ####################
        //     operators     //
        // ####################

        std::vector<SPA*> spatials{&spatial, &spatial};
        std::vector<AnalyticalFunctions<DIM>> src_term;
        src_term.emplace_back(AnalyticalFunctions<DIM>(user_func_source_term2));
        src_term.emplace_back(AnalyticalFunctions<DIM>(user_func_source_term));
        // OPE oper(spatials, params, src_term);
        OPE oper(spatials, params, TimeScheme::EulerImplicit, src_term);
        oper.overload_mobility(Parameters(Parameter("mob", mob)));
        oper.overload_nl_solver(
            NLSolverType::NEWTON,
            Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                       Parameter("rel_tol", 1.e-8), Parameter("abs_tol", 1.e-12)));
        const auto& solver = HypreSolverType::HYPRE_GMRES;
        const auto& precond = HyprePreconditionerType::HYPRE_ILU;
        oper.overload_solver(solver, Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
        oper.overload_preconditioner(precond, Parameters(Parameter("type", 1)));

        auto pst = PST(&spatial, p_pst);
        PB problem1(oper, vars, pst);

        // Coupling 1
        auto cc = Coupling("CahnHilliard Coupling", problem1);

        // ###########################################
        // ###########################################
        //            Time-integration              //
        // ###########################################
        // ###########################################
        const double t_initial = 0.0;
        const double t_final = 1.;
        const double dt = 1.;
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
