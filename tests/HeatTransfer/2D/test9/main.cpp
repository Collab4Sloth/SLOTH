/**
 * @file main.cpp
 * @brief 2D advection-diffusion problem (MMS) with a gradient velocity v = alpha grad(a),
 *        a being an auxiliary variable
 * @version 0.1
 * @date 2026-10-19
 *
 * @copyright Copyright (c) 2026
 *
 * rho cp dT/dt + div(v T) - div(k grad T) = f   on [0,1]^2,  rho = cp = k = 1
 *   v = alpha grad(a), alpha = 1, a = x + (x^2 + y^2)/2  ->  v = (1 + x, y)   (div v = 2)
 * Same problem as test7, the velocity being given by the gradient of the auxiliary variable a
 * (VelocityPrefactor coefficient + parameter "velocity_gradient"). a is quadratic: with Q2
 * elements, grad(a_h) = v exactly and the results are those of test7 (order 3). With Q1 elements,
 * grad(a_h) is piecewise constant: the O(h) error on the velocity (in particular the normal
 * velocity on the boundaries) limits the convergence of T to order 1.
 *   T = t (1 + cos(pi x) sin(pi y)) (exact solution)
 * Boundaries (all Neumann for the diffusion, g = k grad(T).n):
 *   left  (x = 0): v.n = -1     -> inflow, T_in = exact solution, g = 0
 *   right (x = 1): v.n = 2      -> outflow,                       g = 0
 *   lower (y = 0): v.n = 0      -> wall,                          g = -pi t cos(pi x)
 *   upper (y = 1): v.n = 1      -> outflow,                       g = -pi t cos(pi x)
 * The Inflow coefficient is attached to all the boundaries: the inflow/outflow switch is
 * done pointwise from the sign of v.n.
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "./Coefficient.hpp"
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
  using namespace Sloth2D;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const int refinement_level = 0;

  std::vector<std::string> vect_elem{"InlineSquareWithQuadrangles"};
  std::vector<int> vect_order{2, 1};
  std::vector<int> vect_NN{80, 40, 20, 10};
  for (const auto& elem_type : vect_elem) {
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

        auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                           Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
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
        Coefficient density(Glossary::Concentration, 1.);
        Coefficient heat_capacity(Glossary::Cp, 1.);
        Coefficient conductivity(Glossary::Conductivity, 1.);
        Coefficient neumann(Glossary::Neumann, Scheme::Implicit, NeumannCoefficient());
        neumann.set_bdr_index_coef(std::vector<int>{0, 2});

        // Advection: v = alpha grad(a) (a: auxiliary variable, see below) and inflow value
        Coefficient alpha(Glossary::VelocityPrefactor, 1.);
        auto params = Parameters(Parameter("velocity_gradient", "a"));
        Coefficient inflow(Glossary::Inflow, Scheme::Explicit, InflowCoefficient());
        inflow.set_bdr_index_coef(std::vector<int>{0, 1, 2, 3});

        // Manufactured source term
        Coefficient source(Glossary::SourceTerm, Scheme::Explicit, SourceCoefficient());

        // ####################
        //     variables     //
        // ####################

        auto user_func = std::function<double(const mfem::Vector&, double)>(
            [](const mfem::Vector& x, double time) {
              const auto xx = x[0];
              const auto yy = x[1];
              return time * (1. + std::cos(M_PI * xx) * std::sin(M_PI * yy));
            });
        auto T_analytical = AnalyticalFunctions<DIM>(user_func);

        auto heat_vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 0, T_analytical));

        // Auxiliary variable a = x + (x^2 + y^2)/2 (not solved, grad(a) = (1 + x, y))
        auto a_func = std::function<double(const mfem::Vector&, double)>(
            [](const mfem::Vector& x, [[maybe_unused]] double time) {
              return x[0] + 0.5 * (x[0] * x[0] + x[1] * x[1]);
            });
        auto a_initial_condition = AnalyticalFunctions<DIM>(a_func);
        auto a_var = VAR(&spatial, bcs, "a", Glossary::Temperature, 2, a_initial_condition);
        a_var.set_additional_information("a");
        auto aux_vars = VARS(a_var);

        // ###########################################
        // ###########################################
        //      Post-processing                     //
        // ###########################################
        // ###########################################
        const std::string& main_folder_path =
            "Saves_order_" + std::to_string(order_fe) + "_Nx" + std::to_string(NN);
        const auto& level_of_detail = 1;
        const auto& frequency = 1;
        const std::string& calculation_path = "Problem1";
        auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                                Parameter("calculation_path", calculation_path),
                                Parameter("frequency", frequency),
                                Parameter("level_of_detail", level_of_detail),
                                Parameter("enable_compute_energies", false));
        auto pst = PST(&spatial, p_pst);

        // ####################
        //     Problems      //
        // ####################
        Coefficients coef_pb(density, heat_capacity, conductivity, neumann, alpha, inflow, source);
        std::vector<SPA*> spatials{&spatial};
        TransientOPE oper(spatials, {"Fourier", "Advection", "SourceTerm"}, params,
                          TimeScheme::EulerImplicit, "HeatTimeDerivative");

        oper.overload_nl_solver(
            NLSolverType::NEWTON,
            Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                       Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
        TransientPB Heat_pb("Heat", oper, heat_vars, {coef_pb}, pst, aux_vars);

        // Coupling 1
        auto cc = Coupling("Advection-diffusion", Heat_pb);
        // ###########################################
        // ###########################################
        //            Time-integration              //
        // ###########################################
        // ###########################################
        // T is linear in time: one implicit Euler step is exact in time
        const auto& t_initial = 0.0;
        const auto& t_final = 1.;
        const auto& dt = 1.;
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
