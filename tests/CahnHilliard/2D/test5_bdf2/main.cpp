/**
 * Copyright CEA (c) 2026
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief MMS convergence test IN TIME: implicit Euler vs BDF2.
 *
 * Variant of CahnHilliard/2D/test5 (Zhang et al., JCP 236, 2013) with a non-polynomial time
 * dependence, so that the time discretization error does not vanish:
 *   phi(x, y, t) = exp(t) sin(pi x),   mu = phi^3 - phi - lambda Delta phi
 * (with (t + 1) instead of exp(t), as in test5, both schemes are exact in time).
 *
 * The solution depends on x only (Neumann conditions on y = 0 and y = 1): the mesh is fixed with
 * 160 Q2 elements in x and only 2 in y. Q2 converges at order 3 in L2 (test5 reference), so the
 * spatial error stays well below the smallest BDF2 error.
 * The time step is divided by 2 from 0.1 to 0.0125.
 * Expected slopes of the L2 error on phi at t = 1: EulerImplicit: 1, BDF2: 2 (see
 * tests/tools/time_convergence.py).
 *
 * @version 0.1
 * @date 2026-10-29
 */
#include <iomanip>
#include <sstream>

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
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  using namespace Sloth2D;

  // ##############################
  //  Fixed spatial discretization
  // ##############################
  const int refinement_level = 0;
  const std::string elem_type = "InlineSquareWithQuadrangles";
  const int order_fe = 2;
  const int NN = 160;  // elements in x
  const int NY = 2;    // elements in y (solution independent of y)

  // ##############################
  //  Time schemes and time steps
  // ##############################
  const std::vector<std::pair<std::string, TimeScheme::value>> vect_scheme{
      {"EulerImplicit", TimeScheme::EulerImplicit}, {"BDF2", TimeScheme::BDF2}};
  const std::vector<double> vect_dt{0.1, 0.05, 0.025, 0.0125};

  for (const auto& [scheme_name, scheme] : vect_scheme) {
    for (const auto dt : vect_dt) {
      const double lx = 1.;
      const double ly = 1.;
      const std::tuple<int, int, double, double>& tuple_of_dimensions =
          std::make_tuple(NN, NY, lx, ly);
      SPA spatial(elem_type, order_fe, refinement_level, tuple_of_dimensions);

      // ##############################
      //     Boundary conditions     //
      // ##############################
      // phi = mu = 0 on x = 0 and x = 1 for all t (sin(pi x) = 0): time-independent conditions
      auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 0.),
                         Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
      auto bcs_phi = BCS(&spatial, boundaries);
      auto boundaries_mu = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 0.),
                            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
      auto bcs_mu = BCS(&spatial, boundaries_mu);

      // ####################
      //     parameters    //
      // ####################
      const double epsilon(1.0);
      const double sigma(1.);
      const double mob(1.);
      const double lambda = 1.;
      const double omega = 1.;
      auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                               Parameter("lambda", lambda), Parameter("omega", omega));

      // ####################
      //     coefficients  //
      // ####################
      Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
      Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, Fw(omega));
      Coefficient capillary(Glossary::Capillary, lambda);
      Coefficient mobility(Glossary::Mobility, mob);

      // ####################
      //     variables     //
      // ####################
      // phi analytical solution (mms)
      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, double time) {
            return std::exp(time) * std::sin(M_PI * x[0]);
          });

      // mu analytical solution (mms): mu = c^3 - c - Delta c, Delta c = -pi^2 c
      auto mu_user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, double time) {
            const double c = std::exp(time) * std::sin(M_PI * x[0]);
            const double delta_c = -M_PI * M_PI * c;
            return c * c * c - c - delta_c;
          });

      auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto phi_analytical_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto mu_initial_condition = AnalyticalFunctions<DIM>(mu_user_func_solution);
      const std::string& var_name_1 = "phi";
      const std::string& var_name_2 = "mu";
      auto v1 = VAR(&spatial, bcs_phi, var_name_1, Glossary::PhaseField, 2, phi_initial_condition,
                    phi_analytical_condition);
      auto v2 =
          VAR(&spatial, bcs_mu, var_name_2, Glossary::ChemicalPotential, 2, mu_initial_condition);
      auto vars = VARS(v1, v2);

      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // Folder parsed by time_convergence.py: Saves_<scheme>_dt_<dt>
      std::ostringstream dt_str;
      dt_str << std::setprecision(10) << dt;
      const std::string& main_folder_path = "Saves_" + scheme_name + "_dt_" + dt_str.str();

      const int level_of_detail = 1;
      const int frequency = 1000000;  // no field output, only the specialized values
      std::string calculation_path = "Problem1";
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
      Coefficient no_source(Glossary::SourceTerm, 0.);
      // Source term evaluated at the time of the step (t^{n+1} for both schemes)
      Coefficient source(Glossary::SourceTerm, Scheme::Explicit, SourceCoefficient());
      Coefficients coef_pb0(double_well, capillary, mobility, grad_energy, no_source);
      Coefficients coef_pb1(double_well, capillary, mobility, grad_energy, source);
      std::vector<SPA*> spatials{&spatial, &spatial};

      TransientOPE oper(spatials, {"CahnHilliard", "SourceTerm"}, params, scheme,
                        "SplitTimeDerivative");

      oper.overload_nl_solver(
          NLSolverType::NEWTON,
          Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                     Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
      const auto& solver = HypreSolverType::HYPRE_GMRES;
      const auto& precond = HyprePreconditionerType::HYPRE_ILU;
      oper.overload_solver(solver, Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
      oper.overload_preconditioner(precond, Parameters(Parameter("type", 1)));

      auto pst = PST(&spatial, p_pst);
      TransientPB problem1(oper, vars, {coef_pb0, coef_pb1}, pst);

      auto cc = Coupling("CahnHilliard Coupling", problem1);

      // ###########################################
      //            Time-integration              //
      // ###########################################
      const double t_initial = 0.0;
      const double t_final = 1.;
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
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  return 0;
}
