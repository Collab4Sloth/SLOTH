/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief MMS convergence test in time for Allen-Cahn (periodic square): implicit Euler vs BDF2.
 *
 * Transient variant of AllenCahn/2D_periodic/test3 (steady MMS):
 *   phi(x, y, t) = g(t) sin^2(2x) cos^2(3y),   g(t) = 1 + sin(2 pi t) / 2,
 *   d(phi)/dt = -M (W'(phi) - lambda Delta phi) + S,   W(phi) = phi^2 (1 - phi)^2,
 * on [0, pi]^2 with periodic conditions, t in [0, 1] (S: see coefficient.json).
 * The oscillating g gives time errors large enough with respect to the spatial error: the
 * Laplacian strongly damps the time errors of this high-frequency profile.
 *
 * Fixed mesh: Q3, 48 x 48 (to use a coarser mesh than with Q2)
 * Time step divided by 2 from 0.1 to 0.0125.
 * Expected slopes of the L2 error on phi at t = 1 :  EulerImplicit: 1, BDF2: 2
 *
 * @version 0.1
 * @date 2026-10-29
 *
 * @copyright Copyright (c) 2026
 */
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>
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
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  using namespace Sloth2D;

  // ##############################
  //  Fixed spatial discretization
  // ##############################
  const int refinement_level = 0;
  const double L = M_PI;
  const int order = 3;
  const int NN = 48;

  // ##############################
  //  Time schemes and time steps
  // ##############################
  const std::vector<std::pair<std::string, TimeScheme::value>> vect_scheme{
      {"EulerImplicit", TimeScheme::EulerImplicit}, {"BDF2", TimeScheme::BDF2}};
  const std::vector<double> vect_dt{0.1, 0.05, 0.025, 0.0125};

  for (const auto& [scheme_name, scheme] : vect_scheme) {
    for (const auto dt : vect_dt) {
      //---------------------------------------
      // Profiling start
      Profiling::getInstance().enable();
      //---------------------------------------
      mfem::Vector x_translation({L, 0.0});
      mfem::Vector y_translation({0.0, L});
      std::vector<mfem::Vector> translations = {x_translation, y_translation};
      SPA spatial("InlineSquareWithQuadrangles", order, refinement_level,
                  std::make_tuple(NN, NN, L, L), translations);
      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                         Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
      auto bcs = BCS(&spatial, boundaries);

      // ####################
      //     parameters    //
      // ####################
      const double mob(1.);
      const double lambda = 1.;
      const double omega = 1.;
      Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
      Coefficient double_well_imp(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
      Coefficient capillary(Glossary::Capillary, lambda);
      Coefficient mobility(Glossary::Mobility, mob);
      // Source term evaluated at the time of the step (t^{n+1} for both schemes)
      Coefficient source(Glossary::SourceTerm, Scheme::Explicit, SourceCoefficient());
      Coefficients coef_ac(double_well_imp, capillary, mobility, grad_energy, source);

      // ####################
      //     variables     //
      // ####################
      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double g = 1. + 0.5 * std::sin(2. * M_PI * time);
            return g * std::pow(std::sin(2 * v[0]), 2) * std::pow(std::cos(3 * v[1]), 2);
          });
      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);

      auto vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, initial_condition,
                           analytical_solution));

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

      auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                              Parameter("calculation_path", calculation_path),
                              Parameter("frequency", frequency),
                              Parameter("level_of_detail", level_of_detail),
                              Parameter("force_clean_output_dir", false));

      // ####################
      //     operators     //
      // ####################
      std::vector<SPA*> spatials{&spatial};
      TransientOPE oper(spatials, {"AllenCahn", "SourceTerm"}, scheme, "TimeDerivative");
      oper.overload_nl_solver(
          NLSolverType::NEWTON,
          Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                     Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
      oper.overload_solver(HypreSolverType::HYPRE_GMRES,
                           Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
      oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU,
                                   Parameters(Parameter("type", 1)));
      auto pst = PST(&spatial, p_pst);
      TransientPB problem1("AllenCahn", oper, vars, {coef_ac}, pst);

      auto cc = Coupling("Default Coupling", problem1);

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
