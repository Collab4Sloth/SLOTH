/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief MMS convergence test in time for a partitioned coupling Allen-Cahn / Cahn-Hilliard
 *        (implicit Euler vs BDF2), with the structure of Studies/ostwald and Studies/sintering:
 *        each problem sees the variable of the other one as an auxiliary variable.
 *
 * Free energy f(c, eta) = c^2 (1-c)^2 + eta^2 (1-eta)^2 + K c^2 eta^2   (K = -1)
 *   Cahn-Hilliard (c, mu):  dc/dt = Delta mu + S_c,   mu = df/dc - Delta c   (eta: auxiliary)
 *   Allen-Cahn (eta):        deta/dt = -(df/deta - Delta eta) + S_eta        (c: auxiliary)
 *
 * Manufactured solution on [0, 1]^2 (Dirichlet on x = 0 and x = 1, Neumann on y = 0, y = 1):
 *   c = (1 + sin(2 pi t)/2) sin(pi x),   eta = (1 - cos(2 pi t)/2) sin(pi x).
 *
 * Coupling("...", ac_pb, ch_pb): Allen-Cahn is solved first and sees c^n (lag of one step),
 * Cahn-Hilliard is solved second and sees eta^{n+1} (Gauss-Seidel, no coupling iteration).
 * The lag of c is a splitting error O(dt).
 * Without treatment, BDF2 is only of order ~1.
 * With BDF2, the lagged auxiliary variable is extrapolated to t^{n+1}, (1 + w) c^n - w c^{n-1}.
 * Expected implicit Euler order 1, BDF2 order 2 on both eta and c.
 *
 * Fixed mesh: Q2, 160 x 2 elements (solution independent of y), spatial error below the
 * smallest BDF2 error. Time steps from 0.05 to 0.00625.
 *
 * @version 0.1
 * @date 2026-10-29
 *
 * Copyright CEA (c) 2026
 */
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "./CoupledCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Normal);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  using namespace Sloth2D;

  // ##############################
  //  Fixed spatial discretization
  // ##############################
  const std::string mesh_type = "InlineSquareWithQuadrangles";
  const int order_fe = 2;
  const int refinement_level = 0;
  const int nx = 160;  // elements in x
  const int ny = 2;    // elements in y (solution independent of y)
  const double lx = 1.;
  const double ly = 1.;
  const double K = -1.;  // coupling constant (also in coefficients.json)

  // ##############################
  //  Time schemes and time steps
  // ##############################
  const std::vector<std::pair<std::string, TimeScheme::value>> vect_scheme{
      {"EulerImplicit", TimeScheme::EulerImplicit}, {"BDF2", TimeScheme::BDF2}};
  const std::vector<double> vect_dt{0.05, 0.025, 0.0125, 0.00625};

  // ##############################
  //  Manufactured solution
  // ##############################
  const auto g1 = [](double t) { return 1. + 0.5 * std::sin(2. * M_PI * t); };
  const auto g2 = [](double t) { return 1. - 0.5 * std::cos(2. * M_PI * t); };
  auto c_solution = std::function<double(const mfem::Vector&, double)>(
      [g1](const mfem::Vector& x, double t) { return g1(t) * std::sin(M_PI * x[0]); });
  auto eta_solution = std::function<double(const mfem::Vector&, double)>(
      [g2](const mfem::Vector& x, double t) { return g2(t) * std::sin(M_PI * x[0]); });
  // mu = df/dc - Delta c = 2c(1-c)(1-2c) + 2K c eta^2 + pi^2 c
  auto mu_solution = std::function<double(const mfem::Vector&, double)>(
      [g1, g2, K](const mfem::Vector& x, double t) {
        const double c = g1(t) * std::sin(M_PI * x[0]);
        const double eta = g2(t) * std::sin(M_PI * x[0]);
        return 2. * c * (1. - c) * (1. - 2. * c) + 2. * K * c * eta * eta + M_PI * M_PI * c;
      });

  for (const auto& [scheme_name, scheme] : vect_scheme) {
    for (const auto dt : vect_dt) {
      Profiling::getInstance().enable();
      SPA spatial(mesh_type, order_fe, refinement_level, std::make_tuple(nx, ny, lx, ly));

      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 0.),
                         Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 0.)};
      auto ch_bcs = BCS(&spatial, boundaries);
      auto ac_bcs = BCS(&spatial, boundaries);

      // ##############################
      //     Coefficients            //
      // ##############################
      const double lambda = 1.;
      const double mob = 1.;
      // Cahn-Hilliard
      Coefficient ch_grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
      Coefficient ch_energy(Glossary::FreeEnergy, Scheme::Implicit, Gc());
      Coefficient ch_capillary(Glossary::Capillary, lambda);
      Coefficient ch_mobility(Glossary::Mobility, mob);
      Coefficient ch_no_source(Glossary::SourceTerm, 0.);
      Coefficient ch_source(Glossary::SourceTerm, Scheme::Explicit, SourceC());
      Coefficients ch_coef0(ch_energy, ch_capillary, ch_mobility, ch_grad_energy, ch_no_source);
      Coefficients ch_coef1(ch_energy, ch_capillary, ch_mobility, ch_grad_energy, ch_source);
      // Allen-Cahn
      Coefficient ac_grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
      Coefficient ac_energy(Glossary::FreeEnergy, Scheme::Implicit, Geta());
      Coefficient ac_capillary(Glossary::Capillary, lambda);
      Coefficient ac_mobility(Glossary::Mobility, mob);
      Coefficient ac_source(Glossary::SourceTerm, Scheme::Explicit, SourceEta());
      Coefficients ac_coef(ac_energy, ac_capillary, ac_mobility, ac_grad_energy, ac_source);

      // ##############################
      //     Variables               //
      // ##############################
      auto ch_v1 = VAR(&spatial, ch_bcs, "c", Glossary::PhaseField, 2,
                       AnalyticalFunctions<DIM>(c_solution), AnalyticalFunctions<DIM>(c_solution));
      ch_v1.set_additional_information("c");
      // Bounds of the extrapolated state (BDF2, c lagged for Allen-Cahn): they contain the exact
      // solution (c in [0, 1.5]), so the results must be unchanged
      ch_v1.set_extrapolation_bounds(0., 2.);
      auto ch_v2 = VAR(&spatial, ch_bcs, "mu", Glossary::ChemicalPotential, 2,
                       AnalyticalFunctions<DIM>(mu_solution));
      ch_v2.set_additional_information("mu");
      auto ch_vars = VARS(ch_v1, ch_v2);

      auto ac_v1 =
          VAR(&spatial, ac_bcs, "eta", Glossary::PhaseField, 2,
              AnalyticalFunctions<DIM>(eta_solution), AnalyticalFunctions<DIM>(eta_solution));
      ac_v1.set_additional_information("eta");
      auto ac_vars = VARS(ac_v1);

      // ##############################
      //     Post-processing         //
      // ##############################
      // Folder parsed by time_convergence.py: Saves_<scheme>_dt_<dt>/<problem>
      std::ostringstream dt_str;
      dt_str << std::setprecision(10) << dt;
      const std::string main_folder_path = "Saves_" + scheme_name + "_dt_" + dt_str.str();
      const int level_of_detail = 1;
      const int frequency = 1000000;  // no field output, only the specialized values
      auto ch_p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                                 Parameter("calculation_path", "CahnHilliard"),
                                 Parameter("frequency", frequency),
                                 Parameter("level_of_detail", level_of_detail),
                                 Parameter("enable_save_specialized_at_iter", true));
      auto ac_p_pst =
          Parameters(Parameter("main_folder_path", main_folder_path),
                     Parameter("calculation_path", "AllenCahn"), Parameter("frequency", frequency),
                     Parameter("level_of_detail", level_of_detail),
                     Parameter("enable_save_specialized_at_iter", true));

      // ##############################
      //     Operators and problems  //
      // ##############################
      const auto newton =
          Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                     Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-10));
      // Cahn-Hilliard: eta as auxiliary variable
      std::vector<SPA*> ch_spatials{&spatial, &spatial};
      TransientOPE ch_oper(ch_spatials, {"CahnHilliard", "SourceTerm"}, scheme,
                           "SplitTimeDerivative");
      ch_oper.overload_nl_solver(NLSolverType::NEWTON, newton);
      ch_oper.overload_solver(HypreSolverType::HYPRE_GMRES,
                              Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
      ch_oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU,
                                      Parameters(Parameter("type", 1)));
      auto ch_pst = PST(&spatial, ch_p_pst);
      TransientPB ch_pb("CahnHilliard", ch_oper, ch_vars, {ch_coef0, ch_coef1}, ch_pst, ac_vars);

      // Allen-Cahn: c as auxiliary variable
      std::vector<SPA*> ac_spatials{&spatial};
      TransientOPE ac_oper(ac_spatials, {"AllenCahn", "SourceTerm"}, scheme, "TimeDerivative");
      ac_oper.overload_nl_solver(NLSolverType::NEWTON, newton);
      ac_oper.overload_solver(HypreSolverType::HYPRE_GMRES,
                              Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
      ac_oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU,
                                      Parameters(Parameter("type", 1)));
      auto ac_pst = PST(&spatial, ac_p_pst);
      TransientPB ac_pb("AllenCahn", ac_oper, ac_vars, {ac_coef}, ac_pst, ch_vars);

      // Partitioned coupling: Allen-Cahn first, then Cahn-Hilliard
      auto cc = Coupling("AllenCahn/CahnHilliard coupling", ac_pb, ch_pb);

      // ##############################
      //     Time integration        //
      // ##############################
      auto time_params = Parameters(Parameter("initial_time", 0.), Parameter("final_time", 1.),
                                    Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);
      time.solve();
      Profiling::getInstance().print();
    }
  }

  MPI_Finalize();
  return 0;
}
