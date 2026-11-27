/**
 * Copyright CEA (c) 2026
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Dendritic growth of a pure substance: PFHub benchmark 3 (Karma-Rappel model)
 *
 * Phase field phi (-1: liquid, 1: solid) and dimensionless temperature U = (T - Tm) / (L / cp):
 *   tau(n) dphi/dt = - dF/dphi,   dU/dt = D Delta U + 1/2 dphi/dt,
 *   F = int 1/2 (W(n))^2 |grad phi|^2 - phi^2/2 + phi^4/4
 *                + lambda U phi (1 - 2/3 phi^2 + 1/5 phi^4),
 *   W(n) = W0 a(n),  tau(n) = tau0 a(n)^2,  a(n) = 1 + e4 cos(4 theta),  n = grad phi / |grad phi|,
 *   lambda = D tau0 / (0.6267 W0^2).
 * Variational free energy and relaxation time (coefficients.json): the anisotropy is written in
 * the JSON file, with cos(4 theta) = nx^4 - 6 nx^2 ny^2 + ny^4.
 *
 * Partitioned coupling: phase field (AllenCahn + RelaxationTimeDerivative with tau(n) explicit,
 * U auxiliary variable), then temperature (Fourier + LatentHeat with L = 1/2, phi auxiliary
 * variable).
 * One quarter of the dendrite: seed of radius 8 at the origin, U = -0.3 (undercooling), no-flux
 * conditions.
 *
 * AMR driven by phi (Kelly estimator) on a nonconforming mesh shared by phi and U, with
 * rebalancing after each AMR pass (as AllenCahn/2D/test1_amr): coarse mesh refined in the
 * interface, initial refinement of the seed before time-stepping.
 *
 * Default options: reduced domain and time for the non-regression test. Benchmark: square of 960
 * units up to t = 1500, e.g. -L 960 -n 240 -al 3 -t 1500 (see the documentation of the example).
 *
 * @version 0.1
 * @date 2026-11-23
 *
 */
#include <cmath>
#include <cstdlib>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "./DendriteCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

struct TestParameters {
  // Mesh (coarse mesh if AMR)
  double control_length = 40.;
  int control_nb_elements = 10;
  int control_refinement_level = 0;

  // AMR
  bool control_amr = true;
  int control_amr_level = 3;
  double control_amr_max_elem_error = 1.e-3;

  // Physics (the anisotropy and lambda are given in coefficients.json)
  double control_diffusivity = 10.;
  double control_undercooling = -0.3;
  double control_seed_radius = 8.;

  // Non-linear solver
  double control_rel_tol = 1.e-10;
  double control_abs_tol = 1.e-12;
  int control_print_level = -1;

  // Time
  bool control_bdf2 = false;
  double control_final_time = 1.;
  double control_time_step = 0.1;

  // Post-processing
  const char* control_output = "Saves";
  int control_save_frequency = 10;
  int control_specialized_frequency = 1;
  int control_level_of_detail = 1;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.control_length, "-L", "--length", "Size of the square domain.");
  args.AddOption(&p.control_nb_elements, "-n", "--nb_elements",
                 "Number of elements in each direction (coarse mesh if AMR).");
  args.AddOption(&p.control_refinement_level, "-rl", "--refinement_level",
                 "Uniform refinement level of the mesh.");
  args.AddOption(&p.control_amr, "-amr", "--amr", "-no-amr", "--no-amr", "Enables AMR.");
  args.AddOption(&p.control_amr_level, "-al", "--amr_level", "Maximum AMR level.");
  args.AddOption(&p.control_amr_max_elem_error, "-ae", "--amr_max_elem_error",
                 "Maximum element error for AMR.");
  args.AddOption(&p.control_diffusivity, "-D", "--diffusivity", "Thermal diffusivity D.");
  args.AddOption(&p.control_undercooling, "-U0", "--undercooling", "Initial undercooling.");
  args.AddOption(&p.control_seed_radius, "-r0", "--seed_radius", "Radius of the seed.");
  args.AddOption(&p.control_rel_tol, "-rtol", "--rel_tol", "Relative tolerance (Newton).");
  args.AddOption(&p.control_abs_tol, "-atol", "--abs_tol", "Absolute tolerance (Newton).");
  args.AddOption(&p.control_print_level, "-pl", "--print_level",
                 "Print level of the non-linear solver (-1: quiet).");
  args.AddOption(&p.control_bdf2, "-bdf2", "--bdf2", "-euler", "--euler",
                 "Time scheme of both problems: BDF2 or implicit Euler.");
  args.AddOption(&p.control_final_time, "-t", "--final_time", "Final time of the simulation.");
  args.AddOption(&p.control_time_step, "-dt", "--time_step",
                 "Constant time-step of the simulation.");
  args.AddOption(&p.control_output, "-o", "--output", "Output folder.");
  args.AddOption(&p.control_save_frequency, "-f", "--save_frequency", "Save frequency.");
  args.AddOption(&p.control_specialized_frequency, "-sf", "--specialized_frequency",
                 "Frequency of the specialized values (energies, integral, isovalue).");
  args.AddOption(&p.control_level_of_detail, "-lod", "--level_of_detail",
                 "Level of detail of the outputs.");
  args.Parse();
  if (!args.Good()) {
    if (mfem::Mpi::WorldRank() == 0) args.PrintUsage(mfem::out);
    std::exit(EXIT_FAILURE);
  }
  if (mfem::Mpi::WorldRank() == 0) args.PrintOptions(mfem::out);
}

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Verbose);
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);

  using namespace Sloth2D;
  {  // scope: every SLOTH/MFEM object is destroyed before MPI_Finalize
    // ###########################################
    //         Spatial Discretization           //
    // ###########################################
    const std::string mesh_type = "InlineSquareWithQuadrangles";
    const int order_fe = 1;
    // Nonconforming mesh (required for derefinement and rebalancing), shared by phi and U
    SPA spatial_phi(mesh_type, order_fe, p.control_refinement_level,
                    std::make_tuple(p.control_nb_elements, p.control_nb_elements, p.control_length,
                                    p.control_length),
                    true);
    SPA spatial_u(spatial_phi.get_mesh(), order_fe);  // own finite element space, shared mesh

    // No-flux conditions on all the boundaries
    auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                       Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
    auto bcs_phi = BCS(&spatial_phi, boundaries);
    auto bcs_u = BCS(&spatial_u, boundaries);

    // ###########################################
    //            Variables                     //
    // ###########################################
    // Quarter of a solid seed at the origin, interface phi = 0 at r = r0
    const double r0 = p.control_seed_radius;
    auto seed = std::function<double(const mfem::Vector&, double)>(
        [r0](const mfem::Vector& x, [[maybe_unused]] double time) {
          const double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
          return -std::tanh((r - r0) / std::sqrt(2.));
        });
    auto phi =
        VAR(&spatial_phi, bcs_phi, "phi", Glossary::PhaseField, 2, AnalyticalFunctions<DIM>(seed));
    phi.set_additional_information("phi");
    auto pf_vars = VARS(phi);

    // Uniform undercooling. Initial conditions given as AnalyticalFunctions: projected again on
    // the mesh after each pre-refinement cycle of the AMR
    const double u0 = p.control_undercooling;
    auto undercooling = std::function<double(const mfem::Vector&, double)>(
        [u0]([[maybe_unused]] const mfem::Vector& x, [[maybe_unused]] double time) { return u0; });
    auto u = VAR(&spatial_u, bcs_u, "U", Glossary::Temperature, 2,
                 AnalyticalFunctions<DIM>(undercooling));
    u.set_additional_information("U");
    auto heat_vars = VARS(u);

    // ###########################################
    //            Coefficients                  //
    // ###########################################
    // Phase field: variational free energy f(phi, grad phi; U) and relaxation time tau(n).
    // tau(n) is explicit (normal of the previous step): its derivative with respect to grad phi
    // behaves as 1 / |grad phi| far from the interface, where grad phi vanishes, and makes the
    // Jacobian too ill-conditioned for the linear solver.
    Coefficient free_energy(Glossary::FreeEnergy, Scheme::Implicit, DendriteFreeEnergy());
    Coefficient relaxation_time(Glossary::RelaxationTime, Scheme::Explicit,
                                AnisotropicRelaxationTime());
    Coefficient mobility(Glossary::Mobility, 1.);
    Coefficients pf_coef(free_energy, mobility, relaxation_time);

    // Temperature: dU/dt = D Delta U + 1/2 dphi/dt
    Coefficient density(Glossary::Concentration, 1.);
    Coefficient heat_capacity(Glossary::Cp, 1.);
    Coefficient conductivity(Glossary::Conductivity, p.control_diffusivity);
    Coefficient latent_heat(Glossary::LatentHeat, 0.5);
    Coefficients heat_coef(density, heat_capacity, conductivity, latent_heat);

    // ###########################################
    //      Post-processing                     //
    // ###########################################
    const std::string main_folder_path = p.control_output;
    // Solid fraction: integral of phi; free energy: Density; tip position: isovalue phi = 0
    std::map<std::string, std::tuple<double, double>> integral_phi = {{"phi", {-1.1, 1.1}}};
    std::map<std::string, double> iso_phi = {{"phi", 0.}};
    auto pf_pst =
        PST(&spatial_phi,
            Parameters(Parameter("main_folder_path", main_folder_path),
                       Parameter("calculation_path", std::string("PhaseField")),
                       Parameter("frequency", p.control_save_frequency),
                       Parameter("specialized_frequency", p.control_specialized_frequency),
                       Parameter("level_of_detail", p.control_level_of_detail),
                       Parameter("integral_to_compute", integral_phi),
                       Parameter("iso_val_to_compute", iso_phi)));
    auto heat_pst = PST(
        &spatial_u, Parameters(Parameter("main_folder_path", main_folder_path),
                               Parameter("calculation_path", std::string("Temperature")),
                               Parameter("frequency", p.control_save_frequency),
                               Parameter("specialized_frequency", p.control_specialized_frequency),
                               Parameter("level_of_detail", p.control_level_of_detail),
                               Parameter("enable_compute_energies", false)));

    // ###########################################
    //            Problems                      //
    // ###########################################
    std::vector<SPA*> spatials_phi{&spatial_phi};
    std::vector<SPA*> spatials_u{&spatial_u};
    const auto newton = Parameters(
        Parameter("description", "Newton solver "), Parameter("print_level", p.control_print_level),
        Parameter("rel_tol", p.control_rel_tol), Parameter("abs_tol", p.control_abs_tol));

    // Same time scheme for both problems (BDF2: U extrapolated to t^{n+1} in the phase-field
    // problem; the latent heat remains first order in time)
    const auto time_scheme = p.control_bdf2 ? TimeScheme::BDF2 : TimeScheme::EulerImplicit;
    TransientOPE pf_oper(spatials_phi, {"AllenCahn"}, time_scheme, "RelaxationTimeDerivative");
    pf_oper.overload_nl_solver(NLSolverType::NEWTON, newton);
    pf_oper.overload_solver(HypreSolverType::HYPRE_GMRES);
    pf_oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU);
    TransientPB pf_pb("PhaseField", pf_oper, pf_vars, {pf_coef}, pf_pst, heat_vars);

    // AMR on phi (Kelly estimator) with rebalancing after each AMR pass; U follows the mesh
    mfem::ConstantCoefficient amr_coef{1.0};
    mfem::DiffusionIntegrator amr_integ{amr_coef};
    SlothErrorEstimators estimator(ErrorEstimatorType::KELLY, &amr_integ);
    SingleVariableAMR<VARS> amr(*spatial_phi.get_mesh(), spatial_phi.is_nc_simplices(), "phi");
    auto amr_params =
        Parameters(Parameter("max_elem_error", p.control_amr_max_elem_error),
                   Parameter("amr_max_level", p.control_amr_level), Parameter("nc_limit", 1),
                   Parameter("max_preref_cycles", p.control_amr_level),
                   Parameter("scale_down_factor", 0.25), Parameter("rebalance", true));
    amr.SetCriteria(&estimator, amr_params);
    if (p.control_amr) pf_pb.set_amr(&amr);

    TransientOPE heat_oper(spatials_u, {"Fourier", "LatentHeat"}, time_scheme,
                           "HeatTimeDerivative");
    heat_oper.overload_nl_solver(NLSolverType::NEWTON, newton);
    TransientPB heat_pb("Temperature", heat_oper, heat_vars, {heat_coef}, heat_pst, pf_vars);

    // Phase field first, then temperature (latent heat of the current step)
    auto cc = Coupling("Dendritic growth", pf_pb, heat_pb);

    // ###########################################
    //            Time-integration              //
    // ###########################################
    auto time_params =
        Parameters(Parameter("initial_time", 0.), Parameter("final_time", p.control_final_time),
                   Parameter("time_step", p.control_time_step));
    auto time = TimeDiscretization(time_params, cc);
    time.solve();

    //---------------------------------------
    // Profiling stop
    //---------------------------------------
    Profiling::getInstance().print();
  }
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  return 0;
}
