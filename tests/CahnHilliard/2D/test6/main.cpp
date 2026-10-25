/**
 * @file main.cpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Rising bubble advected by an imposed Poiseuille flow (Cahn-Hilliard + SplitAdvection),
 *        AMR version
 * @version 0.1
 * @date 2026-10-18
 *
 * Copyright CEA (c) 2026
 *
 * Test case of T. Schnetz's internship, adapted from C. Introïni's PhD:
 *
 *   dphi/dt = div( M grad(mu) ) - div( phi v )
 *   mu      = omega W'(phi) - lambda lap(phi),  W(phi) = phi^2 (1 - phi)^2
 *
 * Domain [0, R] x [0, H], R = 5 rb, H = 8 rb. Bubble (phi = 0) of radius rb centred at
 * (R/2, 3 rb) in a matrix (phi = 1). Imposed velocity:
 *
 *   v = vmax (1 - ((x - R/2)/(R/2))^2) e_y,  vmax = 0.05
 *
 * (v = 0 on the lateral walls, div(v) = 0). Two configurations:
 *  - "periodic": periodic in x and y, the bubble re-enters from the bottom,
 *  - "open"    : periodic in x, inflow (phi = 1) at the bottom, free outflow at the top
 *                (Inflow coefficient on boundaries 0 and 2, AdvectionBoundaryCH).
 * In both cases the integral of phi is saved: constant for "periodic", it increases once the
 * bubble leaves the domain for "open" (phi = 1 replaces phi = 0).
 *
 * An MPI_Problem is coupled to the Cahn-Hilliard problem to save the MPI rank field
 * (Saves_<config>/ProblemMPI) and check the load balancing along the AMR passes
 * (mesh rebalanced after each AMR pass: "rebalance" AMR parameter).
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
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
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  /////////////////////////
  using namespace Sloth2D;

  // ##############################
  //     Geometry and physics    //
  // ##############################
  // Physical data (rising bubble, Bo = 10, Mo = 1e-2)
  const double sigma = 7.e-2;  // interfacial energy [N/m]
  const double rho1 = 1.;      // gas density [kg/m3]
  const double rho2 = 1.e3;    // liquid density [kg/m3]
  const double gravity = 9.81;
  const double bond = 10.;
  const double rb = 0.5 * std::sqrt(sigma * bond / (gravity * (rho2 - rho1)));  // 4.226e-3 m
  const double R = 5. * rb;   // width (must be consistent with xc in coefficient.json)
  const double H = 8. * rb;   // height
  const double xb = 0.5 * R;  // bubble centre
  const double yb = 3. * rb;

  const double epsilon = rb / 7.;  // interface thickness
  const double mob = 1.e-7;        // mobility
  const double lambda = 1.5 * sigma * epsilon;
  const double omega = 12. * sigma / epsilon;

  const std::vector<std::string> vect_config{"periodic", "open"};
  for (const auto& config : vect_config) {
    const bool open = (config == "open");
    // ###########################################
    // ###########################################
    //         Spatial Discretization           //
    // ###########################################
    // ###########################################
    // ##############################
    //           Meshing           //
    // ##############################
    const int order_fe = 1;
    const int refinement_level = 0;
    const int nx = 15;  // coarse mesh, refined by AMR in the interface
    const int ny = 24;
    const auto tuple_of_dimensions = std::make_tuple(nx, ny, R, H);

    std::vector<mfem::Vector> translations;
    translations.emplace_back(mfem::Vector({R, 0.0}));
    if (!open) {
      translations.emplace_back(mfem::Vector({0.0, H}));
    }
    SPA spatial_phi("InlineSquareWithQuadrangles", order_fe, refinement_level, tuple_of_dimensions,
                    translations, true);
    SPA spatial_mu(spatial_phi.get_mesh(), order_fe, true);  // own fespace, shared periodic mesh
    // ##############################
    //     Boundary conditions     //
    // ##############################
    const std::string bc_y = open ? "Neumann" : "Periodic";
    auto boundaries_phi = {Boundary("lower", 0, bc_y), Boundary("right", 1, "Periodic"),
                           Boundary("upper", 2, bc_y), Boundary("left", 3, "Periodic")};
    auto boundaries_mu = {Boundary("lower", 0, bc_y), Boundary("right", 1, "Periodic"),
                          Boundary("upper", 2, bc_y), Boundary("left", 3, "Periodic")};
    auto bcs_phi = BCS(&spatial_phi, boundaries_phi);
    auto bcs_mu = BCS(&spatial_mu, boundaries_mu);

    // ###########################################
    // ###########################################
    //            Physical models               //
    // ###########################################
    // ###########################################
    // ####################
    //     parameters    //
    // ####################
    auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                             Parameter("lambda", lambda), Parameter("omega", omega));
    // ####################
    //     coefficients  //
    // ####################
    Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
    Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
    Coefficient capillary(Glossary::Capillary, lambda);
    Coefficient mobility(Glossary::Mobility, mob);

    // Imposed velocity, carried by phi (block 0)
    Coefficient velocity_x(Glossary::Velocity_X, 0.);
    Coefficient velocity_y(Glossary::Velocity_Y, Scheme::Explicit, PoiseuilleVelocity());
    // Open configuration: matrix (phi = 1) enters at the bottom, free outflow at the top
    Coefficient inflow(Glossary::Inflow, 1.);
    inflow.set_bdr_index_coef(std::vector<int>{0, 2});

    Coefficients coef_phi =
        open ? Coefficients(double_well, capillary, mobility, grad_energy, velocity_x, velocity_y,
                            inflow)
             : Coefficients(double_well, capillary, mobility, grad_energy, velocity_x, velocity_y);
    Coefficients coef_mu(double_well, capillary, mobility, grad_energy);
    // ####################
    //     variables     //
    // ####################
    auto bubble = std::function<double(const mfem::Vector&, double)>(
        [&](const mfem::Vector& x, [[maybe_unused]] double time) {
          const double r = std::hypot(x[0] - xb, x[1] - yb);
          return 0.5 + 0.5 * std::tanh(2. * (r - rb) / epsilon);
        });
    auto phi_initial_condition = AnalyticalFunctions<DIM>(bubble);
    auto mu_initial_condition =
        AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
            [&](const mfem::Vector& x, [[maybe_unused]] double time) { return 0.; }));

    const std::string& var_name_1 = "phi";
    const std::string& var_name_2 = "mu";
    auto v1 =
        VAR(&spatial_phi, bcs_phi, var_name_1, Glossary::PhaseField, 2, phi_initial_condition);
    // mu = 0: equilibrium tanh profile (omega W'(phi) = lambda phi'' since lambda/omega = eps^2/8)
    auto v2 =
        VAR(&spatial_mu, bcs_mu, var_name_2, Glossary::ChemicalPotential, 2, mu_initial_condition);
    auto vars = VARS(v1, v2);

    // ###########################################
    // ###########################################
    //      Post-processing                     //
    // ###########################################
    // ###########################################
    const std::string& main_folder_path = "Saves_" + config;
    const int level_of_detail = 1;
    const int frequency = 60;  // every 0.06 s
    std::string calculation_path = "CahnHilliard";
    std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
        {var_name_1, {-10.1, 10.1}}};
    bool enable_save_specialized_at_iter = true;
    auto p_pst =
        Parameters(Parameter("main_folder_path", main_folder_path),
                   Parameter("calculation_path", calculation_path),
                   Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
                   Parameter("integral_to_compute", map_threshold_integral),
                   Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
    // ####################
    //     operators     //
    // ####################
    SPAS spatials{&spatial_phi, &spatial_mu};
    TransientOPE oper(spatials, {"CahnHilliard", "SplitAdvection"}, params,
                      TimeScheme::EulerImplicit, "SplitTimeDerivative");
    oper.overload_nl_solver(NLSolverType::NEWTON,
                            Parameters(Parameter("description", "Newton solver "),
                                       Parameter("print_level", 1), Parameter("rel_tol", 1.e-10),
                                       Parameter("abs_tol", 1.e-14), Parameter("iter_max", 50)));
    const auto& solver = HypreSolverType::HYPRE_GMRES;
    const auto& precond = HyprePreconditionerType::HYPRE_ILU;
    oper.overload_solver(solver, Parameters(Parameter("tol", 1.e-12), Parameter("kDim", 100)));
    oper.overload_preconditioner(precond);

    auto pst = PST(&spatial_phi, p_pst);
    TransientPB problem1(oper, vars, {coef_phi, coef_mu}, pst);

    // AMR on phi (Kelly estimator), initial refinement of the bubble before time-stepping.
    // Very low derefinement threshold (1e-8) to preserve the integral of phi
    mfem::ConstantCoefficient amr_coef{1.0};
    mfem::DiffusionIntegrator amr_integ{amr_coef};
    SlothErrorEstimators estimator(ErrorEstimatorType::KELLY, &amr_integ);

    SingleVariableAMR<VARS> amr(*spatial_phi.get_mesh(), spatial_phi.is_nc_simplices(), 0);
    auto amr_params =
        Parameters(Parameter("max_elem_error", 1.e-6), Parameter("amr_max_level", 3),
                   Parameter("nc_limit", 1), Parameter("max_preref_cycles", 3),
                   Parameter("scale_down_factor", 1.e-2), Parameter("rebalance", true));
    amr.SetCriteria(&estimator, amr_params);
    problem1.set_amr(&amr);

    // MPI problem: stores the MPI rank owning each dof, to check the load balancing
    // after each AMR pass. The variable lives on the AMR mesh (spatial_phi), so it is
    // resynchronized with the mesh after each (de)refinement and rebalancing. It is not a
    // pilot variable of the AMR: the reference results of problem1 are unchanged.
    auto mpi_init = AnalyticalFunctions<DIM>(std::function<double(const mfem::Vector&, double)>(
        [&](const mfem::Vector& x, [[maybe_unused]] double time) { return 0.; }));
    auto vars_mpi = VARS(VAR(&spatial_phi, bcs_phi, "MPI rank", Glossary::MPI, 2, mpi_init));
    auto p_pst_mpi = Parameters(Parameter("main_folder_path", main_folder_path),
                                Parameter("calculation_path", std::string("ProblemMPI")),
                                Parameter("frequency", frequency),
                                Parameter("level_of_detail", level_of_detail));
    auto pst_mpi = PST(&spatial_phi, p_pst_mpi);
    PB_MPI problem_mpi(vars_mpi, pst_mpi);

    // Coupling 1
    auto cc = Coupling("CahnHilliard Coupling", problem1, problem_mpi);

    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const double t_initial = 0.0;
    const double t_final = 0.06;  // 0.6;
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
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
