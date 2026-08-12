/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief 2D spinodal decomposition solved by Cahn-Hilliard equations (AMR version)
 * @version 0.1
 * @date 2026-08-11
 *
 * Copyright CEA (c) 2026
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "./CahnHilliardCoefficients.hpp"
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
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  using PST = Test<DIM>::PST;
  /////////////////////////

  using OPE = TransientOperator<FECollection, DIM>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string mesh_type =
      "InlineSquareWithQuadrangles";  // type of mesh // "InlineSquareWithTriangles"
  const int order_fe = 1;             // finite element order
  const int refinement_level = 0;     // number of levels of uniform refinement
  const int nx = 10;
  const int ny = 10;
  const double lx = 200;
  const double ly = 200;
  const std::tuple<int, int, double, double>& tuple_of_dimensions =
      std::make_tuple(nx, ny, lx, ly);  // Number of elements and maximum length in each direction

  SPA spatial_phi(mesh_type, order_fe, refinement_level, tuple_of_dimensions, true);
  SPA spatial_mu(spatial_phi.get_mesh(), order_fe);  // fespace_ indépendant, mesh_ partagé
  SPA spatial_mpi(spatial_phi.get_mesh(), order_fe);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                     Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  auto bcs_phi = BCS(&spatial_phi, boundaries);
  auto boundaries_mu = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                        Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  auto bcs_mu = BCS(&spatial_mu, boundaries_mu);
  auto boundaries_mpi = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                         Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  auto bcs_mpi = BCS(&spatial_mpi, boundaries_mu);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  const double mob(5.);
  const double lambda = 2.;
  auto params = Parameters(Parameter("lambda", lambda));
  // ####################
  //     coefficients  //
  // ####################

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradEnergy());
  Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, DoubleWell());
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  // ####################
  //     variables     //
  // ####################

  auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& x, [[maybe_unused]] double time) {
        double co = 0.5;
        double epsilon = 0.01;
        double xx = x[0];
        double yy = x[1];

        double sol =
            co + epsilon * (std::cos(0.105 * xx) * std::cos(0.11 * yy) +
                            (std::cos(0.13 * xx) * std::cos(0.087 * yy)) *
                                (std::cos(0.13 * xx) * std::cos(0.087 * yy)) +
                            (std::cos(0.025 * xx - 0.15 * yy) * std::cos(0.07 * xx - 0.02 * yy)));

        return sol;
      });
  auto ch_initial_condition_mu =
      std::function<double(const mfem::Vector&, double)>([&](const mfem::Vector& x, double time) {
        double phi = user_func_solution(x, time);
        double sol = 2 * phi * (1 - phi) * (1 - 2 * phi);

        return sol;
      });

  auto phi_initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto mu_initial_condition = AnalyticalFunctions<DIM>(ch_initial_condition_mu);
  const std::string& var_name_1 = "phi";
  const std::string& var_name_2 = "mu";
  auto v1 = VAR(&spatial_phi, bcs_phi, var_name_1, Glossary::PhaseField, 2, phi_initial_condition);
  auto v2 =
      VAR(&spatial_mu, bcs_mu, var_name_2, Glossary::ChemicalPotential, 2, mu_initial_condition);
  auto vars = VARS(v1, v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const int frequency = 100;
  std::string calculation_path = "CahnHilliard";
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
      {var_name_1, {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail),
                 Parameter("integral_to_compute", map_threshold_integral),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  Coefficients coef_pb1(double_well, capillary, mobility, grad_energy);
  std::vector<SPA*> spatials{&spatial_phi, &spatial_mu};
  OPE oper(spatials, {"CahnHilliard"}, params, TimeScheme::EulerImplicit, "SplitTimeDerivative");
  oper.overload_nl_solver(NLSolverType::NEWTON,
                          Parameters(Parameter("description", "Newton solver "),
                                     Parameter("print_level", -1), Parameter("rel_tol", 1.e-12),
                                     Parameter("abs_tol", 1.e-12), Parameter("iter_max", 1000)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_ILU;
  oper.overload_solver(solver);
  oper.overload_preconditioner(precond);

  auto pst = PST(&spatial_phi, p_pst);
  PB problem1(oper, vars, {coef_pb1, coef_pb1}, pst);

  // AMR
  mfem::ConstantCoefficient amr_coef{1.0};
  mfem::DiffusionIntegrator amr_integ{amr_coef};
  SlothErrorEstimators estimator(ErrorEstimatorType::KELLY, &amr_integ);

  SingleVariableAMR<VARS> amr(*spatial_phi.get_mesh(), spatial_phi.is_nc_simplices(), 0);

  auto amr_params = Parameters(Parameter("max_elem_error", 1.e-4), Parameter("amr_max_level", 4),
                               Parameter("nc_limit", 0), Parameter("max_preref_cycles", 4));
  amr.SetCriteria(/*estimator*/ &estimator, amr_params);
  problem1.set_amr(&amr);

  // Coupling 1
  auto cc = Coupling("CahnHilliard Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 3.;  // 5.e4;
  const double dt = 1.;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

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
