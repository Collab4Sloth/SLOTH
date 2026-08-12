/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Sintering test taken from "Programming Phase-field modeling - Biner 2017" (AMR version)
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

#include "./SinteringCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Quiet);

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
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
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

  const int refinement_level = 0;  // number of levels of uniform refinement
  const int nx = 25;
  const int ny = 25;
  const double lx = 50;
  const double ly = 50;
  const std::tuple<int, int, double, double>& tuple_of_dimensions =
      std::make_tuple(nx, ny, lx, ly);  // Number of elements and maximum length in each direction

  SPA spatial_c(mesh_type, order_fe, refinement_level, tuple_of_dimensions, true);
  SPA spatial_mu(spatial_c.get_mesh(), 1);
  SPA spatial_eta1(spatial_c.get_mesh(), 1);
  SPA spatial_eta2(spatial_c.get_mesh(), 1);
  // ##############################
  //     Boundary conditions     //
  // ##############################

  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                     Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  // CahnHilliard
  auto ch_bcs_c = BoundaryConditions<FECollection, DIM>(&spatial_c, boundaries);
  auto ch_bcs_mu = BoundaryConditions<FECollection, DIM>(&spatial_mu, boundaries);
  // allenCahn
  auto ac_bcs_eta1 = BoundaryConditions<FECollection, DIM>(&spatial_eta1, boundaries);
  auto ac_bcs_eta2 = BoundaryConditions<FECollection, DIM>(&spatial_eta2, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     Coefficients    //
  // ####################

  // CahnHilliard
  const double ch_lambda = 5.;
  Coefficient ch_grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradC());
  Coefficient ch_energy(Glossary::FreeEnergy, Scheme::Implicit, Gc());
  Coefficient ch_capillary(Glossary::Capillary, ch_lambda);
  Coefficient ch_mobility(Glossary::Mobility, Scheme::Explicit, D());
  Coefficients ch_coef(ch_energy, ch_capillary, ch_mobility, ch_grad_energy);

  // AllenCahn
  const double ac_mob(10.);
  const double ac_lambda = 2.;
  Coefficient ac_grad_energy(Glossary::GradEnergy, Scheme::Implicit, Gradeta());
  Coefficient ac_energy(Glossary::FreeEnergy, Scheme::Implicit, Geta());
  Coefficient ac_capillary(Glossary::Capillary, ac_lambda);
  Coefficient ac_mobility(Glossary::Mobility, ac_mob);
  Coefficients ac_coef(ac_energy, ac_capillary, ac_mobility, ac_grad_energy);

  // ####################
  //     variables     //
  // ####################

  // CahnHilliard
  auto ch_initial_condition =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx = x[0] - 25;
        double yy1 = x[1] - 20;
        double yy2 = x[1] - 35;
        double rr1 = std::sqrt(xx * xx + yy1 * yy1);
        double rr2 = std::sqrt(xx * xx + yy2 * yy2);

        double sol = 0.0;
        if (rr1 < 10 || rr2 < 5) {
          sol = 1;
        }

        return sol;
      });
  auto ch_initial_condition_mu =
      std::function<double(const mfem::Vector&, double)>([&](const mfem::Vector& x, double time) {
        double phi = ch_initial_condition(x, time);
        double sol = 2 * phi * (1 - phi) * (1 - 2 * phi);

        return sol;
      });

  auto ch_initial_condition_1 = AnalyticalFunctions<DIM>(ch_initial_condition);
  auto mu_initial_condition = AnalyticalFunctions<DIM>(ch_initial_condition_mu);
  auto ch_v1 = VAR(&spatial_c, ch_bcs_c, "c", Glossary::PhaseField, 2, ch_initial_condition_1);
  ch_v1.set_additional_information("c");
  auto ch_v2 =
      VAR(&spatial_mu, ch_bcs_mu, "mu", Glossary::ChemicalPotential, 2, mu_initial_condition);
  ch_v2.set_additional_information("mu");
  auto ch_vars = VARS(ch_v1, ch_v2);

  // AllenCahn

  auto ac_initial_condition_1 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx = x[0] - 25;
        double yy1 = x[1] - 20;
        double yy2 = x[1] - 35;
        double rr1 = std::sqrt(xx * xx + yy1 * yy1);
        double rr2 = std::sqrt(xx * xx + yy2 * yy2);

        double sol = 0.0;
        if (rr1 < 10) {
          sol = 1;
        }

        return sol;
      });
  auto ac_initial_condition_2 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx = x[0] - 25;
        double yy1 = x[1] - 20;
        double yy2 = x[1] - 35;
        double rr1 = std::sqrt(xx * xx + yy1 * yy1);
        double rr2 = std::sqrt(xx * xx + yy2 * yy2);

        double sol = 0.0;
        if (rr2 < 5) {
          sol = 1;
        }
        return sol;
      });

  auto ac_v1 = VAR(&spatial_eta1, ac_bcs_eta1, "eta_1", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_1));
  ac_v1.set_additional_information("eta");
  auto ac_v2 = VAR(&spatial_eta2, ac_bcs_eta2, "eta_2", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_2));
  ac_v2.set_additional_information("eta");
  auto ac_vars = VARS(ac_v1, ac_v2);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const auto& frequency = 1;

  const std::vector<int> iterations_list = {0, 1, 3, 5};
  const std::vector<double> times_list = {0.35, 0.45};

  std::string calculation_path = "CahnHilliard";
  const double threshold = 10.;
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {{"phi", {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  auto ch_p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", "CahnHilliard"),

                 Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
                 Parameter("integral_to_compute", map_threshold_integral),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));

  auto ac_p_pst = Parameters(
      Parameter("main_folder_path", main_folder_path), Parameter("calculation_path", "AllenCahn"),

      Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  // ####################
  //     operators     //
  // ####################

  // CahnHilliard
  std::vector<SPA*> spatials{&spatial_c, &spatial_mu};
  OPE ch_oper(spatials, {"CahnHilliard"}, TimeScheme::EulerImplicit, "SplitTimeDerivative");
  ch_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));
  const auto& ch_solver = HypreSolverType::HYPRE_GMRES;
  const auto& ch_precond = HyprePreconditionerType::HYPRE_ILU;
  ch_oper.overload_solver(ch_solver);
  ch_oper.overload_preconditioner(ch_precond);

  auto ch_pst = PST(&spatial_c, ch_p_pst);
  PB ch_pb("CahnHilliard", ch_oper, ch_vars, {ch_coef, ch_coef}, ch_pst, ac_vars);

  std::vector<SPA*> ac_spatials{&spatial_eta1, &spatial_eta2};
  OPE ac_oper(ac_spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  ac_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_ILU;
  ac_oper.overload_solver(solver);
  ac_oper.overload_preconditioner(precond);
  auto ac_pst = PST(&spatial_eta1, ac_p_pst);
  PB ac_pb("AllenCahn", ac_oper, ac_vars, {ac_coef, ac_coef}, ac_pst, ch_vars);
  // AMR
  /////////////////////////////////
  ///  AC
  mfem::ConstantCoefficient amr_coef_ac{1.0};
  mfem::DiffusionIntegrator amr_integ_ac{amr_coef_ac};
  SlothErrorEstimators estimator_ac(ErrorEstimatorType::KELLY, &amr_integ_ac);

  MultiVariableMaxAMR<VARS> amr_ac(*spatial_eta1.get_mesh(), spatial_eta1.is_nc_simplices());

  auto amr_params = Parameters(Parameter("max_elem_error", 1.e-5), Parameter("amr_max_level", 3),
                               Parameter("nc_limit", 0), Parameter("max_preref_cycles", 3));

  amr_ac.SetCriteria(/*estimator*/ &estimator_ac, amr_params);
  ac_pb.set_amr(&amr_ac);
  /////////////////////////////////
  // AMR
  // Coupling 1
  auto cc = Coupling("CahnHilliard/AllenCahn Coupling", ac_pb, ch_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 2.e-4;  // 2.
  const double dt = 1.e-4;
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
