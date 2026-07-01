/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Ostwald ripening test taken from PFHub website
 * @version 0.1
 * @date 2026-05-18
 *
 * Copyright CEA (c) 2028
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "./OstwaldCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

struct TestParameters {
  // PDEs
  bool control_restart = false;

  // Time
  double control_time_step = 1.e-2;
  double control_initial_time = 0.;
  double control_final_time = 1.e-1;  // 2.e4;
  int control_iteration_restart = 0;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.control_restart, "-r", "--restart", "-nr", "--no-restart",
                 "Enables restart calculation.");
  args.AddOption(&p.control_initial_time, "-i", "--initial_time",
                 "Initial time of the simulation.");
  args.AddOption(&p.control_final_time, "-t", "--final_time", "Final time of the simulation.");
  args.AddOption(&p.control_time_step, "-dt", "--time_step",
                 "Constant time-step of the simulation.");
  args.AddOption(&p.control_iteration_restart, "-s", "--iteration-restart",
                 "Iteration to restart the calculation.");

  args.Parse();

  if (!args.Good()) {
    if (mfem::Mpi::WorldRank() == 0) {
      args.PrintUsage(mfem::out);
      std::exit(EXIT_FAILURE);
    }
  }
  if (mfem::Mpi::WorldRank() == 0) args.PrintOptions(mfem::out);
}

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

  // ################ //
  // ################ //
  //   Read options   //
  // ################ //
  // ################ //
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);
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
  const std::string& main_folder_path = "Saves";
  const std::string& gf_folder_path = "Restart";
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string mesh_type =
      "InlineSquareWithQuadrangles";  // type of mesh // "InlineSquareWithTriangles"
  const int order_fe = 1;             // finite element order

  const int refinement_level = 0;  // number of levels of uniform refinement
  const int nx = 100;
  const int ny = 100;
  const double lx = 200;
  const double ly = 200;
  mfem::Vector x_translation({lx, 0.0});
  mfem::Vector y_translation({0.0, ly});
  std::vector<mfem::Vector> translations = {x_translation, y_translation};
  const std::tuple<int, int, double, double>& tuple_of_dimensions =
      std::make_tuple(nx, ny, lx, ly);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions, translations);
  // ##############################
  //     Boundary conditions     //
  // ##############################

  auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                     Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
  // CahnHilliard
  auto ch_bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);
  // allenCahn
  auto ac_bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     Coefficients    //
  // ####################

  // CahnHilliard
  const double ch_mob(5.);
  const double ch_lambda = 3.;
  Coefficient ch_grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradC());
  Coefficient ch_energy(Glossary::FreeEnergy, Scheme::Implicit, Gc());
  Coefficient ch_capillary(Glossary::Capillary, ch_lambda);
  Coefficient ch_mobility(Glossary::Mobility, ch_mob);
  Coefficients ch_coef(ch_energy, ch_capillary, ch_mobility, ch_grad_energy);

  // AllenCahn
  const double ac_mob(5.);
  const double ac_lambda = 3.;
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
        double co = 0.5;
        double epsilon = 0.05;
        double xx = x[0];
        double yy = x[1];

        double sol =
            co + epsilon * (std::cos(0.105 * xx) * std::cos(0.11 * yy) +
                            (std::cos(0.13 * xx) * std::cos(0.087 * yy)) *
                                (std::cos(0.13 * xx) * std::cos(0.087 * yy)) +
                            (std::cos(0.025 * xx - 0.15 * yy) * std::cos(0.07 * xx - 0.02 * yy)));

        return sol;
      });

  auto ch_initial_condition_1 = AnalyticalFunctions<DIM>(ch_initial_condition);
  double mu_initial_condition = 0.0;

  auto ch_v1 = VAR(&spatial, ch_bcs, "c", Glossary::PhaseField, 2, ch_initial_condition);
  if (p.control_restart) {
    std::string file_restart = "c_" + std::to_string(p.control_iteration_restart) + ".gf";

    ch_v1 = VAR(&spatial, ch_bcs, "c", Glossary::PhaseField, 2,
                std::make_tuple(file_restart, gf_folder_path));
  }
  ch_v1.set_additional_information("c");

  auto ch_v2 = VAR(&spatial, ch_bcs, "mu", Glossary::ChemicalPotential, 2, 0.);
  if (p.control_restart) {
    std::string file_restart = "mu_" + std::to_string(p.control_iteration_restart) + ".gf";

    ch_v2 = VAR(&spatial, ch_bcs, "mu", Glossary::PhaseField, 2,
                std::make_tuple(file_restart, gf_folder_path));
  }
  ch_v2.set_additional_information("mu");

  auto ch_vars = VARS(ch_v1, ch_v2);

  // AllenCahn
  auto ac_initial_condition_1 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double i = 1.0;
        double epsilon_eta = 0.1;
        double psi = 1.5;
        double xx = x[0];
        double yy = x[1];

        double sol =
            epsilon_eta *
            std::pow(
                (std::cos(xx * i * 0.01 - 4) * std::cos(yy * (0.007 + 0.01 * i)) +
                 std::cos(xx * (0.11 + 0.01 * i)) * std::cos(yy * (0.11 + 0.01 * i)) +
                 psi * std::pow((std::cos(xx * (0.001 * i + 0.046) + yy * (0.0405 + 0.001 * i)) *
                                 std::cos(xx * (0.031 + 0.001 * i) - yy * (0.004 + 0.001 * i))),
                                2.0)),
                2.0);

        return sol;
      });
  auto ac_initial_condition_2 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double i = 2.0;
        double epsilon_eta = 0.1;
        double psi = 1.5;
        double xx = x[0];
        double yy = x[1];

        double sol =
            epsilon_eta *
            std::pow(
                (std::cos(xx * i * 0.01 - 4) * std::cos(yy * (0.007 + 0.01 * i)) +
                 std::cos(xx * (0.11 + 0.01 * i)) * std::cos(yy * (0.11 + 0.01 * i)) +
                 psi * std::pow((std::cos(xx * (0.001 * i + 0.046) + yy * (0.0405 + 0.001 * i)) *
                                 std::cos(xx * (0.031 + 0.001 * i) - yy * (0.004 + 0.001 * i))),
                                2.0)),
                2.0);

        return sol;
      });
  auto ac_initial_condition_3 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double i = 3.0;
        double epsilon_eta = 0.1;
        double psi = 1.5;
        double xx = x[0];
        double yy = x[1];

        double sol =
            epsilon_eta *
            std::pow(
                (std::cos(xx * i * 0.01 - 4) * std::cos(yy * (0.007 + 0.01 * i)) +
                 std::cos(xx * (0.11 + 0.01 * i)) * std::cos(yy * (0.11 + 0.01 * i)) +
                 psi * std::pow((std::cos(xx * (0.001 * i + 0.046) + yy * (0.0405 + 0.001 * i)) *
                                 std::cos(xx * (0.031 + 0.001 * i) - yy * (0.004 + 0.001 * i))),
                                2.0)),
                2.0);

        return sol;
      });
  auto ac_initial_condition_4 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double i = 4.0;
        double epsilon_eta = 0.1;
        double psi = 1.5;
        double xx = x[0];
        double yy = x[1];

        double sol =
            epsilon_eta *
            std::pow(
                (std::cos(xx * i * 0.01 - 4) * std::cos(yy * (0.007 + 0.01 * i)) +
                 std::cos(xx * (0.11 + 0.01 * i)) * std::cos(yy * (0.11 + 0.01 * i)) +
                 psi * std::pow((std::cos(xx * (0.001 * i + 0.046) + yy * (0.0405 + 0.001 * i)) *
                                 std::cos(xx * (0.031 + 0.001 * i) - yy * (0.004 + 0.001 * i))),
                                2.0)),
                2.0);

        return sol;
      });

  auto ac_v1 = VAR(&spatial, ac_bcs, "eta_1", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_1));

  if (p.control_restart) {
    std::string file_restart = "eta_1_" + std::to_string(p.control_iteration_restart) + ".gf";

    ac_v1 = VAR(&spatial, ac_bcs, "eta_1", Glossary::PhaseField, 2,
                std::make_tuple(file_restart, gf_folder_path));
  }
  ac_v1.set_additional_information("eta");

  auto ac_v2 = VAR(&spatial, ac_bcs, "eta_2", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_2));
  if (p.control_restart) {
    std::string file_restart = "eta_2_" + std::to_string(p.control_iteration_restart) + ".gf";
    ac_v2 = VAR(&spatial, ac_bcs, "eta_2", Glossary::PhaseField, 2,
                std::make_tuple(file_restart, gf_folder_path));
  }
  ac_v2.set_additional_information("eta");

  auto ac_v3 = VAR(&spatial, ac_bcs, "eta_3", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_3));
  if (p.control_restart) {
    std::string file_restart = "eta_3_" + std::to_string(p.control_iteration_restart) + ".gf";

    ac_v3 = VAR(&spatial, ac_bcs, "eta_3", Glossary::PhaseField, 2,
                std::make_tuple(file_restart, gf_folder_path));
  }
  ac_v3.set_additional_information("eta");

  auto ac_v4 = VAR(&spatial, ac_bcs, "eta_4", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_4));
  if (p.control_restart) {
    std::string file_restart = "eta_4_" + std::to_string(p.control_iteration_restart) + ".gf";

    ac_v4 = VAR(&spatial, ac_bcs, "eta_4", Glossary::PhaseField, 2,
                std::make_tuple(file_restart, gf_folder_path));
  }
  ac_v4.set_additional_information("eta");

  auto ac_vars = VARS(ac_v1, ac_v2, ac_v3, ac_v4);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const int level_of_detail = 1;
  const auto& frequency = 2000;

  std::vector<int> iterations_list_save_gf = {5,       100000,  200000,  400000,  600000, 800000,
                                              1000000, 1200000, 1400000, 1600000, 1800000};
  if (p.control_restart) {
    int it = p.control_iteration_restart;
    iterations_list_save_gf.erase(
        std::remove_if(iterations_list_save_gf.begin(), iterations_list_save_gf.end(),
                       [it](int v) { return v <= it; }),
        iterations_list_save_gf.end());
  }

  std::string calculation_path = "CahnHilliard";
  const double threshold = 10.;
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {{"phi", {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  auto ch_p_pst = Parameters(
      Parameter("main_folder_path", main_folder_path), Parameter("gf_folder_path", gf_folder_path),
      Parameter("calculation_path", "CahnHilliard"),
      Parameter("iterations_list_save_gf", iterations_list_save_gf),
      Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
      Parameter("integral_to_compute", map_threshold_integral),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));

  auto ac_p_pst = Parameters(
      Parameter("main_folder_path", main_folder_path), Parameter("gf_folder_path", gf_folder_path),
      Parameter("calculation_path", "AllenCahn"),
      Parameter("iterations_list_save_gf", iterations_list_save_gf),
      Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
      Parameter("enable_compute_energies", false),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  // ####################
  //     operators     //
  // ####################

  // CahnHilliard
  std::vector<SPA*> spatials{&spatial, &spatial};
  OPE ch_oper(spatials, {"CahnHilliard"}, TimeScheme::EulerImplicit, "SplitTimeDerivative");
  ch_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));

  auto ch_pst = PST(&spatial, ch_p_pst);
  PB ch_pb("CahnHilliard", ch_oper, ch_vars, {ch_coef, ch_coef}, ch_pst, ac_vars);

  std::vector<SPA*> ac_spatials{&spatial, &spatial, &spatial, &spatial};
  OPE ac_oper(ac_spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  ac_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "LBFGS solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_ILU;
  ac_oper.overload_solver(solver);
  ac_oper.overload_preconditioner(precond);
  auto ac_pst = PST(&spatial, ac_p_pst);
  PB ac_pb("AllenCahn", ac_oper, ac_vars, {ac_coef, ac_coef, ac_coef, ac_coef}, ac_pst, ch_vars);

  // Coupling 1
  auto cc = Coupling("CahnHilliard/AllenCahn Coupling", ac_pb, ch_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = p.control_initial_time;
  const double t_final = p.control_final_time;
  const double dt = p.control_time_step;
  auto time_params_ = Parameters(Parameter("initial_time", t_initial),
                                 Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time_params = time_params_;
  if (p.control_restart) {
    auto restart_parameters =
        Parameters(Parameter("initial_iteration", p.control_iteration_restart));
    time_params = time_params_ + restart_parameters;
  }
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
