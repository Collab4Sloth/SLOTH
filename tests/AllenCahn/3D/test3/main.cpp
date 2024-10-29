/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Spinodal decomposition in a 3D domain
 * @version 0.1
 * @date 2024-09-3
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  const auto DIM = 3;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::SemiImplicit,
                                         ThermodynamicsPotentials::F, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
  using PB = Problem<OPE, VAR, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################

  const int NN = 32;

  auto refinement_level = 0;

  auto L = 2.;
  SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithHexaedres", 1, refinement_level,
                                                   std::make_tuple(NN, NN, NN, L, L, L));

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("rear", 0, "Neumann", 0.),  Boundary("lower", 1, "Neumann", 0.),
                     Boundary("right", 2, "Neumann", 0.), Boundary("upper", 3, "Neumann", 0.),
                     Boundary("left", 4, "Neumann", 0.),  Boundary("front", 5, "Neumann", 0.)};
  auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  const auto& eps = 0.02;
  const auto& epsilon(eps);
  const auto& mob(1.);
  const auto& lambda = 1.;
  const auto& omega = 1. / (epsilon * epsilon);
  auto params = Parameters(Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################

  auto user_func =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& v, double time) {
        std::random_device rd;   // Seed for the random number generator
        std::mt19937 gen(rd());  // Mersenne Twister random number generator
        std::uniform_real_distribution<> dis(-1.0, 1.0);
        const auto func = 0.01 * dis(gen);
        return func;
      });
  auto initial_condition = AnalyticalFunctions<DIM>(user_func);
  auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, initial_condition));

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 100;
  std::string calculation_path = "Problem1";
  auto p_pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail),
                 Parameter("enable_save_specialized_at_iter", true),
                 Parameter("force_clean_output_dir", true));
  auto pst = PST(&spatial, p_pst);

  // ####################
  //     operator     //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, params, TimeScheme::EulerImplicit);

  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  auto nl_params = Parameters(Parameter("description", "Newton Algorithm"),
                              Parameter("print_level", 1), Parameter("abs_tol", 1.e-8));
  auto s_params = Parameters(Parameter("description", " solver "), Parameter("print_level", 0));
  auto p_params = Parameters(Parameter("description", " preconditionner"), Parameter("type", 0));

  oper.overload_nl_solver(NLSolverType::NEWTON, nl_params);
  oper.overload_solver(HypreSolverType::HYPRE_GMRES, s_params, HyprePreconditionerType::HYPRE_ILU,
                       p_params);
  oper.overload_mass_solver(HypreSolverType::HYPRE_GMRES, s_params,
                            HyprePreconditionerType::HYPRE_SMOOTHER, p_params);

  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  // ####################
  //     Problem       //
  // ####################
  PB problem1(oper, vars, pst, convergence);

  // ####################
  //     Coupling      //
  // ####################
  auto cc = Coupling("Default Coupling", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 4.e-6; //0.012;
  const auto& dt = 2.e-6;
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
