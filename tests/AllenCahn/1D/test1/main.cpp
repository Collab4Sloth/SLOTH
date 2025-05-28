/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 1D AllenCahn problem along a radius
 * @version 0.1
 * @date 2024-05-23
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

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
  /////////////////////////
  const int DIM = 1;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////
  using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using NLFI2 = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Explicit,
                                          ThermodynamicsPotentials::W, Mobility::Constant>;
  using NLFI3 = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::SemiImplicit,
                                          ThermodynamicsPotentials::W, Mobility::Constant>;

  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
  using OPE2 = AllenCahnOperator<FECollection, DIM, NLFI2>;
  using OPE3 = AllenCahnOperator<FECollection, DIM, NLFI3>;

  using PB = Problem<OPE, VARS, PST>;
  using PB2 = Problem<OPE2, VARS, PST>;
  using PB3 = Problem<OPE3, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string& mesh_type = "InlineLineWithSegments";  // type of mesh
  const int order_fe = 1;                                   // finite element order
  const int refinement_level = 0;  // number of levels of uniform refinement
  const std::tuple<int, double>& tuple_of_dimensions =
      std::make_tuple(30, 1.e-3);  // Number of elements and maximum length

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  //  Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                           Parameter("lambda", lambda), Parameter("omega", omega));
  // ####################
  //     variables     //
  // ####################
  const auto& center_x = 0.;
  const auto& a_x = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = 5.e-4;

  auto user_func = std::function<double(const mfem::Vector&, double)>(
      [center_x, a_x, radius, thickness](const mfem::Vector& x, double time) {
        const auto xx = a_x * (x[0] - center_x);
        const auto r = xx;
        const auto func = 0.5 + 0.5 * std::tanh(2. * (r - radius) / thickness);
        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func);
  auto analytical_solution = AnalyticalFunctions<DIM>(
      AnalyticalFunctionsType::from("HyperbolicTangent"), center_x, a_x, epsilon, radius);
  auto vars = VARS(VAR(&spatial, bcs, "phi1", 2, initial_condition, analytical_solution));
  auto vars2 = VARS(VAR(&spatial, bcs, "phi2", 2, initial_condition, analytical_solution));
  auto vars3 = VARS(VAR(&spatial, bcs, "phi3", 2, initial_condition, analytical_solution));
  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  // ####################
  //     operators     //
  // ####################
  std::string calculation_path = "Problem1";
  auto p_pst1 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));

  // Problem 1:
  const auto crit_cvg_1 = 1.e-12;
  OPE oper(&spatial, params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
  auto pst = PST(&spatial, p_pst1);
  PB problem1(oper, vars, pst, convergence);

  // Problem 2:
  const auto crit_cvg_2 = 1.e-12;
  calculation_path = "Problem2";
  auto p_pst2 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  OPE2 oper2(&spatial, params, TimeScheme::EulerExplicit);
  oper2.overload_mobility(Parameters(Parameter("mob", mob)));
  PhysicalConvergence convergence2(ConvergenceType::RELATIVE_MAX, crit_cvg_2);
  auto pst2 = PST(&spatial, p_pst2);
  PB2 problem2(oper2, vars2, pst2, convergence2);

  // Problem 3:
  calculation_path = "Problem3";
  auto p_pst3 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  const auto crit_cvg_3 = 1.e-12;
  OPE3 oper3(&spatial, params, TimeScheme::RungeKutta4);
  oper3.overload_mobility(Parameters(Parameter("mob", mob)));
  PhysicalConvergence convergence3(ConvergenceType::RELATIVE_MAX, crit_cvg_3);
  auto pst3 = PST(&spatial, p_pst3);
  PB3 problem3(oper3, vars3, pst3, convergence3);

  // Coupling 1
  auto cc = Coupling("coupling 1 ", problem1, problem2, problem3);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 50.0;
  const auto& dt = 0.01;
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
