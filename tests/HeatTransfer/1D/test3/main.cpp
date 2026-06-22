/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D axisymmetric heat transfer problem
 * @version 0.1
 * @date 2024-09-3
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Debug);
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
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  //
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
  const std::string& mesh_type = "InlineLineWithSegments";  // type of mesh
  const int refinement_level = 0;  // number of levels of uniform refinement

  const double pellet_radius = 0.00465;
  std::vector<int> vect_NN{30, 60, 90, 120};
  for (const auto NN : vect_NN) {
    const int order = 1;  // finite element order
    const int n_elements = NN;

    SPA spatial(mesh_type, order, refinement_level, std::make_tuple(n_elements, pellet_radius));

    // ##############################
    //     Boundary conditions     //
    // ##############################
    const double T_initial = 750.0;
    auto boundaries = {Boundary("left", 0, "Neumann"),
                       Boundary("right", 1, "Dirichlet", T_initial)};
    auto Tbcs = BCS(&spatial, boundaries);
    // ###########################################
    // ###########################################
    //            Physical models               //
    // ###########################################
    // ###########################################
    // ####################
    //     parameters    //
    // ####################
    // Heat
    const auto& rho(1.e3);
    const auto& cp(10.);
    const auto& cond(2.);
    Coefficient density(Glossary::Concentration, rho);
    Coefficient heat_capacity(Glossary::Cp, cp);
    Coefficient conductivity(Glossary::Conductivity, cond);

    // ############################
    //     variables IC + SRC    //
    // ###########################
    // Heat
    auto pl = 26.e3;

    auto src_func = std::function<double(const mfem::Vector&, double)>(
        [pl, pellet_radius]([[maybe_unused]] const mfem::Vector& vcoord,
                            [[maybe_unused]] double time) {
          const double rr = vcoord[0];
          const double func = pl * rr / (M_PI * pellet_radius * pellet_radius);

          return func;
        });

    auto T_analytical = std::function<double(const mfem::Vector&, double)>(
        [pl, pellet_radius, T_initial, cond]([[maybe_unused]] const mfem::Vector& vcoord,
                                             [[maybe_unused]] double time) {
          const double rr = vcoord[0] * vcoord[0];
          const double func = T_initial + pl * (pellet_radius * pellet_radius - rr) /
                                              (4.0 * M_PI * pellet_radius * pellet_radius * cond);

          return func;
        });

    auto heat_vars = VARS(VAR(&spatial, Tbcs, "T", Glossary::Temperature, 2, T_initial,
                              AnalyticalFunctions<DIM>(T_analytical)));
    // ###########################################
    // ###########################################
    //      Post-processing                     //
    // ###########################################
    // ###########################################

    const std::string& main_folder_path =
        "Saves_order_" + std::to_string(order) + "_Nx" + std::to_string(NN);
    const auto& level_of_detail = 1;
    const auto& frequency = 1;
    // Heat
    const std::string& calculation_path = "Problem1";
    auto p_pst = Parameters(
        Parameter("main_folder_path", main_folder_path),
        Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
        Parameter("level_of_detail", level_of_detail), Parameter("enable_compute_energies", false));
    auto pst = PST(&spatial, p_pst);

    // ####################
    //     Problems      //
    // ####################

    // Heat:
    Coefficients coef_pb(density, heat_capacity, conductivity);
    std::vector<AnalyticalFunctions<DIM> > src_term;
    src_term.emplace_back(AnalyticalFunctions<DIM>(src_func));
    std::vector<SPA*> spatials{&spatial};
    OPE oper(spatials, {"Fourier"}, TimeScheme::EulerImplicit, "HeatTimeDerivative", src_term);

    oper.overload_nl_solver(
        NLSolverType::NEWTON,
        Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                   Parameter("rel_tol", 1.e-9), Parameter("abs_tol", 1.e-9)));

    auto T_cvg = PhysicalConvergence(ConvergenceType::ABSOLUTE_MAX, 1.e-16);
    auto CVG = Convergence(T_cvg);

    PB Heat_pb("Heat", oper, heat_vars, {coef_pb}, CVG, pst);
    Heat_pb.setGeometry(Geometry::Axisymmetric);
    // Coupling 1
    auto cc = Coupling("Heat transfer", Heat_pb);
    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const auto& t_initial = 0.0;
    const auto& t_final = 10.;
    const auto& dt = 0.01;
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
