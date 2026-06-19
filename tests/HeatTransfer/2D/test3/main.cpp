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

#include "./Robin.hpp"
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
  const int DIM = 2;
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
  const std::string& mesh_type = "InlineSquareWithQuadrangles";  // type of mesh
  const int order_fe = 1;                                        // finite element order
  const int refinement_level = 1;  // number of levels of uniform refinement
  const double L = 1.0;
  std::vector<int> vect_NN{30, 60, 120, 240};
  for (const auto NN : vect_NN) {
    const int order = 1;  // finite element order
    const int n_elements = NN;

    SPA spatial(mesh_type, order_fe, refinement_level,
                std::make_tuple(n_elements, n_elements, L, L));

    // ##############################
    //     Boundary conditions     //
    // ##############################
    auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                       Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
    auto Xboundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                        Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
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
    const auto& rho(1.);
    const auto& cp(1.);
    const auto& cond(1.);
    Coefficient density(Glossary::Concentration, rho);
    Coefficient heat_capacity(Glossary::Cp, cp);
    Coefficient conductivity(Glossary::Conductivity, cond);

    // ############################
    //     variables IC + SRC    //
    // ###########################
    // Heat
    double qtop = 1.0;
    auto T_analytical = std::function<double(const mfem::Vector&, double)>(
        [qtop]([[maybe_unused]] const mfem::Vector& vcoord, [[maybe_unused]] double time) {
          const double r = vcoord[0];
          const double z = vcoord[1];

          return std::exp(-time) * std::cos(M_PI * r * r) * std::cos(M_PI * z);
        });

    auto heat_vars =
        VARS(VAR(&spatial, Tbcs, "T", Glossary::Temperature, 2,
                 AnalyticalFunctions<DIM>(T_analytical), AnalyticalFunctions<DIM>(T_analytical)));
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
    Coefficients coef_pb(density, heat_capacity, conductivity);  //, robin_b);
    std::vector<SPA*> spatials{&spatial};

    auto src_func = std::function<double(const mfem::Vector&, double)>(
        [rho, cp, cond]([[maybe_unused]] const mfem::Vector& vcoord, [[maybe_unused]] double time) {
          const double r = vcoord[0];
          const double z = vcoord[1];

          const double crr = std::cos(M_PI * r * r);
          const double srr = std::sin(M_PI * r * r);
          const double cz = std::cos(M_PI * z);

          return r * std::exp(-time) *
                 (-rho * cp * crr * cz +
                  cond * (4.0 * M_PI * srr * cz + 4.0 * M_PI * M_PI * r * r * crr * cz +
                          M_PI * M_PI * crr * cz));
        });

    std::vector<AnalyticalFunctions<DIM> > src_term;
    src_term.emplace_back(AnalyticalFunctions<DIM>(src_func));

    OPE oper(spatials, {"Fourier"}, TimeScheme::EulerImplicit, "HeatTimeDerivative", src_term);

    oper.overload_nl_solver(
        NLSolverType::NEWTON,
        Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                   Parameter("rel_tol", 1.e-9), Parameter("abs_tol", 1.e-9)));

    PB Heat_pb("Heat", oper, heat_vars, {coef_pb}, pst);
    Heat_pb.setGeometry(Geometry::Axisymmetric);
    // Coupling 1
    auto cc = Coupling("Heat transfer", Heat_pb);
    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const auto& t_initial = 0.0;
    const auto& t_final = 0.005;
    const auto& dt = 0.001;
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
