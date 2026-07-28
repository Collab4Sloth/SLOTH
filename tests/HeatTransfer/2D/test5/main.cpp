/**
 * @file main.cpp
 * @author mh286406 (marine.harel@cea.fr)
 * @brief 2D heat transfer problem
 * @version 0.1
 * @date 2026-07-17
 *
 * @copyright Copyright (c) 2026
 *
 */
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
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
  const int refinement_level = 0;

  std::vector<std::string> vect_elem{"InlineSquareWithQuadrangles"};
  std::vector<int> vect_order{2, 1};
  std::vector<int> vect_NN{160, 80, 40, 20};
  for (const auto& elem_type : vect_elem) {
    for (const auto order_fe : vect_order) {
      for (const auto NN : vect_NN) {
        const int nx = NN;
        const int ny = NN;
        const double lx = 1.;
        const double ly = 1.;

        const std::tuple<int, int, double, double>& tuple_of_dimensions = std::make_tuple(
            nx, ny, lx, ly);  // Number of elements and maximum length in each direction

        SPA spatial(elem_type, order_fe, refinement_level, tuple_of_dimensions);

        // ##############################
        //     Boundary conditions     //
        // ##############################

        auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet"),
                          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet")};
        auto bcs = BCS(&spatial, boundaries);

        auto Xboundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
          Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
        auto Xbcs = BCS(&spatial, Xboundaries);

        // ###########################################
        // ###########################################
        //            Physical models               //
        // ###########################################
        // ###########################################
        // ####################
        //     parameters    //
        // ####################
        // Heat

        Coefficient density(Glossary::Concentration, 1.);
        Coefficient heat_capacity(Glossary::Cp, 1.);
        Coefficient conductivity(Glossary::Conductivity, 1.);
        Coefficient neumann(Glossary::Neumann, Scheme::Implicit, NeumannCoefficient());
        Coefficient dirichlet_left(Glossary::Dirichlet, Scheme::Implicit, DirichletCoefficient());
  Coefficient dirichlet_right(Glossary::Dirichlet, Scheme::Implicit, DirichletCoefficient(-1));
        neumann.set_bdr_index_coef(std::vector<int>{0,2});
        dirichlet_left.set_bdr_index_coef(std::vector<int>{3});
        dirichlet_right.set_bdr_index_coef(std::vector<int>{1});

        // ####################
        //     variables     //
        // ####################

        auto user_func = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, double time) {
            const auto xx = x[0];
            const auto yy = x[1];
            const auto func = time * std::cos(M_PI*xx) * std::sin(M_PI*yy);
            return func;
          });
        auto T_analytical = AnalyticalFunctions<DIM>(user_func);

  auto heat_vars = VARS(VAR(&spatial, bcs, "T", Glossary::Temperature, 2, 0, T_analytical));

        // Coord
        auto xcoord = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& vcoord, double time) { return vcoord[0]; });
        auto ycoord = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& vcoord, double time) { return vcoord[1]; });
        auto XC = VAR(&spatial, Xbcs, "XCOORD", Glossary::Coordinate, 2,
                    AnalyticalFunctions<DIM>(xcoord));
        XC.set_additional_information("XCOORD");
        auto YC = VAR(&spatial, Xbcs, "YCOORD", Glossary::Coordinate, 2,
                    AnalyticalFunctions<DIM>(ycoord));
        YC.set_additional_information("YCOORD");
        auto coord = VARS(XC, YC);

        // ###########################################
        // ###########################################
        //      Post-processing                     //
        // ###########################################
        // ###########################################
        const std::string& main_folder_path =
            "Saves_order_" + std::to_string(order_fe) + "_Nx" + std::to_string(NN);
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

        auto user_func_source_term = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& x, [[maybe_unused]] double time) {
            const auto xx = x[0];
            const auto yy = x[1];
            const auto func = (1+2*M_PI*M_PI*time) * std::cos(M_PI*xx) * std::sin(M_PI*yy);
            return func;
          });

        std::vector<AnalyticalFunctions<DIM> > src_term;
        src_term.emplace_back(AnalyticalFunctions<DIM>(user_func_source_term));
        
        // Heat
        Coefficients coef_pb(density, heat_capacity, conductivity, neumann, dirichlet_left, dirichlet_right);
        std::vector<SPA*> spatials{&spatial};
  OPE oper(spatials, {"Fourier"}, TimeScheme::EulerImplicit, "HeatTimeDerivative", src_term);

        oper.overload_nl_solver(
            NLSolverType::NEWTON,
            Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                      Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-12)));
        PB Heat_pb("Heat", oper, heat_vars, {coef_pb}, pst, coord);

        // Coupling 1
        auto cc = Coupling("Heat transfer", Heat_pb);
        // ###########################################
        // ###########################################
        //            Time-integration              //
        // ###########################################
        // ###########################################
        const auto& t_initial = 0.0;
        const auto& t_final = 1.;
        const auto& dt = 1.;
        auto time_params = Parameters(Parameter("initial_time", t_initial),
                                      Parameter("final_time", t_final), Parameter("time_step", dt));
        auto time = TimeDiscretization(time_params, cc);

        time.solve();
        //---------------------------------------
        // Profiling stop
        //---------------------------------------
        Profiling::getInstance().print();
      }
    }
  }
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
