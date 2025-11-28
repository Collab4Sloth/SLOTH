/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief Test for external loading f = (r,z,t)
 * @version 0.1
 * @date 2025-11-28
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

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
  //
  //---------------------------------------
  /////////////////////////
  const int DIM = 3;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  int refinement_level = 0;
  double L = 1.;
  int NN = 32;
  int order = 1;
  int level_of_storage = 2;
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------
  SPA spatial("InlineSquareWithTetraedres", order, refinement_level,
              std::make_tuple(NN, NN, 8, L, L, L));
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("rear", 0, "Neumann", 0.),    Boundary("lower", 1, "Neumann", 0.),
                     Boundary("right", 2, "Dirichlet", 1.), Boundary("upper", 3, "Neumann", 0.),
                     Boundary("left", 4, "Dirichlet", 0.),  Boundary("front", 5, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  auto fun1 = std::function<double(const std::vector<double>&, double)>(
      [](const std::vector<double>& vcoord, double time) {
        double m = 9.;
        double res = 1.;
        for (std::size_t i = 0; i < vcoord.size(); i++) {
          res *= std::pow(std::sin(M_PI * m * vcoord[i]), 2.);
        }
        res *= time;
        return res;
      });
  auto fun2 = std::function<double(const std::vector<double>&, double)>(
      [](const std::vector<double>& vcoord, double time) {
        double m = 3.;
        double res = 1.;
        for (std::size_t i = 0; i < vcoord.size(); i++) {
          res *= std::pow(std::cos(M_PI * m * vcoord[i]), 2.);
        }
        res *= time;
        return res;
      });
  auto fun3 = std::function<double(const std::vector<double>&, double)>(
      [](const std::vector<double>& vcoord, double time) {
        double res = 1.;
        for (std::size_t i = 0; i < vcoord.size(); i++) {
          res *= vcoord[i];
        }
        return res;
      });
  std::vector<std::function<double(const std::vector<double>&, double)>> funcs;
  funcs.push_back(fun1);
  funcs.push_back(fun2);
  funcs.push_back(fun3);
  auto param = Parameter("LoadingFunctions", funcs);
  auto params = Parameters(param);
  // ####################
  //     variables     //
  // ####################
  auto user_func_solution =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& v, double time) {
        const auto func = 0.;
        return func;
      });

  auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto var = VAR(&spatial, bcs, "x1", 2, initial_condition);
  auto var2 = VAR(&spatial, bcs, "x2", 2, initial_condition);
  auto var3 = VAR(&spatial, bcs, "x3", 2, initial_condition);
  auto vars = VARS(var, var2, var3);

  // COORDS
  auto xcoord = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& vcoord, double time) { return vcoord[0]; });
  auto ycoord = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& vcoord, double time) { return vcoord[1]; });
  auto zcoord = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& vcoord, double time) { return vcoord[2]; });
  auto XC = VAR(&spatial, bcs, "XCOORD", level_of_storage, AnalyticalFunctions<DIM>(xcoord));
  XC.set_additional_information("XCOORD");
  auto YC = VAR(&spatial, bcs, "YCOORD", level_of_storage, AnalyticalFunctions<DIM>(ycoord));
  YC.set_additional_information("YCOORD");
  auto ZC = VAR(&spatial, bcs, "ZCOORD", level_of_storage, AnalyticalFunctions<DIM>(zcoord));
  ZC.set_additional_information("ZCOORD");
  auto coord = VARS(XC, YC, ZC);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  std::string calculation_path = "Problem1";
  auto p_pst1 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst = PST(&spatial, p_pst1);
  // ####################
  //     operators     //
  // ####################

  // Problem 1:
  std::vector<SPA*> spatials{&spatial};
  ExternalLoadingProblem<FittedLoading, CoordinateSystem::Axisymmetric, VARS, PST> problem1(
      "Fitted external loading example", params, vars, pst, coord);

  // Coupling 1
  auto cc = Coupling("Test", problem1);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 5.;
  const auto& dt = 1.;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

  // time.get_tree();
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