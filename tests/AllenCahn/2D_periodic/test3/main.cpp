/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief Allen-Cahn problem solved in a square using MMS
 * @version 0.1
 * @date 2024-07-30
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include <map>
#include <memory>
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
  int size = mfem::Mpi::WorldSize();
  int rank = mfem::Mpi::WorldRank();
  mfem::Hypre::Init();
  //
  // ---------------------------------------
  const auto DIM = 2;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI, SteadyPhaseFieldOperatorBase>;

  using PB = Problem<OPE, VAR, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  auto L = M_PI;
  std::vector<int> vect_order{1, 2};
  std::vector<int> vect_NN{8, 16, 32};

  for (const auto& order : vect_order) {
    for (const auto& NN : vect_NN) {
      //---------------------------------------
      // Profiling start
      Profiling::getInstance().enable();
      //---------------------------------------
      mfem::Vector x_translation({L, 0.0});
      mfem::Vector y_translation({0.0, L});
      std::vector<mfem::Vector> translations = {x_translation, y_translation};
      SpatialDiscretization<FECollection, DIM> spatial("InlineSquareWithQuadrangles", order,
                                                       refinement_level,
                                                       std::make_tuple(NN, NN, L, L), translations);
      // ##############################
      //     Boundary conditions     //
      // ##############################
      auto boundaries = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                         Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
      auto bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

      // ###########################################
      // ###########################################
      //            Physical models               //
      // ###########################################
      // ###########################################
      // ####################
      //     parameters    //
      // ####################
      //  Interface thickness
      const auto& epsilon(1.);
      // Interfacial energy
      const auto& sigma(1.);
      // Two-phase mobility
      const auto& mob(1.);
      const auto& lambda = 1.;
      const auto& omega = 1.;
      auto params = Parameters(Parameter("epsilon", epsilon), Parameter("mobility", mob),
                               Parameter("sigma", sigma), Parameter("lambda", lambda),
                               Parameter("omega", omega));
      // ####################,
      //     variables     //
      // ####################

      auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double x = v[0];
            const double y = v[1];
            const auto func = std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2);
            return func;
          });
      auto user_func_source_term = std::function<double(const mfem::Vector&, double)>(
          [](const mfem::Vector& v, double time) {
            const double x = v[0];
            const double y = v[1];
            const double H = 1.;
            const double a = 1.;
            const auto func =
                2 * H * (-2 * std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) + 1) *
                    (-std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) + 1) *
                    std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) -
                a * (18 * std::pow(std::sin(2 * x), 2) * std::pow(std::sin(3 * y), 2) -
                     26 * std::pow(std::sin(2 * x), 2) * std::pow(std::cos(3 * y), 2) +
                     8 * std::pow(std::cos(2 * x), 2) * std::pow(std::cos(3 * y), 2));
            return func;
          });

      auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
      auto analytical_solution = AnalyticalFunctions<DIM>(user_func_solution);

      auto vars = VAR(Variable<FECollection, DIM>(&spatial, bcs, "phi", 2, initial_condition,
                                                  analytical_solution));

      // ###########################################
      // ###########################################
      //      Post-processing                     //
      // ###########################################
      // ###########################################
      const std::string& main_folder_path =
          "Saves_order_" + std::to_string(order) + "_Nx" + std::to_string(NN);
      const auto& level_of_detail = 1;
      const auto& frequency = 1;
      // ####################
      //     operators     //
      // ####################

      // Problem 1:
      const auto crit_cvg_1 = 1.e-12;
      auto src_term = AnalyticalFunctions<DIM>(user_func_source_term);
      OPE oper(&spatial, params, src_term);
      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);
      auto pst = PST(main_folder_path, "Problem1", &spatial, frequency, level_of_detail);
      PB problem1("PSteady AllenCahn", oper, vars, pst, convergence);

      // Coupling 1
      auto cc = Coupling("Default Coupling", problem1);

      // ###########################################
      // ###########################################
      //            Time-integration              //
      // ###########################################
      // ###########################################
      const auto& t_initial = 0.0;
      const auto& t_final = 0.25;
      const auto& dt = 0.25;
      auto time_params = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
      auto time = TimeDiscretization(time_params, cc);

      // time.get_tree();
      time.solve();
      //---------------------------------------
      // Profiling stop
      //---------------------------------------
      Profiling::getInstance().print();
    }
  }

  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
