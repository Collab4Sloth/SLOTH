/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief  Homogeneous nucleation from PFHub website
 * @version 0.1
 * @date 2026-05-07
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
  using BCS = Test<DIM>::BCS;
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

  std::vector<std::tuple<std::string, double>> vect_radius{
      std::make_tuple("099", 0.99),
      std::make_tuple("1", 1),
      std::make_tuple("101", 101),
  };
  for (const auto& tup_radius : vect_radius) {
    const int NN = 256;
    const int order_fe = 1;          // finite element order
    const int refinement_level = 0;  // number of levels of uniform refinement
    const double L = 100.;
    mfem::Vector x_translation({L, 0.0});
    mfem::Vector y_translation({0.0, L});
    std::vector<mfem::Vector> translations = {x_translation, y_translation};
    const std::tuple<int, int, double, double>& tuple_of_dimensions =
        std::make_tuple(NN, NN, L, L);  // Number of elements and maximum length in each direction

    SPA spatial("InlineSquareWithQuadrangles", order_fe, refinement_level, tuple_of_dimensions,
                translations);
    // ##############################
    //     Boundary conditions     //
    // ##############################
    auto boundaries_phi = {Boundary("lower", 0, "Periodic"), Boundary("right", 1, "Periodic"),
                           Boundary("upper", 2, "Periodic"), Boundary("left", 3, "Periodic")};
    auto bcs_phi = BCS(&spatial, boundaries_phi);

    // ###########################################
    // ###########################################
    //            Physical models               //
    // ###########################################
    // ###########################################
    // ####################
    //     coefficients  //
    // ####################
    const double rstar = 5.0;
    const double rro = rstar * std::get<1>(tup_radius);
    const double DeltaF = -std::sqrt(2.0) / (6.0 * rstar);
    //  Interface thickness
    const double epsilon(1.);
    // Interfacial energy
    const double sigma(1.);
    // Two-phase mobility
    const double mob(1.);
    const double lambda = 1.;
    const double omega = 1.;

    Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
    Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
    Coefficient interpolation(Glossary::InterpolationFunction, Scheme::Implicit, H());
    Coefficient capillary(Glossary::Capillary, lambda);
    Coefficient mobility(Glossary::Mobility, mob);
    // ####################
    //     variables     //
    // ####################

    auto initial_condition = std::function<double(const mfem::Vector&, double)>(
        [rro](const mfem::Vector& x, double time) {
          double xx = x[0] - 50.0;
          double yy = x[1] - 50.0;

          double rr = std::sqrt(xx * xx + yy * yy);

          double sol = 0.5 * (1.0 - std::tanh((rr - rro) / std::sqrt(2.0)));
          return sol;
        });

    const std::string& var_name_1 = "phi";
    auto v1 = VAR(&spatial, bcs_phi, var_name_1, Glossary::PhaseField, 2,
                  AnalyticalFunctions<DIM>(initial_condition));
    auto vars = VARS(v1);

    // ###########################################
    // ###########################################
    //      Post-processing                     //
    // ###########################################
    // ###########################################

    const std::string& main_folder_path = "Saves_" + std::get<0>(tup_radius);
    const int level_of_detail = 1;
    const int frequency = 200;
    std::string calculation_path = "Problem1";
    bool enable_save_specialized_at_iter = true;
    std::map<std::string, std::tuple<double, double>> map_threshold_integral = {
        {var_name_1, {-0.1, 1.1}}};
    auto p_pst =
        Parameters(Parameter("main_folder_path", main_folder_path),
                   Parameter("calculation_path", calculation_path),
                   Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
                   Parameter("integral_to_compute", map_threshold_integral),
                   Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));

    // ####################
    //     operators     //
    // ####################
    auto ac_param = Parameters(Parameter("melting_factor", DeltaF));
    // Problem 1:
    std::vector<SPA*> spatials{&spatial};
    OPE oper(spatials, {"AllenCahn", "MeltingConstant"}, ac_param, TimeScheme::EulerImplicit,
             "TimeDerivative");
    oper.overload_nl_solver(NLSolverType::NEWTON,
                            Parameters(Parameter("description", "Newton solver "),
                                       Parameter("print_level", 1), Parameter("rel_tol", 1.e-10),
                                       Parameter("abs_tol", 1.e-14), Parameter("iter_max", 1000)));
    const auto& solver = HypreSolverType::HYPRE_GMRES;
    const auto& precond = HyprePreconditionerType::HYPRE_ILU;
    oper.overload_solver(solver);
    oper.overload_preconditioner(precond);
    auto pst = PST(&spatial, p_pst);
    Coefficients coef_pb1(double_well, capillary, mobility, grad_energy, interpolation);

    PB problem1(oper, vars, {coef_pb1}, pst);

    // Coupling 1
    auto cc = Coupling("AllenCahn Coupling", problem1);

    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const double t_initial = 0.0;
    const double t_final = 0.1;  // 100.;
    const double dt = 1.e-2;
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
