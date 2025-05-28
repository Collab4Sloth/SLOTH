/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D coupling problem
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
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  /////////////////////////

  // ALLEN-CAHN
  using NLFI = AllenCahnTemperatureMeltingNLFormIntegrator<
      VARS, ThermodynamicsPotentialDiscretization::Implicit, ThermodynamicsPotentials::W,
      Mobility::Constant, ThermodynamicsPotentials::H>;
  using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
  using PB = Problem<OPE, VARS, PST>;

  // Heat
  using NLFI2 =
      HeatNLFormIntegrator<VARS, CoefficientDiscretization::Explicit, Conductivity::Constant>;
  using OPE2 = HeatOperator<FECollection, DIM, NLFI2, Density::Constant, HeatCapacity::Constant>;
  using PB2 = Problem<OPE2, VARS, PST>;

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  SPA spatial("GMSH", 1, refinement_level, "camembert2D.msh", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  // ALLEN-CAHN
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Neumann", 0.),
                     Boundary("upper", 1, "Neumann", 0.)};
  auto Tboundaries = {Boundary("lower", 0, "Neumann", 0.),
                      Boundary("external", 2, "Dirichlet", 1073.15),
                      Boundary("upper", 1, "Neumann", 0.)};
  auto bcs = BCS(&spatial, boundaries);
  auto Tbcs = BCS(&spatial, Tboundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     parameters    //
  // ####################
  // ALLEN-CAHN
  //    Melting factor
  const auto& dH(7.e3);
  const auto& tf(1500.);
  // Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-5);
  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto ac_params =
      Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                 Parameter("lambda", lambda), Parameter("omega", omega),
                 Parameter("melting_temperature", tf), Parameter("melting_enthalpy", dH));
  // Heat
  const auto& rho(10.e4);
  const auto& cp(400.);
  const auto& cond(2.);

  // ############################
  //     variables IC + SRC    //
  // ###########################
  // ALLEN-CAHN
  const auto& pellet_radius = 0.00465;
  const auto& center_x = 0.;
  const auto& center_y = 0.;
  const auto& a_x = 1.;
  const auto& a_y = 1.;
  const auto& thickness = 5.e-5;
  const auto& radius = 1.e-2 * pellet_radius;

  auto ac_ic = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x,
                                        center_y, a_x, a_y, thickness, radius);
  auto ac_vars = VARS(VAR(&spatial, bcs, "phi", 2, ac_ic));

  // Heat
  auto temp = VAR(&spatial, Tbcs, "T", 2, 1073.15);
  temp.set_additional_information("T");
  auto heat_vars = VARS(temp);
  auto pl = 15.e4;
  auto src_func = std::function<double(const mfem::Vector&, double)>(
      [pl, pellet_radius](const mfem::Vector& vcoord, double time) {
        const double radius = std::sqrt(vcoord[0] * vcoord[0] + vcoord[1] * vcoord[1]);
        // auto chi = 90.;  // inverse neutron diffusion length (0.9cm−1 ->90m-1).
        // auto chia = chi * pellet_radius;
        // auto I1_chia = std::cyl_bessel_i(1, chia);
        // auto chir =
        //     chi * (pellet_radius - radius);  //  radius;  //  (Radius[_nnodes - 1] - Radius[i]);1
        // auto I0_chir = std::cyl_bessel_i(0, chir);
        const auto bess = 1.;  // chia * I0_chir / (2. * I1_chia);

        const auto func = pl * bess / (M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
  // Allen-Cahn
  std::string calculation_path = "AllenCahn";
  auto p_pst1 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst = PST(&spatial, p_pst1);
  // Heat
  calculation_path = "Heat";
  auto p_pst2 =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto pst2 = PST(&spatial, p_pst2);

  // ####################
  //     Probelms      //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  // AllenCahn:
  OPE oper(&spatial, ac_params, TimeScheme::EulerImplicit);
  oper.overload_mobility(Parameters(Parameter("mob", mob)));

  PB allencahn_pb("AllenCahn", oper, ac_vars, pst, convergence, heat_vars);

  // Heat:
  auto source_term = AnalyticalFunctions<DIM>(src_func);
  OPE2 oper2(&spatial, TimeScheme::EulerImplicit, source_term);
  oper2.overload_density(Parameters(Parameter("rho", rho)));
  oper2.overload_heat_capacity(Parameters(Parameter("cp", cp)));
  oper2.overload_conductivity(Parameters(Parameter("lambda", cond)));

  oper2.overload_nl_solver(NLSolverType::NEWTON,
                           Parameters(Parameter("description", "Newton solver "),
                                      Parameter("print_level", 1), Parameter("abs_tol", 1.e-9)));

  PB2 Heat_pb("Heat", oper2, heat_vars, pst2, convergence);

  // MPI_Problem<VAR, PST> mpi(vars3, pst3, convergence);

  // Coupling 1
  auto cc = Coupling("AC-Heat coupling", allencahn_pb, Heat_pb);
  // auto cc = Coupling("AC-Heat coupling", Heat_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 0.5;  // 100.;
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
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return 0;
}
