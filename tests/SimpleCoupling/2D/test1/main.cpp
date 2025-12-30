/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D coupling problem Thermal-AC
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

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Debug);
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
  using AC_OPE = PhaseFieldOperator<FECollection, DIM>;
  using AC_PB = Problem<AC_OPE, VARS, PST>;

  // Heat
  using HEAT_OPE = DiffusionOperator<FECollection, DIM>;
  using HEAT_PB = Problem<HEAT_OPE, VARS, PST>;

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
  auto ac_vars = VARS(VAR(&spatial, bcs, "phi", Glossary::PhaseField, 2, ac_ic));

  // Heat
  auto temp = VAR(&spatial, Tbcs, "T", Glossary::Temperature, 2, 1073.15);
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

  // ####################
  //     coefficients  //
  // ####################

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);

  Coefficient density(Glossary::Concentration, rho);
  Coefficient heat_capacity(Glossary::Cp, cp);
  Coefficient conductivity(Glossary::Conductivity, cond);
  Coefficient interpolation(Glossary::InterpolationFunction, Scheme::Implicit, H());

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
  auto p_pst2 = Parameters(
      Parameter("main_folder_path", main_folder_path),
      Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
      Parameter("level_of_detail", level_of_detail), Parameter("enable_compute_energies", false));
  auto pst2 = PST(&spatial, p_pst2);

  // ####################
  //     Probelms      //
  // ####################

  // AllenCahn:
  Coefficients ac_coef(double_well, capillary, mobility, interpolation, grad_energy);
  std::vector<SPA*> spatials{&spatial};
  AC_OPE oper(spatials, {"AllenCahn", "MeltingTemperature"}, ac_params, TimeScheme::EulerImplicit,
              "TimeDerivative");

  AC_PB allencahn_pb("AllenCahn", oper, ac_vars, {ac_coef}, pst, heat_vars);

  // Heat:
  Coefficients heat_coef(density, heat_capacity, conductivity);
  std::vector<AnalyticalFunctions<DIM> > src_term;
  src_term.emplace_back(AnalyticalFunctions<DIM>(src_func));
  HEAT_OPE oper_heat(spatials, {"Fourier"}, TimeScheme::EulerImplicit, "HeatTimeDerivative",
                     src_term);
  oper_heat.overload_nl_solver(
      NLSolverType::NEWTON, Parameters(Parameter("description", "Newton solver "),
                                       Parameter("print_level", 1), Parameter("rel_tol", 1.e-11),
                                       Parameter("abs_tol", 1.e-11), Parameter("iter_max", 1000)));

  HEAT_PB heat_pb("Heat", oper_heat, heat_vars, {heat_coef}, pst2);

  // Coupling 1
  auto cc = Coupling("AC-Heat coupling", allencahn_pb, heat_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 20.;
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
