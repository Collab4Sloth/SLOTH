/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D Inter-diffusion test for a binary system
 * @version 0.1
 * @date 2025-03-13
 *
 * Copyright (c) 2025
 *
 */

//---------------------------------------
// Headers
//---------------------------------------

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
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  // Common aliases
  //---------------------------------------
  const int DIM = 2;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  //---------------------------------------
  // Meshing & Boundary Conditions
  //---------------------------------------
  const int refinement_level = 0;
  auto length = 1.e-3;
  auto nb_fe = 30;

  SPA spatial("InlineSquareWithQuadrangles", 1, refinement_level,
              std::make_tuple(nb_fe, nb_fe, length, length));

  std::vector<SPA*> spatials{&spatial};
  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto boundaries = {Boundary("lower", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.),
                     Boundary("upper", 2, "Neumann", 0.), Boundary("left", 3, "Neumann", 0.)};

  auto bcs = BCS(&spatial, boundaries);

  //---------------------------------------
  // Multiphysics coupling scheme
  //---------------------------------------
  const int level_of_storage = 2;
  const auto& diffusionCoeff(1.e-8);
  const auto& dt = 0.05;

  auto user_func_init = std::function<double(const mfem::Vector&, double)>(
      [length, diffusionCoeff, dt](const mfem::Vector& x, [[maybe_unused]] double time) {
        const auto xx = x[0];
        auto func =
            0.5 * (1 + std::erf((xx - length / 2) / std::sqrt(4 * diffusionCoeff * 5. * dt)));

        return func;
      });

  auto user_func_analytical = std::function<double(const mfem::Vector&, double)>(
      [length, dt, diffusionCoeff](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        auto func =
            0.5 *
            (1 + std::erf((xx - length / 2) / std::sqrt(4 * diffusionCoeff * (time + 5. * dt))));

        return func;
      });
  auto initial_compo = AnalyticalFunctions<DIM>(user_func_init);
  auto analytical_compo = AnalyticalFunctions<DIM>(user_func_analytical);
  //==========================================
  //======      Fickian diffusion       ======
  //==========================================
  Coefficient D(Glossary::Diffusivity, diffusionCoeff);
  Coefficients coef_fick(D);
  auto diffu_vars = VARS(VAR(&spatial, bcs, "c", Glossary::MoleFraction, level_of_storage,
                             initial_compo, analytical_compo));
  TransientOperator<FECollection, DIM> diffu_oper(spatials, {"Fick"}, TimeScheme::RungeKutta4,
                                                  "TimeDerivative");

  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const auto& stabCoeff(1.e-7);
  auto td_parameters = Parameters(Parameter("last_component", "B"),
                                  Parameter("ScaleCoefficientsByTemperature", false),
                                  Parameter("EnableDiffusionChemicalPotentials", true));
  auto mobA = VAR(&spatial, bcs, "Ma", Glossary::InterDiffusionMobility, 2, diffusionCoeff);
  mobA.set_additional_information("SOLID", "A", "mob");
  auto mobB = VAR(&spatial, bcs, "Mb", Glossary::InterDiffusionMobility, 2, diffusionCoeff);
  mobB.set_additional_information("SOLID", "B", "mob");
  auto mobilities = VARS(mobA, mobB);

  auto inter_mob = VAR(&spatial, bcs, "MA", Glossary::InterDiffusionMobility, 2, 0.);
  inter_mob.set_additional_information("A", "inter_mob");
  auto MA = VARS(inter_mob);

  Coefficient Dstab(Glossary::Diffusivity, stabCoeff);
  Coefficients coef_inter(Dstab);
  TransientOperator<FECollection, DIM> interdiffu_oper(spatials, {"MassFlux"}, td_parameters,
                                                       TimeScheme::RungeKutta4, "TimeDerivative");

  //==========================================
  //======      CALPHAD Analytical      ======
  //==========================================
  //--- Variables
  // Temperature
  auto temp = VAR(&spatial, bcs, "T", Glossary::Temperature, level_of_storage, 1. / Physical::R);
  temp.set_additional_information("K", "T");
  auto heat_vars = VARS(temp);
  // Pressure
  auto pres = VAR(&spatial, bcs, "pressure", Glossary::Pressure, level_of_storage, 1.);
  pres.set_additional_information("Pa", "P");
  auto p_vars = VARS(pres);

  // Initial condition for composition
  auto xa = VAR(&spatial, bcs, "A", Glossary::MoleFraction, 2, initial_compo, analytical_compo);
  xa.set_additional_information("A", "x");
  auto compo_vars = VARS(xa);

  // Chemical potential
  auto mua = VAR(&spatial, bcs, "dmu", Glossary::ChemicalPotential, 2, 0.);
  mua.set_additional_information("A", "dmu");

  auto mu_var = VARS(mua);

  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");
  auto calphad_parameters = Parameters(description_calphad);

  //==========================================
  //==========================================
  //--- Post-Processing
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "MobilitiesA";
  const auto& frequency = 1;
  auto pst_parameters_mob = Parameters(Parameter("main_folder_path", main_folder_path),
                                       Parameter("calculation_path", calculation_path),
                                       Parameter("frequency", frequency));
  auto mob_pst_a = PST(&spatial, pst_parameters_mob);

  calculation_path = "InterDiffusion";
  auto pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_compute_energies", false));
  auto interdiffu_pst = PST(&spatial, pst_parameters);
  calculation_path = "Calphad";
  auto cc_pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                      Parameter("calculation_path", calculation_path),
                                      Parameter("frequency", frequency));
  auto cc_pst = PST(&spatial, cc_pst_parameters);

  calculation_path = "Diffusion";
  auto diffu_pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_compute_energies", false));
  auto diffu_pst = PST(&spatial, diffu_pst_parameters);

  //-----------------------
  // Problems
  //-----------------------
  Calphad_Problem<AnalyticalIdealSolution<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, mu_var, cc_pst, heat_vars, p_vars, compo_vars);

  auto pp_parameters =
      Parameters(Parameter("Description", "A Mobilities"), Parameter("first_component", "A"),
                 Parameter("last_component", "B"), Parameter("primary_phase", "SOLID"));
  Property_problem<InterDiffusionCoefficient, VARS, PST> a_interdiffusion_mobilities(
      "A inter-diffusion mobilities", pp_parameters, MA, mob_pst_a, compo_vars, heat_vars,
      mobilities);

  Problem<TransientOperator<FECollection, DIM>, VARS, PST> interdiffu_problem(
      "InterDiffusion", interdiffu_oper, compo_vars, {coef_inter}, interdiffu_pst, mu_var, MA);

  Problem<TransientOperator<FECollection, DIM>, VARS, PST> diffu_problem(
      "Fick", diffu_oper, diffu_vars, {coef_fick}, diffu_pst);

  //-----------------------
  // Coupling
  //-----------------------
  auto main_coupling =
      Coupling("Main coupling", cc_problem, a_interdiffusion_mobilities, interdiffu_problem);
  auto checking_coupling = Coupling("Checking coupling", diffu_problem);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& t_initial = 0.0;
  const auto& t_final = 10.0;
  auto time_parameters = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_parameters, main_coupling, checking_coupling);

  time.solve();
  //---------------------------------------
  // Profiling stop
  //---------------------------------------
  Profiling::getInstance().print();

  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  mfem::Mpi::Finalize();
}
