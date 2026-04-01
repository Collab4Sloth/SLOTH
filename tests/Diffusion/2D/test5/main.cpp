/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D Inter-diffusion test for a ternary
 * @version 0.1
 * @date 2025-03-28
 *
 * @copyright Copyright (c) 2025
 *
 */

//---------------------------------------
// Headers
//---------------------------------------

#include <string>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

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
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                     Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};

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
            0.5 * (1 + 0.5 * std::erf((xx - length / 2) / std::sqrt(4 * diffusionCoeff * 5. * dt)));

        return func;
      });

  auto user_func_init_b = std::function<double(const mfem::Vector&, double)>(
      []([[maybe_unused]] const mfem::Vector& x, [[maybe_unused]] double time) {
        auto func = 0.1;

        return func;
      });

  auto initial_compo = AnalyticalFunctions<DIM>(user_func_init);
  auto initial_compo_b = AnalyticalFunctions<DIM>(user_func_init_b);

  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const auto& stabCoeff(1.e-7);
  auto td_parameters = Parameters(Parameter("last_component", "C"),
                                  Parameter("ScaleCoefficientsByTemperature", false),
                                  Parameter("EnableDiffusionChemicalPotentials", true));
  auto mobA = VAR(&spatial, bcs, "Ma", Glossary::InterDiffusionMobility, 2, diffusionCoeff);
  mobA.set_additional_information("SOLID", "A", "mob");
  auto mobB = VAR(&spatial, bcs, "Mb", Glossary::InterDiffusionMobility, 2, diffusionCoeff);
  mobB.set_additional_information("SOLID", "B", "mob");
  auto mobC = VAR(&spatial, bcs, "Mc", Glossary::InterDiffusionMobility, 2, diffusionCoeff);
  mobC.set_additional_information("SOLID", "C", "mob");
  auto mobilities = VARS(mobA, mobB, mobC);

  auto mo1 = VAR(&spatial, bcs, "Mo1", Glossary::InterDiffusionMobility, 2, 0.);
  mo1.set_additional_information("A", "inter_mob");
  auto mo2 = VAR(&spatial, bcs, "Mo2", Glossary::InterDiffusionMobility, 2, 0.);
  mo2.set_additional_information("B", "inter_mob");
  auto MA = VARS(mo1, mo2);

  auto mu1 = VAR(&spatial, bcs, "Mu1", Glossary::InterDiffusionMobility, 2, 0.);
  mu1.set_additional_information("B", "inter_mob");
  auto mu2 = VAR(&spatial, bcs, "Mu2", Glossary::InterDiffusionMobility, 2, 0.);
  mu2.set_additional_information("A", "inter_mob");
  auto MB = VARS(mu1, mu2);

  Coefficient explicit_time_A(Glossary::ExplicitTime_A, 1.0);
  Coefficient explicit_time_B(Glossary::ExplicitTime_B, 1.0);
  Coefficient Dstab(Glossary::Diffusivity, stabCoeff);
  Coefficients coef_inter(Dstab, explicit_time_A, explicit_time_B);
  TransientOperator<FECollection, DIM> interdiffu_oper(spatials, {"MassFlux"}, td_parameters,
                                                       TimeScheme::RungeKutta4, "TimeDerivative");

  TransientOperator<FECollection, DIM> interdiffu_oper_b(spatials, {"MassFlux"}, td_parameters,
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
  auto xa = VAR(&spatial, bcs, "A", Glossary::MoleFraction, 2, initial_compo);
  xa.set_additional_information("A", "x");
  auto var_a = VARS(xa);
  auto xb = VAR(&spatial, bcs, "B", Glossary::MoleFraction, 2, initial_compo_b);
  xb.set_additional_information("B", "x");
  auto var_b = VARS(xb);

  // Chemical potential
  auto mua = VAR(&spatial, bcs, "muA", Glossary::ChemicalPotential, 2, 0.);
  mua.set_additional_information("A", "dmu");
  auto mub = VAR(&spatial, bcs, "muB", Glossary::ChemicalPotential, 2, 0.);
  mub.set_additional_information("B", "dmu");

  auto mu_var = VARS(mua, mub);

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

  calculation_path = "MobilitiesB";
  auto pst_parameters_mob_b = Parameters(Parameter("main_folder_path", main_folder_path),
                                         Parameter("calculation_path", calculation_path),
                                         Parameter("frequency", frequency));
  auto mob_pst_b = PST(&spatial, pst_parameters_mob_b);

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

  calculation_path = "InterDiffusion_b";
  auto diffu_pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_compute_energies", false));
  auto interdiffu_pst_b = PST(&spatial, diffu_pst_parameters);

  //======================
  // Calphad
  //======================
  Calphad_Problem<AnalyticalIdealSolution<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, mu_var, cc_pst, heat_vars, p_vars, var_a, var_b);

  auto ppa_parameters =
      Parameters(Parameter("Description", "A Mobilities"), Parameter("first_component", "A"),
                 Parameter("last_component", "C"), Parameter("primary_phase", "SOLID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> a_interdiffusion_mobilities(
      "A inter-diffusion mobilities", ppa_parameters, MA, mob_pst_a, var_a, var_b, heat_vars,
      mobilities);

  Problem<TransientOperator<FECollection, DIM>, VARS, PST> interdiffu_problem_a(
      "Interdiffusion A", interdiffu_oper, var_a, {coef_inter}, interdiffu_pst, mu_var, MA);

  auto ppb_parameters =
      Parameters(Parameter("Description", "B Mobilities"), Parameter("first_component", "B"),
                 Parameter("last_component", "C"), Parameter("primary_phase", "SOLID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> b_interdiffusion_mobilities(
      "B inter-diffusion mobilities", ppb_parameters, MB, mob_pst_b, var_a, var_b, heat_vars,
      mobilities);

  Problem<TransientOperator<FECollection, DIM>, VARS, PST> interdiffu_problem_b(
      "Interdiffusion B", interdiffu_oper_b, var_b, {coef_inter}, interdiffu_pst_b, mu_var, MB);

  //-----------------------
  // Coupling
  //-----------------------
  auto cc_coupling = Coupling("Calphad coupling", cc_problem);
  auto diffusion_coupling =
      Coupling("Diffusion coupling", a_interdiffusion_mobilities, b_interdiffusion_mobilities,
               interdiffu_problem_a, interdiffu_problem_b);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& t_initial = 0.0;
  const auto& t_final = 10.0;
  auto time_parameters = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_parameters, cc_coupling, diffusion_coupling);

  time.solve();

  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  mfem::Mpi::Finalize();
}
