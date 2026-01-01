/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D Inter-diffusion test for a dummy ternary system
 * @version 0.1
 * @date 2025-03-27
 *
 * @copyright Copyright (c) 2025
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
  setVerbosity(Verbosity::Debug);

  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //---------------------------------------
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
  const int fe_order = 1;

  SPA spatial("GMSH", fe_order, refinement_level, "camembert2D.msh", false);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto interdiffu_bcs =
      BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Neumann", 0.),
          Boundary("upper", 1, "Neumann", 0.));
  auto thermal_bcs =
      BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Dirichlet", 835.),
          Boundary("upper", 1, "Neumann", 0.));
  auto calphad_bcs =
      BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Neumann", 0.),
          Boundary("upper", 1, "Neumann", 0.));
  auto pressure_bcs =
      BCS(&spatial, Boundary("lower", 0, "Dirichlet", 5.e6),
          Boundary("external", 2, "Dirichlet", 5.e6), Boundary("upper", 1, "Dirichlet", 5.e6));

  //---------------------------------------
  // Multiphysics coupling scheme
  //---------------------------------------
  const int level_of_storage = 2;

  //==========================================
  //======      CALPHAD from TDB        ======
  //==========================================
  //--- Variables
  // Temperature
  auto parabolic_temp = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& vcoord, [[maybe_unused]] double time) {
        const double Text = 835.;
        const double pellet_radius = 4.65e-3;
        const double puissance = 26087.78539541;
        const double rr = vcoord[0] * vcoord[0] + vcoord[1] * vcoord[1];
        const auto func = Text + puissance * (pellet_radius * pellet_radius - rr) /
                                     (4. * M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });
  auto temp =
      VAR(&spatial, thermal_bcs, "T", Glossary::Temperature, level_of_storage, parabolic_temp);
  temp.set_additional_information("K", "T");
  auto heat_vars = VARS(temp);
  // Pressure
  auto pres = VAR(&spatial, pressure_bcs, "pressure", Glossary::Pressure, level_of_storage, 50.e5);
  pres.set_additional_information("Pa", "P");
  auto p_vars = VARS(pres);

  // Initial condition for composition
  const double Nmol = 2.995;
  const double initial_compo_a = 1.995 / Nmol;
  const double initial_compo_b = 0.8 / Nmol;
  const double initial_compo_c = 1. - initial_compo_a - initial_compo_b;

  auto xa =
      VAR(&spatial, interdiffu_bcs, "A", Glossary::MoleFraction, level_of_storage, initial_compo_a);
  xa.set_additional_information("A", "x");
  auto xa_vars = VARS(xa);

  auto xb =
      VAR(&spatial, interdiffu_bcs, "B", Glossary::MoleFraction, level_of_storage, initial_compo_b);
  xb.set_additional_information("B", "x");
  auto xb_vars = VARS(xb);

  auto xc =
      VAR(&spatial, interdiffu_bcs, "C", Glossary::MoleFraction, level_of_storage, initial_compo_c);
  xc.set_additional_information("C", "x");
  auto xc_vars = VARS(xc);

  // Chemical potential
  auto mua = VAR(&spatial, calphad_bcs, "mua", Glossary::ChemicalPotential, level_of_storage, 0.);
  mua.set_additional_information("A", "mu");

  auto mub = VAR(&spatial, calphad_bcs, "mub", Glossary::ChemicalPotential, level_of_storage, 0.);
  mub.set_additional_information("B", "mu");

  auto muc = VAR(&spatial, calphad_bcs, "muc", Glossary::ChemicalPotential, level_of_storage, 0.);
  muc.set_additional_information("C", "mu");

  // Mobilities
  auto mobA =
      VAR(&spatial, calphad_bcs, "Mo", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  mobA.set_additional_information("SOLID", "A", "mob");

  auto mobB =
      VAR(&spatial, calphad_bcs, "Mu", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  mobB.set_additional_information("SOLID", "B", "mob");

  auto mobC =
      VAR(&spatial, calphad_bcs, "Mpu", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  mobC.set_additional_information("SOLID", "C", "mob");

  auto calphad_outputs = VARS(mua, mub, muc, mobA, mobB, mobC);

  // NeuralNetworks

  auto given_phase = Parameter("GivenPhase", "SOLID");

  vTupleStringString CommonNeuralNetwork;
  CommonNeuralNetwork.emplace_back(std::make_tuple("model.pt", "SOLID"));

  auto neural_network_model_mu = Parameter("ChemicalPotentialsNeuralNetwork", CommonNeuralNetwork);
  vTupleStringString MobilitiesNeuralNetwork;
  MobilitiesNeuralNetwork.emplace_back(std::make_tuple("model.pt", "SOLID"));

  vTupleStringInt MobilitiesNeuralNetworkIndex;
  MobilitiesNeuralNetworkIndex.emplace_back(std::make_tuple("SOLID", 3));

  auto neural_network_model_mob = Parameter("MobilitiesNeuralNetwork", CommonNeuralNetwork);
  auto index_neural_network_model_mob =
      Parameter("MobilitiesNeuralNetworkIndex", MobilitiesNeuralNetworkIndex);
  auto input_composition_factor = Parameter("InputCompositionFactor", Nmol);
  std::vector<std::string> composition_order{"A", "B", "C"};
  auto input_composition_order = Parameter("InputCompositionOrder", composition_order);

  auto own_mobility_model = Parameter("OwnMobilityModel", false);

  auto calphad_parameters =
      Parameters(neural_network_model_mu, neural_network_model_mob, index_neural_network_model_mob,
                 input_composition_factor, input_composition_order, given_phase);

  auto M11 =
      VAR(&spatial, calphad_bcs, "M11", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M11.set_additional_information("A", "inter_mob");
  auto M12 =
      VAR(&spatial, calphad_bcs, "M12", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M12.set_additional_information("B", "inter_mob");
  auto M13 =
      VAR(&spatial, calphad_bcs, "M13", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M13.set_additional_information("C", "inter_mob");

  auto MA = VARS(M11, M12, M13);

  auto M21 =
      VAR(&spatial, calphad_bcs, "M21", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M21.set_additional_information("B", "inter_mob");
  auto M22 =
      VAR(&spatial, calphad_bcs, "M22", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M22.set_additional_information("A", "inter_mob");
  auto M23 =
      VAR(&spatial, calphad_bcs, "M23", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M23.set_additional_information("C", "inter_mob");

  auto MB = VARS(M21, M22, M23);
  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const double& stabCoeff(1.e-7);

  auto td_parameters = Parameters(Parameter("last_component", "C"),
                                  Parameter("ScaleCoefficientsByTemperature", true));

  Coefficient Dstab(Glossary::Diffusivity, stabCoeff);
  Coefficients coef_inter(Dstab);

  //--- Operator definition
  std::vector<SPA*> spatials{&spatial};

  DiffusionOperator<FECollection, DIM> interdiffu_oper_a(
      spatials, {"MassFlux"}, td_parameters, TimeScheme::EulerImplicit, "TimeDerivative");
  interdiffu_oper_a.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("abs_tol", 1.e-20)));

  DiffusionOperator<FECollection, DIM> interdiffu_oper_b(
      spatials, {"MassFlux"}, td_parameters, TimeScheme::EulerImplicit, "TimeDerivative");
  interdiffu_oper_b.overload_nl_solver(
      NLSolverType::NEWTON, Parameters(Parameter("description", "Newton solver "),
                                       Parameter("rel_tol", 1.e-20), Parameter("abs_tol", 1.e-20)));

  //==========================================
  //==========================================
  //--- Post-Processing
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "MobilitiesO";
  const auto& frequency = 1;
  auto pst_parameters_mob = Parameters(Parameter("main_folder_path", main_folder_path),
                                       Parameter("calculation_path", calculation_path),
                                       Parameter("frequency", frequency));
  auto mob_pst_a = PST(&spatial, pst_parameters_mob);

  calculation_path = "MobilitiesU";
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

  //-----------------------
  // Problems
  //-----------------------
  // Calphad
  Calphad_Problem<CalphadInformedNeuralNetwork<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, calphad_outputs, cc_pst, heat_vars, p_vars, xa_vars, xb_vars, xc_vars);

  //======================
  // A
  //======================
  auto ppa_parameters =
      Parameters(Parameter("Description", "A Mobilities"), Parameter("first_component", "A"),
                 Parameter("last_component", "C"), Parameter("primary_phase", "SOLID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> A_interdiffusion_mobilities(
      "A inter-diffusion mobilities", ppa_parameters, MA, mob_pst_a, xa_vars, xb_vars, heat_vars,
      calphad_outputs);

  Problem<DiffusionOperator<FECollection, DIM>, VARS, PST> interdiffu_problem_a(
      "Interdiffusion A", interdiffu_oper_a, xa_vars, {coef_inter}, interdiffu_pst, calphad_outputs,
      MA, heat_vars);

  //======================
  // B
  //======================
  auto ppb_parameters =
      Parameters(Parameter("Description", "B Mobilities"), Parameter("first_component", "B"),
                 Parameter("last_component", "C"), Parameter("primary_phase", "SOLID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> B_interdiffusion_mobilities(
      "B inter-diffusion mobilities", ppb_parameters, MB, mob_pst_b, xa_vars, xb_vars, heat_vars,
      calphad_outputs);

  Problem<DiffusionOperator<FECollection, DIM>, VARS, PST> interdiffu_problem_b(
      "Interdiffusion B", interdiffu_oper_b, xb_vars, {coef_inter}, interdiffu_pst_b,
      calphad_outputs, MB, heat_vars);

  //-----------------------
  // Coupling
  //-----------------------
  auto cc_coupling = Coupling("Calphad coupling", cc_problem);
  auto diffusion_coupling =
      Coupling("Diffusion coupling", A_interdiffusion_mobilities, B_interdiffusion_mobilities,
               interdiffu_problem_a, interdiffu_problem_b);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& dt = 0.5;
  const auto& t_initial = 0.0;
  const auto& t_final = 10.0;
  auto time_parameters = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_parameters, cc_coupling, diffusion_coupling);

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
