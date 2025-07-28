/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 1D Inter-diffusion test for a ternary system in a two-phase system
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

#include "boost/math/special_functions/bessel.hpp"
#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Quiet);

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
  const int DIM = 1;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PSTCollection = Test<DIM>::PSTCollection;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  using BCS = Test<DIM>::BCS;
  //---------------------------------------
  // Meshing & Boundary Conditions
  //---------------------------------------
  const int refinement_level = 0;
  const int fe_order = 1;

  const std::string& mesh_type = "InlineLineWithSegments";  // type of mesh
  const int NN = 100;
  const double pellet_radius = 6.07e-3;

  const std::tuple<int, double>& tuple_of_dimensions =
      std::make_tuple(NN, pellet_radius);  // Number of elements and maximum length

  SPA spatial(mesh_type, fe_order, refinement_level, tuple_of_dimensions);

  // ##############################
  //     Boundary conditions     //
  // ##############################
  auto interdiffu_bcs =
      BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("upper", 1, "Neumann", 0.));
  auto thermal_bcs =
      BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("upper", 1, "Dirichlet", 700.));
  auto calphad_bcs =
      BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("upper", 1, "Neumann", 0.));
  auto pressure_bcs = BCS(&spatial, Boundary("lower", 0, "Dirichlet", 5.e6),
                          Boundary("upper", 1, "Dirichlet", 5.e6));

  //---------------------------------------
  // Multiphysics coupling scheme
  //---------------------------------------
  const int level_of_storage = 2;

  std::vector<SPA*> spatials{&spatial};
  //==========================================
  //======      HEAT TRANSFER           ======
  //==========================================
  using LHS_NLFI = TimeNLFormIntegrator<VARS>;
  using TH_NLFI =
      HeatNLFormIntegrator<VARS, CoefficientDiscretization::Explicit, Conductivity::Constant>;
  using TH_OPE =
      HeatOperator<FECollection, DIM, TH_NLFI, LHS_NLFI, Density::Constant, HeatCapacity::Constant>;
  using TH_PB = Problem<TH_OPE, VARS, PST>;

  auto temp = VAR(&spatial, thermal_bcs, "T", level_of_storage, 700.);
  temp.set_additional_information("K", "T");
  auto heat_vars = VARS(temp);

  auto src_func = std::function<double(const mfem::Vector&, double)>(
      [pellet_radius](const mfem::Vector& vcoord, double time) {
        const double pl = 10.e4;

        const double radius = std::sqrt(vcoord[0] * vcoord[0]);
        auto chi = 90.;  // inverse neutron diffusion length (0.9cm−1 ->90m-1).
        auto chia = chi * pellet_radius;
        auto I1_chia = boost::math::cyl_bessel_i(1, chia);
        auto chir = chi * radius;  //  (pellet_radius - radius);
        auto I0_chir = boost::math::cyl_bessel_i(0, chir);
        const auto bess = chia * I0_chir / (2. * I1_chia);
        const auto func = pl * bess / (M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });
  const auto& rho(32.e3);  // mol/m3
  const auto& cp(60.);     // J/mol/K
  const auto& cond(2.7);   //  W/m/K

  std::vector<AnalyticalFunctions<DIM>> source_term;
  source_term.emplace_back(AnalyticalFunctions<DIM>(src_func));
  TH_OPE th_operator(spatials, TimeScheme::EulerImplicit, source_term);
  th_operator.overload_density(Parameters(Parameter("rho", rho)));
  th_operator.overload_heat_capacity(Parameters(Parameter("cp", cp)));
  th_operator.overload_conductivity(Parameters(Parameter("lambda", cond)));

  th_operator.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-6), Parameter("abs_tol", 1.e-6)));
  th_operator.overload_solver(HypreSolverType::HYPRE_GMRES);
  th_operator.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU);

  // Interface thickness
  const auto& epsilon(5.e-4);
  // Interfacial energy
  const auto& sigma(6.e-2);
  // Two-phase mobility
  const auto& mob(1.e-4);
  //==========================================
  //======      CALPHAD from TDB        ======
  //==========================================
  //--- Variables
  auto xcoord = std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& vcoord, double time) { return vcoord[0]; });
  auto XC =
      VAR(&spatial, calphad_bcs, "XCOORD", level_of_storage, AnalyticalFunctions<DIM>(xcoord));
  XC.set_additional_information("XCOORD");
  auto coord = VARS(XC);
  // Pressure
  auto pres = VAR(&spatial, pressure_bcs, "pressure", level_of_storage, 50.e5);
  pres.set_additional_information("Pa", "P");
  auto p_vars = VARS(pres);

  // Initial condition for composition
  const double Nmol = 3.005;
  const double initial_compo_o = 2.005 / Nmol;
  const double initial_compo_u = 0.8 / Nmol;
  const double initial_compo_pu = 1. - initial_compo_o - initial_compo_u;

  auto xo = VAR(&spatial, interdiffu_bcs, "O", level_of_storage, initial_compo_o);
  xo.set_additional_information("O", "x");
  auto xo_vars = VARS(xo);

  auto xu = VAR(&spatial, interdiffu_bcs, "U", level_of_storage, initial_compo_u);
  xu.set_additional_information("U", "x");
  auto xu_vars = VARS(xu);

  auto xpu = VAR(&spatial, interdiffu_bcs, "PU", level_of_storage, initial_compo_pu);
  xpu.set_additional_information("PU", "x");
  auto xpu_vars = VARS(xpu);

  // Chemical potential
  auto muo = VAR(&spatial, calphad_bcs, "muO", level_of_storage, 0.);
  muo.set_additional_information("O", "mu");
  auto mu_var = VARS(muo);

  auto muu = VAR(&spatial, calphad_bcs, "muU", level_of_storage, 0.);
  muu.set_additional_information("U", "mu");
  auto muu_var = VARS(muu);

  auto mupu = VAR(&spatial, calphad_bcs, "muPU", level_of_storage, 0.);
  mupu.set_additional_information("PU", "mu");
  auto mupu_var = VARS(mupu);

  // Mobilities
  auto mobO = VAR(&spatial, calphad_bcs, "Mo", level_of_storage, 0.);
  mobO.set_additional_information("SOLID", "O", "mob");

  auto mobU = VAR(&spatial, calphad_bcs, "Mu", level_of_storage, 0.);
  mobU.set_additional_information("SOLID", "U", "mob");

  auto mobPU = VAR(&spatial, calphad_bcs, "Mpu", level_of_storage, 0.);
  mobPU.set_additional_information("SOLID", "PU", "mob");

  // MOB LIQUID

  auto lmobO = VAR(&spatial, calphad_bcs, "Mo", level_of_storage, 1.e-8);
  lmobO.set_additional_information("LIQUID", "O", "mob");

  auto lmobU = VAR(&spatial, calphad_bcs, "Mu", level_of_storage, 1.e-9);
  lmobU.set_additional_information("LIQUID", "U", "mob");

  auto lmobPU = VAR(&spatial, calphad_bcs, "Mpu", level_of_storage, 1.e-15);
  lmobPU.set_additional_information("LIQUID", "PU", "mob");

  auto mob_liquid = VARS(lmobO, lmobU, lmobPU);

  // Driving forces
  auto dgm_s = VAR(&spatial, calphad_bcs, "DGM_s", level_of_storage, 0.);
  dgm_s.set_additional_information("SOLID", "dgm");
  auto dgm_l = VAR(&spatial, calphad_bcs, "DGM_l", level_of_storage, 0.);
  dgm_l.set_additional_information("LIQUID", "dgm");

  auto nuc_l = VAR(&spatial, calphad_bcs, "NUC_l", level_of_storage, 0.);
  nuc_l.set_additional_information("LIQUID", "nucleus");

  // Diffusion chemical potential
  auto dmu_opu = VAR(&spatial, calphad_bcs, "dmu_opu", level_of_storage, 0.);
  dmu_opu.set_additional_information("O", "dmu");
  auto dmu_upu = VAR(&spatial, calphad_bcs, "dmu_upu", level_of_storage, 0.);
  dmu_upu.set_additional_information("U", "dmu");

  // Mole fraction of phases
  auto xph_l = VAR(&spatial, calphad_bcs, "xph_l", level_of_storage, 0.);
  xph_l.set_additional_information("LIQUID", "xph");

  // Element molar fraction by phase

  auto xo_s = VAR(&spatial, interdiffu_bcs, "xsO", level_of_storage, initial_compo_o);
  xo_s.set_additional_information("O", "SOLID", "xp");
  auto xu_s = VAR(&spatial, interdiffu_bcs, "xsU", level_of_storage, initial_compo_u);
  xu_s.set_additional_information("U", "SOLID", "xp");
  auto xpu_s = VAR(&spatial, interdiffu_bcs, "xsPU", level_of_storage, initial_compo_pu);
  xpu_s.set_additional_information("PU", "SOLID", "xp");

  auto xo_l = VAR(&spatial, interdiffu_bcs, "xlO", level_of_storage, initial_compo_o);
  xo_l.set_additional_information("O", "LIQUID", "xp");
  auto xu_l = VAR(&spatial, interdiffu_bcs, "xlU", level_of_storage, initial_compo_u);
  xu_l.set_additional_information("U", "LIQUID", "xp");
  auto xpu_l = VAR(&spatial, interdiffu_bcs, "xlPU", level_of_storage, initial_compo_pu);
  xpu_l.set_additional_information("PU", "LIQUID", "xp");

  // Gibbs energy
  auto gl = VAR(&spatial, calphad_bcs, "g_l", level_of_storage, 0.);
  gl.set_additional_information("LIQUID", "g");
  auto gs = VAR(&spatial, calphad_bcs, "g_s", level_of_storage, 0.);
  gs.set_additional_information("SOLID", "g");

  auto calphad_outputs = VARS(muo, muu, mupu, mobO, mobU, mobPU, dgm_s, dgm_l, dmu_opu, dmu_upu,
                              xph_l, xo_s, xu_s, xpu_s, xo_l, xu_l, xpu_l, nuc_l, gs, gl);

  auto phi = VAR(&spatial, calphad_bcs, "phi", level_of_storage, 1.);
  phi.set_additional_information("SOLID", "phi");
  auto var_phi = VARS(phi);
  // TDB file
  auto description_calphad =
      Parameter("description", "Calphad description for a U-O-Pu ternary system");
  auto element_removed_from_ic = Parameter("element_removed_from_ic", "PU");
  vTuple2StringDouble map_unsuspended_phases = {{"SOLID", "entered", -1}};
  auto unsuspended_phases = Parameter("unsuspended_phases", map_unsuspended_phases);

  auto KKS_secondary_phase = Parameter("KKS_secondary_phase", "LIQUID");
  auto KKS_temperature_increment = Parameter("KKS_temperature_increment", 1.);
  auto KKS_composition_increment = Parameter("KKS_composition_increment", 1.e-7);
  auto KKS_seed = Parameter("KKS_seed", 0.5);
  auto KKS_seed_radius = Parameter("KKS_seed_radius", 1.e-4);
  auto KKS_threshold = Parameter("KKS_threshold", 5.e-2);
  auto KKS_temperature_threshold = Parameter("KKS_temperature_threshold", 2500.);
  auto KKS_freeze_nucleation = Parameter("KKS_freeze_nucleation", true);
  auto KKS_nucleation_started = Parameter("KKS_nucleation_started", false);
  auto KKS_enable_specialized = Parameter("KKS_enable_save_specialized", false);
  auto KKS_nucleation_strategy = Parameter("KKS_nucleation_strategy", "GivenMeltingTemperature");
  auto KKS_given_melting_temperature = Parameter("KKS_given_melting_temperature", 3000.0);
  auto KKS_mobility = Parameter("KKS_mobility", mob);

  auto enable_KKS = Parameter("enable_KKS", true);
  auto KKS_parameters =
      Parameters(KKS_enable_specialized, KKS_secondary_phase, KKS_temperature_increment,
                 KKS_composition_increment, KKS_seed, KKS_seed_radius, KKS_threshold,
                 KKS_temperature_threshold, KKS_freeze_nucleation, KKS_nucleation_started,
                 KKS_mobility, enable_KKS, KKS_nucleation_strategy, KKS_given_melting_temperature);

  // NeuralNetworks

  vTupleStringString CommonNeuralNetwork;
  CommonNeuralNetwork.emplace_back(std::make_tuple("solid_model.pt", "SOLID"));
  CommonNeuralNetwork.emplace_back(std::make_tuple("liquid_model.pt", "LIQUID"));

  auto neural_network_model_mu = Parameter("ChemicalPotentialsNeuralNetwork", CommonNeuralNetwork);
  vTupleStringInt ChemicalPotentialNeuralNetworkIndex;
  ChemicalPotentialNeuralNetworkIndex.emplace_back(std::make_tuple("SOLID", 7));
  ChemicalPotentialNeuralNetworkIndex.emplace_back(std::make_tuple("LIQUID", 4));
  auto index_neural_network_model_mu =
      Parameter("ChemicalPotentialsNeuralNetworkIndex", ChemicalPotentialNeuralNetworkIndex);

  vTupleStringInt MobilitiesNeuralNetworkIndex;
  MobilitiesNeuralNetworkIndex.emplace_back(std::make_tuple("SOLID", 4));

  vTupleStringString MobNeuralNetwork;
  MobNeuralNetwork.emplace_back(std::make_tuple("solid_model.pt", "SOLID"));
  auto neural_network_model_mob = Parameter("MobilitiesNeuralNetwork", MobNeuralNetwork);
  auto index_neural_network_model_mob =
      Parameter("MobilitiesNeuralNetworkIndex", MobilitiesNeuralNetworkIndex);

  auto neural_network_model_nrj = Parameter("EnergiesNeuralNetwork", CommonNeuralNetwork);

  // If the inputs of the model are  moles, not molar fractions (comment in this case)
  auto input_composition_factor = Parameter("InputCompositionFactor", 1.);

  std::vector<std::string> composition_order{"O", "PU", "U"};
  auto input_composition_order = Parameter("InputCompositionOrder", composition_order);

  std::vector<std::string> energy_order{"G", "GM", "H", "HM"};
  auto input_energies_order = Parameter("InputEnergiesOrder", energy_order);

  auto own_mobility_model = Parameter("OwnMobilityModel", false);
  auto own_energy_model = Parameter("OwnEnergyModel", false);

  auto element_removed_from_nn_inputs = Parameter("element_removed_from_nn_inputs", "PU");

  auto calphad_parameters =
      Parameters(element_removed_from_ic, neural_network_model_mu, index_neural_network_model_mu,
                 neural_network_model_mob, index_neural_network_model_mob, neural_network_model_nrj,
                 own_mobility_model, own_energy_model, input_composition_factor,
                 input_composition_order, input_energies_order, element_removed_from_nn_inputs) +
      KKS_parameters;

  auto M11 = VAR(&spatial, calphad_bcs, "M11", level_of_storage, 0.);
  M11.set_additional_information("O", "inter_mob");
  auto M12 = VAR(&spatial, calphad_bcs, "M12", level_of_storage, 0.);
  M12.set_additional_information("U", "inter_mob");

  auto MO = VARS(M11, M12);

  auto M21 = VAR(&spatial, calphad_bcs, "M21", level_of_storage, 0.);
  M21.set_additional_information("U", "inter_mob");
  auto M22 = VAR(&spatial, calphad_bcs, "M22", level_of_storage, 0.);
  M22.set_additional_information("O", "inter_mob");

  auto MU = VARS(M21, M22);
  //==========================================
  //======      Melting                 ======
  //==========================================
  using AC_NLFI =
      AllenCahnCalphadMeltingNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                              ThermodynamicsPotentials::W, Mobility::Constant,
                                              ThermodynamicsPotentials::H>;
  using AC_OPE = PhaseFieldOperator<FECollection, DIM, AC_NLFI, LHS_NLFI>;
  using AC_PB = Problem<AC_OPE, VARS, PST>;

  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto nuc_parameters =
      Parameters(Parameter("primary_phase", "SOLID"), Parameter("secondary_phase", "LIQUID"),
                 Parameter("melting_factor", 1.));
  auto ac_params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                              Parameter("lambda", lambda), Parameter("omega", omega)) +
                   nuc_parameters;

  AC_OPE ac_oper(spatials, ac_params, TimeScheme::EulerImplicit);
  ac_oper.overload_mobility(Parameters(Parameter("mob", mob)));
  ac_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-16)));

  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  using MD = MassDiffusionFluxNLFormIntegrator<VARS>;
  //--- Variables
  const double& stabCoeff(1.e-4);

  auto td_parameters = Parameters(Parameter("ScaleCoefficientsByTemperature", true),
                                  Parameter("EnableDiffusionChemicalPotentials", true));

  //--- Operator definition
  // Operator for InterDiffusion equation on O
  DiffusionOperator<FECollection, DIM, MD, Density::Constant, TimeNLFormIntegrator<VARS>>
      interdiffu_oper_o(spatials, td_parameters, TimeScheme::EulerImplicit);
  interdiffu_oper_o.overload_diffusion(Parameters(Parameter("D", stabCoeff)));
  interdiffu_oper_o.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-11), Parameter("abs_tol", 5.e-14)));
  // Operator for InterDiffusion equation on U
  DiffusionOperator<FECollection, DIM, MD, Density::Constant, TimeNLFormIntegrator<VARS>>
      interdiffu_oper_u(spatials, td_parameters, TimeScheme::EulerImplicit);
  interdiffu_oper_u.overload_diffusion(Parameters(Parameter("D", stabCoeff)));
  interdiffu_oper_u.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-11), Parameter("abs_tol", 5.e-14)));

  //==========================================
  //==========================================
  //--- Post-Processing
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "MobilitiesO";
  bool enable_save_specialized_at_iter = true;
  const auto& frequency = 200;
  auto pst_parameters_mob =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto mob_pst_o = PST(&spatial, pst_parameters_mob);

  calculation_path = "HeatTransfer";
  auto pst_parameters_heat =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto heat_pst = PST(&spatial, pst_parameters_heat);

  calculation_path = "Melting";
  auto pst_parameters_ac =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto ac_pst = PST(&spatial, pst_parameters_ac);

  calculation_path = "MobilitiesU";
  auto pst_parameters_mob_u =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto mob_pst_u = PST(&spatial, pst_parameters_mob_u);

  calculation_path = "InterDiffusion_o";
  auto pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto interdiffu_pst = PST(&spatial, pst_parameters);
  calculation_path = "Calphad";
  auto cc_pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto cc_pst = PST(&spatial, cc_pst_parameters);

  calculation_path = "InterDiffusion_u";
  auto diffu_pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto interdiffu_pst_u = PST(&spatial, diffu_pst_parameters);

  //-----------------------
  // Problems
  //-----------------------
  //==========================================
  //======      HEAT TRANSFER           ======
  //==========================================
  TH_PB th_problem("Heat tranfer", th_operator, heat_vars, heat_pst);

  //==========================================
  //======      CALPHAD                 ======
  //==========================================
  Calphad_Problem<CalphadInformedNeuralNetwork<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, calphad_outputs, cc_pst, heat_vars, p_vars, xo_vars, xu_vars, xpu_vars,
      var_phi, coord);

  //==========================================
  //======      Melting                 ======
  //==========================================
  AC_PB ac_problem("AllenCahn", ac_oper, var_phi, ac_pst, calphad_outputs);
  //======================
  // Oxygen
  //======================
  auto ppo_parameters =
      Parameters(Parameter("Description", "Oxygen Mobilities"), Parameter("first_component", "O"),
                 Parameter("last_component", "Pu"), Parameter("primary_phase", "SOLID"),
                 Parameter("secondary_phase", "LIQUID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> oxygen_interdiffusion_mobilities(
      "Oxygen inter-diffusion mobilities", ppo_parameters, MO, mob_pst_o, xo_vars, xu_vars,
      heat_vars, calphad_outputs, var_phi, mob_liquid);

  Problem<DiffusionOperator<FECollection, DIM, MD, Density::Constant, TimeNLFormIntegrator<VARS>>,
          VARS, PST>
      interdiffu_problem_o("Interdiffusion O", interdiffu_oper_o, xo_vars, interdiffu_pst,
                           calphad_outputs, MO, heat_vars);

  //======================
  // Uranium
  //======================
  auto ppu_parameters =
      Parameters(Parameter("Description", "Oxygen Mobilities"), Parameter("first_component", "U"),
                 Parameter("last_component", "PU"), Parameter("primary_phase", "SOLID"),
                 Parameter("secondary_phase", "LIQUID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> uranium_interdiffusion_mobilities(
      "Uranium inter-diffusion mobilities", ppu_parameters, MU, mob_pst_u, xo_vars, xu_vars,
      heat_vars, calphad_outputs, var_phi, mob_liquid);

  Problem<DiffusionOperator<FECollection, DIM, MD, Density::Constant, TimeNLFormIntegrator<VARS>>,
          VARS, PST>
      interdiffu_problem_u("Interdiffusion U", interdiffu_oper_u, xu_vars, interdiffu_pst_u,
                           calphad_outputs, MU, heat_vars);

  //-----------------------
  // Coupling
  //-----------------------
  auto th_coupling = Coupling("Thermal coupling", th_problem);
  auto cc_coupling = Coupling("Calphad coupling", cc_problem);
  auto ac_coupling = Coupling("Melting coupling", ac_problem);
  auto diffusion_coupling =
      Coupling("Diffusion coupling", oxygen_interdiffusion_mobilities,
               uranium_interdiffusion_mobilities, interdiffu_problem_o, interdiffu_problem_u);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& dt = 5.e-3;
  const auto& t_initial = 0.0;
  const auto& t_final = 20.0;
  auto time_parameters = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_parameters, th_coupling, cc_coupling, ac_coupling,
                                 diffusion_coupling);

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
