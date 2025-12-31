/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 1D Inter-diffusion test for a dummy ternary system in a two-phase system
 * @version 0.1
 * @date 2025-03-27
 *
 * @copyright Copyright (c) 2025
 *
 */

//---------------------------------------
// Headers
//---------------------------------------
#include <cmath>
#include <string>
#include <tuple>
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
  const int nb_elem = 100;
  const double pellet_radius = 6.07e-3;

  const std::tuple<int, double>& tuple_of_dimensions =
      std::make_tuple(nb_elem, pellet_radius);  // Number of elements and maximum length

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
  using TH_OPE = DiffusionOperator<FECollection, DIM>;
  using TH_PB = Problem<TH_OPE, VARS, PST>;

  auto temp = VAR(&spatial, thermal_bcs, "T", Glossary::Temperature, level_of_storage, 700.);
  temp.set_additional_information("K", "T");
  auto heat_vars = VARS(temp);
  const double rho(32.e3);  // mol/m3
  const double cp(60.);     // J/mol/K
  const double cond(2.7);   //  W/m/K

  auto src_func = std::function<double(const mfem::Vector&, double)>(
      [pellet_radius](const mfem::Vector& vcoord, double time) {
        const double pl = 10.e4;

        const double radius = std::sqrt(vcoord[0] * vcoord[0]);
        auto chi = 90.;  // inverse neutron diffusion length (0.9cm−1 ->90m-1).
        auto chia = chi * pellet_radius;
        auto I1_chia = std::cyl_bessel_i(1, chia);
        auto chir = chi * radius;  //  (pellet_radius - radius);
        auto I0_chir = std::cyl_bessel_i(0, chir);
        const auto bess = chia * I0_chir / (2. * I1_chia);
        const auto func = pl * bess / (M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });

  std::vector<AnalyticalFunctions<DIM>> source_term;
  source_term.emplace_back(AnalyticalFunctions<DIM>(src_func));
  Coefficient density(Glossary::Concentration, rho);
  Coefficient heat_capacity(Glossary::Cp, cp);
  Coefficient conductivity(Glossary::Conductivity, cond);
  Coefficients coef_heat(density, heat_capacity, conductivity);

  TH_OPE th_operator(spatials, {"Fourier"}, TimeScheme::EulerImplicit, "HeatTimeDerivative",
                     source_term);
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
  auto XC = VAR(&spatial, calphad_bcs, "XCOORD", Glossary::Coordinate, level_of_storage,
                AnalyticalFunctions<DIM>(xcoord));
  XC.set_additional_information("XCOORD");
  auto coord = VARS(XC);
  // Pressure
  auto pres = VAR(&spatial, pressure_bcs, "pressure", Glossary::Pressure, level_of_storage, 50.e5);
  pres.set_additional_information("Pa", "P");
  auto p_vars = VARS(pres);

  // Initial condition for composition
  const double Nmol = 3.005;
  const double initial_compo_a = 2.005 / Nmol;
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

  // MOB LIQUID

  auto lmobA =
      VAR(&spatial, calphad_bcs, "Mo", Glossary::InterDiffusionMobility, level_of_storage, 1.e-8);
  lmobA.set_additional_information("LIQUID", "A", "mob");

  auto lmobB =
      VAR(&spatial, calphad_bcs, "Mu", Glossary::InterDiffusionMobility, level_of_storage, 1.e-9);
  lmobB.set_additional_information("LIQUID", "B", "mob");

  auto lmobC =
      VAR(&spatial, calphad_bcs, "Mpu", Glossary::InterDiffusionMobility, level_of_storage, 1.e-15);
  lmobC.set_additional_information("LIQUID", "C", "mob");

  auto mob_liquid = VARS(lmobA, lmobB, lmobC);

  // Driving forces
  auto dgm_s = VAR(&spatial, calphad_bcs, "DGM_s", Glossary::DrivingForce, level_of_storage, 0.);
  dgm_s.set_additional_information("SOLID", "dgm");
  auto dgm_l = VAR(&spatial, calphad_bcs, "DGM_l", Glossary::DrivingForce, level_of_storage, 0.);
  dgm_l.set_additional_information("LIQUID", "dgm");

  auto nuc_l = VAR(&spatial, calphad_bcs, "NUC_l", Glossary::Nucleus, level_of_storage, 0.);
  nuc_l.set_additional_information("LIQUID", "nucleus");

  // Diffusion chemical potential
  auto dmu_ac =
      VAR(&spatial, calphad_bcs, "dmu_ac", Glossary::ChemicalPotential, level_of_storage, 0.);
  dmu_ac.set_additional_information("A", "dmu");
  auto dmu_bc =
      VAR(&spatial, calphad_bcs, "dmu_bc", Glossary::ChemicalPotential, level_of_storage, 0.);
  dmu_bc.set_additional_information("B", "dmu");

  // Mole fraction of phases
  auto xph_l = VAR(&spatial, calphad_bcs, "xph_l", Glossary::MoleFraction, level_of_storage, 0.);
  xph_l.set_additional_information("LIQUID", "xph");

  // Element molar fraction by phase

  auto xa_s = VAR(&spatial, interdiffu_bcs, "xsA", Glossary::MoleFraction, level_of_storage,
                  initial_compo_a);
  xa_s.set_additional_information("A", "SOLID", "xp");
  auto xb_s = VAR(&spatial, interdiffu_bcs, "xsB", Glossary::MoleFraction, level_of_storage,
                  initial_compo_b);
  xb_s.set_additional_information("B", "SOLID", "xp");
  auto xc_s = VAR(&spatial, interdiffu_bcs, "xsC", Glossary::MoleFraction, level_of_storage,
                  initial_compo_c);
  xc_s.set_additional_information("C", "SOLID", "xp");

  auto xa_l = VAR(&spatial, interdiffu_bcs, "xlA", Glossary::MoleFraction, level_of_storage,
                  initial_compo_a);
  xa_l.set_additional_information("A", "LIQUID", "xp");
  auto xb_l = VAR(&spatial, interdiffu_bcs, "xlB", Glossary::MoleFraction, level_of_storage,
                  initial_compo_b);
  xb_l.set_additional_information("B", "LIQUID", "xp");
  auto xc_l = VAR(&spatial, interdiffu_bcs, "xlC", Glossary::MoleFraction, level_of_storage,
                  initial_compo_c);
  xc_l.set_additional_information("C", "LIQUID", "xp");

  // Gibbs energy
  auto gl = VAR(&spatial, calphad_bcs, "g_l", Glossary::GibbsEnergy, level_of_storage, 0.);
  gl.set_additional_information("LIQUID", "g");
  auto gs = VAR(&spatial, calphad_bcs, "g_s", Glossary::GibbsEnergy, level_of_storage, 0.);
  gs.set_additional_information("SOLID", "g");

  auto calphad_outputs = VARS(mua, mub, muc, mobA, mobB, mobC, dgm_s, dgm_l, dmu_ac, dmu_bc, xph_l,
                              xa_s, xb_s, xc_s, xa_l, xb_l, xc_l, nuc_l, gs, gl);

  auto phi = VAR(&spatial, calphad_bcs, "phi", Glossary::PhaseField, level_of_storage, 1.);
  phi.set_additional_information("SOLID", "phi");
  auto var_phi = VARS(phi);
  // TDB file
  auto description_calphad =
      Parameter("description", "Calphad description for a U-O-Pu ternary system");
  auto element_removed_from_ic = Parameter("element_removed_from_ic", "C");
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

  std::vector<std::string> composition_order{"A", "B", "C"};
  auto input_composition_order = Parameter("InputCompositionOrder", composition_order);

  std::vector<std::string> energy_order{"G", "GM", "H", "HM"};
  auto input_energies_order = Parameter("InputEnergiesOrder", energy_order);

  auto own_mobility_model = Parameter("OwnMobilityModel", false);
  auto own_energy_model = Parameter("OwnEnergyModel", false);

  auto element_removed_from_nn_inputs = Parameter("element_removed_from_nn_inputs", "C");

  auto calphad_parameters =
      Parameters(element_removed_from_ic, neural_network_model_mu, index_neural_network_model_mu,
                 neural_network_model_mob, index_neural_network_model_mob, neural_network_model_nrj,
                 own_mobility_model, own_energy_model, input_composition_factor,
                 input_composition_order, input_energies_order, element_removed_from_nn_inputs) +
      KKS_parameters;

  auto M11 =
      VAR(&spatial, calphad_bcs, "M11", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M11.set_additional_information("A", "inter_mob");
  auto M12 =
      VAR(&spatial, calphad_bcs, "M12", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M12.set_additional_information("B", "inter_mob");

  auto MA = VARS(M11, M12);

  auto M21 =
      VAR(&spatial, calphad_bcs, "M21", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M21.set_additional_information("B", "inter_mob");
  auto M22 =
      VAR(&spatial, calphad_bcs, "M22", Glossary::InterDiffusionMobility, level_of_storage, 0.);
  M22.set_additional_information("A", "inter_mob");

  auto MB = VARS(M21, M22);
  //==========================================
  //======      Melting                 ======
  //==========================================
  using AC_OPE = PhaseFieldOperator<FECollection, DIM>;
  using AC_PB = Problem<AC_OPE, VARS, PST>;

  const auto& lambda = 3. * sigma * epsilon / 2.;
  const auto& omega = 12. * sigma / epsilon;
  auto nuc_parameters =
      Parameters(Parameter("primary_phase", "SOLID"), Parameter("secondary_phase", "LIQUID"),
                 Parameter("melting_factor", 1.));
  auto ac_params = Parameters(Parameter("epsilon", epsilon), Parameter("sigma", sigma),
                              Parameter("lambda", lambda), Parameter("omega", omega)) +
                   nuc_parameters;

  AC_OPE ac_oper(spatials, {"AllenCahn", "MeltingCalphad"}, ac_params, TimeScheme::EulerImplicit,
                 "TimeDerivative");
  ac_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-12), Parameter("abs_tol", 1.e-16)));

  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const double& stabCoeff(1.e-4);

  auto td_parameters = Parameters(Parameter("ScaleCoefficientsByTemperature", true),
                                  Parameter("EnableDiffusionChemicalPotentials", true));

  //--- Operator definition
  // Operator for InterDiffusion equation on A
  Coefficient Dstab(Glossary::Diffusivity, stabCoeff);
  Coefficients coef_pb(Dstab);
  DiffusionOperator<FECollection, DIM> interdiffu_oper_a(
      spatials, {"MassFlux"}, td_parameters, TimeScheme::EulerImplicit, "TimeDerivative");
  interdiffu_oper_a.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-11), Parameter("abs_tol", 5.e-14)));
  // Operator for InterDiffusion equation on B
  DiffusionOperator<FECollection, DIM> interdiffu_oper_b(
      spatials, {"MassFlux"}, td_parameters, TimeScheme::EulerImplicit, "TimeDerivative");
  interdiffu_oper_b.overload_nl_solver(
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
  auto mob_pst_a = PST(&spatial, pst_parameters_mob);

  calculation_path = "HeatTransfer";
  auto pst_parameters_heat =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter),
                 Parameter("enable_compute_energies", false));
  auto heat_pst = PST(&spatial, pst_parameters_heat);

  calculation_path = "Melting";
  auto pst_parameters_ac =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter),
                 Parameter("enable_compute_energies", false));
  auto ac_pst = PST(&spatial, pst_parameters_ac);

  calculation_path = "MobilitiesB";
  auto pst_parameters_mob_b =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto mob_pst_b = PST(&spatial, pst_parameters_mob_b);

  calculation_path = "InterDiffusion_a";
  auto pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter),
                 Parameter("enable_compute_energies", false));
  auto interdiffu_pst = PST(&spatial, pst_parameters);
  calculation_path = "Calphad";
  auto cc_pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  auto cc_pst = PST(&spatial, cc_pst_parameters);

  calculation_path = "InterDiffusion_b";
  auto diffu_pst_parameters =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter),
                 Parameter("enable_compute_energies", false));
  auto interdiffu_pst_b = PST(&spatial, diffu_pst_parameters);

  //-----------------------
  // Problems
  //-----------------------
  //==========================================
  //======      HEAT TRANSFER           ======
  //==========================================
  TH_PB th_problem("Heat tranfer", th_operator, heat_vars, {coef_heat}, heat_pst);

  //==========================================
  //======      CALPHAD                 ======
  //==========================================
  Calphad_Problem<CalphadInformedNeuralNetwork<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, calphad_outputs, cc_pst, heat_vars, p_vars, xa_vars, xb_vars, xc_vars,
      var_phi, coord);

  //==========================================
  //======      Melting                 ======
  //==========================================

  Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
  Coefficient double_well(Glossary::FreeEnergy, Scheme::Implicit, W(omega));
  Coefficient capillary(Glossary::Capillary, lambda);
  Coefficient mobility(Glossary::Mobility, mob);
  Coefficient interpolation(Glossary::InterpolationFunction, Scheme::Implicit, H());

  Coefficients coef_ac(double_well, capillary, mobility, interpolation, grad_energy);

  AC_PB ac_problem("AllenCahn", ac_oper, var_phi, {coef_ac}, ac_pst, calphad_outputs);
  //======================
  // A
  //======================
  auto ppa_parameters =
      Parameters(Parameter("Description", "A Mobilities"), Parameter("first_component", "A"),
                 Parameter("last_component", "C"), Parameter("primary_phase", "SOLID"),
                 Parameter("secondary_phase", "LIQUID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> A_interdiffusion_mobilities(
      "A inter-diffusion mobilities", ppa_parameters, MA, mob_pst_a, xa_vars, xb_vars, heat_vars,
      calphad_outputs, var_phi, mob_liquid);

  Problem<DiffusionOperator<FECollection, DIM>, VARS, PST> interdiffu_problem_a(
      "Interdiffusion A", interdiffu_oper_a, xa_vars, {coef_pb}, interdiffu_pst, calphad_outputs,
      MA, heat_vars);

  //======================
  // B
  //======================
  auto ppb_parameters =
      Parameters(Parameter("Description", "B Mobilities"), Parameter("first_component", "B"),
                 Parameter("last_component", "C"), Parameter("primary_phase", "SOLID"),
                 Parameter("secondary_phase", "LIQUID"));

  Property_problem<InterDiffusionCoefficient, VARS, PST> B_interdiffusion_mobilities(
      "B inter-diffusion mobilities", ppb_parameters, MB, mob_pst_b, xa_vars, xb_vars, heat_vars,
      calphad_outputs, var_phi, mob_liquid);

  Problem<DiffusionOperator<FECollection, DIM>, VARS, PST> interdiffu_problem_b(
      "Interdiffusion B", interdiffu_oper_b, xb_vars, {coef_pb}, interdiffu_pst_b, calphad_outputs,
      MB, heat_vars);

  //-----------------------
  // Coupling
  //-----------------------
  auto th_coupling = Coupling("Thermal coupling", th_problem);
  auto cc_coupling = Coupling("Calphad coupling", cc_problem);
  auto ac_coupling = Coupling("Melting coupling", ac_problem);
  auto diffusion_coupling =
      Coupling("Diffusion coupling", A_interdiffusion_mobilities, B_interdiffusion_mobilities,
               interdiffu_problem_a, interdiffu_problem_b);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& dt = 5.e-3;
  const auto& t_initial = 0.0;
  const auto& t_final = 40.0;
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
