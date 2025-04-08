/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D Inter-diffusion test for a ternary system
 * @version 0.1
 * @date 2025-03-27
 *
 * @copyright Copyright (c) 2025
 *
 */

//---------------------------------------
// Headers
//---------------------------------------

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
  const int DIM = 2;
  //   const int DIM = 3;
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

  SPA spatial("GMSH", fe_order, refinement_level, "camembert2D.msh", false);
  //   SPA spatial("GMSH", fe_order, refinement_level, "camembert3D.msh", false);
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

  //   auto interdiffu_bcs = BCS(
  //       &spatial, Boundary("InterPelletPlane", 0, "Neumann", 0.),
  //       Boundary("MidPelletPlane", 1, "Neumann", 0.), Boundary("FrontSurface", 3, "Neumann", 0.),
  //       Boundary("BehindSurface", 2, "Neumann", 0.), Boundary("ExternalSurface", 4, "Neumann",
  //       0.));
  //   auto thermal_bcs =
  //       BCS(&spatial, Boundary("InterPelletPlane", 0, "Neumann", 0.),
  //           Boundary("MidPelletPlane", 1, "Neumann", 0.), Boundary("FrontSurface", 3, "Neumann",
  //           0.), Boundary("BehindSurface", 2, "Neumann", 0.), Boundary("ExternalSurface", 4,
  //           "Dirichlet", 835.));
  //   auto calphad_bcs = BCS(
  //       &spatial, Boundary("InterPelletPlane", 0, "Neumann", 0.),
  //       Boundary("MidPelletPlane", 1, "Neumann", 0.), Boundary("FrontSurface", 3, "Neumann", 0.),
  //       Boundary("BehindSurface", 2, "Neumann", 0.), Boundary("ExternalSurface", 4, "Neumann",
  //       0.));
  //   auto pressure_bcs = BCS(&spatial, Boundary("InterPelletPlane", 0, "Dirichlet", 5.e6),
  //                           Boundary("MidPelletPlane", 1, "Dirichlet", 5.e6),
  //                           Boundary("FrontSurface", 3, "Dirichlet", 5.e6),
  //                           Boundary("BehindSurface", 2, "Dirichlet", 5.e6),
  //                           Boundary("ExternalSurface", 4, "Dirichlet", 5.e6));

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
      [](const mfem::Vector& vcoord, double time) {
        const double Text = 835.;
        const double pellet_radius = 4.65e-3;
        const double puissance = 26087.78539541;
        const double conductivity = 2.;
        const double rr = vcoord[0] * vcoord[0] + vcoord[1] * vcoord[1];
        const auto func = Text + puissance * (pellet_radius * pellet_radius - rr) /
                                     (4. * M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });
  auto temp = VAR(&spatial, thermal_bcs, "T", level_of_storage, parabolic_temp);
  temp.set_additional_information("Temperature", "K");
  auto heat_vars = VARS(temp);
  // Pressure
  auto pres = VAR(&spatial, pressure_bcs, "pressure", level_of_storage, 50.e5);
  pres.set_additional_information("Pressure", "Pa");
  auto p_vars = VARS(pres);

  // Initial condition for composition
  const double Nmol = 2.995;
  const double initial_compo_o = 1.995 / Nmol;
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
  mobO.set_additional_information("C1_MO2", "O", "mob");

  auto mobU = VAR(&spatial, calphad_bcs, "Mu", level_of_storage, 0.);
  mobU.set_additional_information("C1_MO2", "U", "mob");

  auto mobPU = VAR(&spatial, calphad_bcs, "Mpu", level_of_storage, 0.);
  mobPU.set_additional_information("C1_MO2", "PU", "mob");

  auto calphad_outputs = VARS(muo, muu, mupu, mobO, mobU, mobPU);

  // NeuralNetworks
  auto neural_network_model_mu = Parameter("ChemicalPotentialsNeuralNetwork", "model.pt");

  vTupleStringString MobilitiesNeuralNetwork;
  MobilitiesNeuralNetwork.emplace_back(std::make_tuple("model.pt", "C1_MO2"));

  vTupleStringInt MobilitiesNeuralNetworkIndex;
  MobilitiesNeuralNetworkIndex.emplace_back(std::make_tuple("C1_MO2", 3));

  auto neural_network_model_mob = Parameter("MobilitiesNeuralNetwork", MobilitiesNeuralNetwork);
  auto index_neural_network_model_mob =
      Parameter("MobilitiesNeuralNetworkIndex", MobilitiesNeuralNetworkIndex);
  auto input_composition_factor = Parameter("InputCompositionFactor", Nmol);
  std::vector<std::string> composition_order{"O", "U", "PU"};
  auto input_composition_order = Parameter("InputCompositionOrder", composition_order);
  auto calphad_parameters =
      Parameters(neural_network_model_mu, neural_network_model_mob, index_neural_network_model_mob,
                 input_composition_factor, input_composition_order);

  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const double& stabCoeff(1.e-7);

  auto td_parameters = Parameters(Parameter("last_component", "PU"),
                                  Parameter("InterdiffusionScalingByTemperature", true));
  //--- Integrator : alias definition for the sake of clarity
  using InterDiffusionIntegrator =
      TernaryInterDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                            Diffusion::Constant>;
  //--- Operator definition
  // Operator for InterDiffusion equation on O
  DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>
      interdiffu_oper_o(&spatial, td_parameters, TimeScheme::EulerImplicit);
  interdiffu_oper_o.overload_diffusion(Parameters(Parameter("D", stabCoeff)));
  interdiffu_oper_o.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("abs_tol", 1.e-20)));

  // Operator for InterDiffusion equation on U
  DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>
      interdiffu_oper_u(&spatial, td_parameters, TimeScheme::EulerImplicit);
  interdiffu_oper_u.overload_diffusion(Parameters(Parameter("D", stabCoeff)));
  interdiffu_oper_u.overload_nl_solver(
      NLSolverType::NEWTON, Parameters(Parameter("description", "Newton solver "),
                                       Parameter("rel_tol", 1.e-20), Parameter("abs_tol", 1.e-20)));

  //==========================================
  //==========================================
  //--- Post-Processing
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "InterDiffusion_O";
  const auto& frequency = 20;
  auto pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                   Parameter("calculation_path", calculation_path),
                                   Parameter("frequency", frequency));
  auto interdiffu_pst_o = PST(&spatial, pst_parameters);

  calculation_path = "InterDiffusion_U";
  auto pst_parameters_u = Parameters(Parameter("main_folder_path", main_folder_path),
                                     Parameter("calculation_path", calculation_path),
                                     Parameter("frequency", frequency));
  auto interdiffu_pst_u = PST(&spatial, pst_parameters_u);

  calculation_path = "Calphad";
  auto cc_pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                      Parameter("calculation_path", calculation_path),
                                      Parameter("frequency", frequency));
  auto cc_pst = PST(&spatial, cc_pst_parameters);

  //--- Physical Convergence
  const double crit_cvg = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg);

  //-----------------------
  // Problems
  //-----------------------
  // Calphad
  Calphad_Problem<CalphadInformedNeuralNetwork<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, calphad_outputs, cc_pst, convergence, heat_vars, p_vars, xo_vars, xu_vars,
      xpu_vars);

  // InterDiffusion : equation on O
  Problem<DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>, VARS,
          PST>
      interdiffu_problem_o("interdiffusion O", interdiffu_oper_o, xo_vars, interdiffu_pst_o,
                           convergence, calphad_outputs, xu_vars, heat_vars);
  // InterDiffusion : equation on U
  Problem<DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>, VARS,
          PST>
      interdiffu_problem_u("interdiffusion U", interdiffu_oper_u, xu_vars, interdiffu_pst_u,
                           convergence, calphad_outputs, xo_vars, heat_vars);

  //-----------------------
  // Coupling
  //-----------------------
  auto main_coupling =
      Coupling("Main coupling", cc_problem, interdiffu_problem_o, interdiffu_problem_u);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& dt = 0.5;
  const auto& t_initial = 0.0;
  const auto& t_final = 1.0;
  auto time_parameters = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_parameters, main_coupling);

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
