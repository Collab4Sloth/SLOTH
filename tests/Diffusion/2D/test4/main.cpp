/**
 * @file main.cpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief 2D Inter-diffusion test for a binary system
 * @version 0.1
 * @date 2025-03-13
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
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //---------------------------------------
  // Common aliases
  //---------------------------------------
  const int DIM = 2;
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
  auto length = 1.e-3;
  auto nb_fe = 30;

  SPA spatial("InlineSquareWithQuadrangles", 1, refinement_level,
              std::make_tuple(nb_fe, nb_fe, length, length));
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
      [length, diffusionCoeff, dt](const mfem::Vector& x, double time) {
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

  auto diffu_vars =
      VARS(VAR(&spatial, bcs, "c", level_of_storage, initial_compo, analytical_compo));
  //--- Integrator : alias definition for the sake of clarity
  using DiffusionIntegrator =
      DiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit, Diffusion::Constant>;

  //--- Operator definition
  DiffusionOperator<FECollection, DIM, DiffusionIntegrator, Density::Constant> diffu_oper(
      &spatial, TimeScheme::RungeKutta4);
  diffu_oper.overload_diffusion(Parameters(Parameter("D", diffusionCoeff)));

  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const auto& stabCoeff(1.e-7);
  auto td_parameters = Parameters(Parameter("last_component", "U"));
  //   Parameter("MO", diffusionCoeff), Parameter("MU", diffusionCoeff),
  auto mobO = VAR(&spatial, bcs, "Mo", 2, diffusionCoeff);
  mobO.set_additional_information("C1_MO2", "O", "mob");
  auto mobU = VAR(&spatial, bcs, "Mu", 2, diffusionCoeff);
  mobU.set_additional_information("C1_MO2", "U", "mob");
  auto mobo_var = VARS(mobO);
  auto mobu_var = VARS(mobU);
  //--- Integrator : alias definition for the sake of clarity
  using InterDiffusionIntegrator =
      BinaryInterDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                           Diffusion::Constant>;
  //--- Operator definition
  DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant> interdiffu_oper(
      &spatial, td_parameters, TimeScheme::RungeKutta4);
  interdiffu_oper.overload_diffusion(Parameters(Parameter("D", stabCoeff)));

  //==========================================
  //======      CALPHAD Analytical      ======
  //==========================================
  //--- Variables
  // Temperature
  auto temp = VAR(&spatial, bcs, "T", level_of_storage, 1. / Physical::R);
  temp.set_additional_information("Temperature", "K");
  auto heat_vars = VARS(temp);
  // Pressure
  auto pres = VAR(&spatial, bcs, "pressure", level_of_storage, 1.);
  pres.set_additional_information("Pressure", "Pa");
  auto p_vars = VARS(pres);

  // Initial condition for composition
  auto xo = VAR(&spatial, bcs, "O", 2, initial_compo, analytical_compo);
  xo.set_additional_information("O", "X");
  auto compo_vars = VARS(xo);

  // Chemical potential
  auto muo = VAR(&spatial, bcs, "muO", 2, 0.);
  muo.set_additional_information("O", "mu");
  auto mu_var = VARS(muo);
  auto muu = VAR(&spatial, bcs, "muU", 2, 0.);
  muu.set_additional_information("U", "mu");
  auto muu_var = VARS(muu);

  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");
  auto calphad_parameters = Parameters(description_calphad);

  //==========================================
  //==========================================
  //--- Post-Processing
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "InterDiffusion";
  const auto& frequency = 50;
  auto pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                   Parameter("calculation_path", calculation_path),
                                   Parameter("frequency", frequency));
  auto interdiffu_pst = PST(&spatial, pst_parameters);
  calculation_path = "Calphad";
  auto cc_pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                      Parameter("calculation_path", calculation_path),
                                      Parameter("frequency", frequency));
  auto cc_pst = PST(&spatial, cc_pst_parameters);

  calculation_path = "Diffusion";
  auto diffu_pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                         Parameter("calculation_path", calculation_path),
                                         Parameter("frequency", frequency));
  auto diffu_pst = PST(&spatial, diffu_pst_parameters);

  //--- Physical Convergence
  const double crit_cvg = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg);

  //-----------------------
  // Problems
  //-----------------------
  Calphad_Problem<AnalyticalIdealSolution<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, mu_var, cc_pst, convergence, heat_vars, p_vars, compo_vars);

  Problem<DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>, VARS,
          PST>
      interdiffu_problem(interdiffu_oper, compo_vars, interdiffu_pst, convergence, mu_var, muu_var,
                         mobo_var, mobu_var);

  Problem<DiffusionOperator<FECollection, DIM, DiffusionIntegrator, Density::Constant>, VARS, PST>
      diffu_problem(diffu_oper, diffu_vars, diffu_pst, convergence);

  //-----------------------
  // Coupling
  //-----------------------
  auto main_coupling = Coupling("Main coupling", cc_problem, interdiffu_problem);
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
  // Finalize MPI
  //---------------------------------------
  mfem::Mpi::Finalize();
}
