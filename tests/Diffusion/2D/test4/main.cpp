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
  setVerbosity(Verbosity::Debug);
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
  //==========================================
  //======      Inter-diffusion         ======
  //==========================================
  //--- Variables
  const auto& stabCoeff(1.e-9);
  const auto& diffusionCoeff(1.e-8);
  auto td_parameters = Parameters(Parameter("M", diffusionCoeff));

  //--- Integrator : alias definition for the sake of clarity
  using NLFI = ThermoDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                               Diffusion::Constant>;

  //--- Operator definition
  //  Interface thickness

  DiffusionOperator<FECollection, DIM, NLFI, Density::Constant> oper(&spatial, td_parameters,
                                                                     TimeScheme::EulerImplicit);
  oper.overload_diffusion(Parameters(Parameter("D", stabCoeff)));

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
  // Composition
  double dst = 1.e-1;
  auto user_func_solution = std::function<double(const mfem::Vector&, double)>(
      [length, diffusionCoeff, dst](const mfem::Vector& x, double time) {
        const auto xx = x[0];
        auto func =
            0.5 * (1 + std::erf((xx - length / 2) / std::sqrt(4 * diffusionCoeff * 5. * dst)));

        return func;
      });

  // Initial condition for composition
  auto initial_condition = AnalyticalFunctions<DIM>(user_func_solution);
  auto xo = VAR(&spatial, bcs, "O", 2, initial_condition);
  xo.set_additional_information("O", "X");
  auto compo_vars = VARS(xo);

  // Chemical potential
  auto muo = VAR(&spatial, bcs, "muO", 2, 1.);
  muo.set_additional_information("O", "mu");
  auto mu_var = VARS(muo);
  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");
  auto calphad_parameters = Parameters(description_calphad);

  //==========================================
  //==========================================
  //--- Post-Processing
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "InterDiffusion";
  const auto& frequency = 1;
  auto pst_parameters = Parameters(Parameter("main_folder_path", main_folder_path),
                                   Parameter("calculation_path", calculation_path),
                                   Parameter("frequency", frequency));
  auto pst = PST(&spatial, pst_parameters);
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
  Calphad_Problem<AnalyticalIdealSolution<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, mu_var, cc_pst, convergence, heat_vars, p_vars, compo_vars);

  Problem<DiffusionOperator<FECollection, DIM, NLFI, Density::Constant>, VARS, PST> td_problem(
      oper, compo_vars, pst, convergence, mu_var);

  //-----------------------
  // Coupling
  //-----------------------
  auto main_coupling = Coupling("Main coupling", cc_problem, td_problem);

  //---------------------------------------
  // Time discretization
  //---------------------------------------
  const auto& t_initial = 0.0;
  const auto& t_final = 30.0;
  const auto& dt = 0.01;
  auto time_parameters = Parameters(Parameter("initial_time", t_initial),
                                    Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_parameters, main_coupling);

  time.solve();

  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  mfem::Mpi::Finalize();
}
