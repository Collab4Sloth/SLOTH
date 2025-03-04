/**
 * @file main.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Coupling Calphad(AnalyticalIdealSolution)/HeatTransfer/Allen-Cahn : covering test
 * @version 0.1
 * @date 2025-01-14
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Quiet);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //
  //---------------------------------------
  // Profiling start
  Profiling::getInstance().enable();
  //---------------------------------------

  const auto DIM = 1;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VARS = Variables<FECollection, DIM>;
  using VAR = Variable<FECollection, DIM>;
  using PB = Calphad_Problem<AnalyticalParaboloidForTwoPhase<mfem::Vector>, VARS, PST>;

  // Allen-Cahn
  using NLFI_AC = AllenCahnDiffusionMeltingNLFormIntegrator<
      VARS, ThermodynamicsPotentialDiscretization::Implicit, ThermodynamicsPotentials::W,
      Mobility::Constant, ThermodynamicsPotentials::H>;
  using OPE_AC = AllenCahnOperator<FECollection, DIM, NLFI_AC>;
  using PB_AC = Problem<OPE_AC, VARS, PST>;

  // Thermodiffusion
  using NLFI_TD1 = ThermoDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                                   Diffusion::Constant>;
  using OPE_TD1 = DiffusionOperator<FECollection, DIM, NLFI_TD1, Density::Constant>;
  using PB_TD1 = Problem<OPE_TD1, VARS, PST>;

  // Thermodiffusion
  using NLFI_TD2 = ThermoDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                                   Diffusion::Constant>;
  using OPE_TD2 = DiffusionOperator<FECollection, DIM, NLFI_TD2, Density::Constant>;
  using PB_TD2 = Problem<OPE_TD2, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  auto refinement_level = 0;
  double L = 2;  // 4.65e-3;
  int NN = 3000;
  //   SpatialDiscretization<FECollection, DIM> spatial("GMSH", 1, refinement_level,
  //                                                    "1d_mesh_non_uniform_extra_fine.msh",
  //                                                    false);

  SpatialDiscretization<FECollection, DIM> spatial("InlineLineWithSegments", 1, refinement_level,
                                                   std::make_tuple(NN, L));

  // ##############################
  //     Boundary conditions     //
  // ##############################
  //   auto Calphadboundaries = {Boundary("lower", 0, "Neumann", 0.),
  //                             Boundary("external", 2, "Neumann", 0.),
  //                             Boundary("upper", 1, "Neumann", 0.)};
  auto Calphadboundaries = {Boundary("Boundary_Left", 0, "Neumann", 0.),
                            Boundary("Boundary_Right", 1, "Neumann", 0.)};
  auto Calphadbcs = BoundaryConditions<FECollection, DIM>(&spatial, Calphadboundaries);
  auto Tboundaries = {Boundary("Boundary_Left", 0, "Neumann", 0.),
                      Boundary("Boundary_Right", 1, "Neumann", 0.)};

  //   auto Tboundaries = {Boundary("lower", 0, "Neumann", 0.),
  //                       Boundary("external", 2, "Dirichlet", 750.),
  //                       Boundary("upper", 1, "Neumann", 0.)};
  auto Tbcs = BoundaryConditions<FECollection, DIM>(&spatial, Tboundaries);
  // ####################
  //     parameters    //
  // ####################
  // Heat
  const auto& rho(1.);
  const auto& cp(1.);
  const auto& cond(2.);

  // ############################
  //     variables IC + SRC    //
  // ###########################
  const auto& pellet_radius = 0.00465;
  // Heat
  auto temp = VAR(&spatial, Tbcs, "T", 2, 750.);
  temp.set_additional_information("Temperature", "K");
  auto heat_vars = VARS(temp);
  // AC
  double sigma = 5.129849310676498e-06;

  double L_phi = 1.2;
  // Parser
  mfem::OptionsParser args(argc, argv);
  double W = 1.2e-3;
  args.AddOption(&W, "-e", "--epsilon", "interface thickness");
  double lambda = 129958.33333333333 * W;  // 0.18713999999999997 / W;  // 155.95; p = o(3)
  double mob = L_phi * lambda / std::pow(W, 2);
  args.AddOption(&mob, "-m", "--mobility", "AC mobility");

  double omega = 12. * sigma / W;  // 1 / lambda;//
  args.AddOption(&omega, "-o", "--omega", "double well heigh");

  double kappa = (3. / 2.) * 5.129849310676498e-06 * W;  //  W * W / lambda;/
  args.AddOption(&kappa, "-k", "--kappa", "coeff de gradient");

  double dt = 1e-10;
  args.AddOption(&dt, "-dt", "--timeStep", "time step");
  args.ParseCheck();
  //---------------------------------------

  std::cout << "mobil = " << mob << std::endl;
  std::cout << "omega = " << omega << std::endl;
  std::cout << "kappa = " << kappa << std::endl;
  std::cout << "sigma = " << sigma << std::endl;
  std::cout << "=========== " << std::endl;

  auto params = Parameters(Parameter("epsilon", 0.), Parameter("sigma", 0.),
                           Parameter("lambda", kappa), Parameter("omega", omega));

  const auto& center_x = 0.;

  double epsilon = W;  // std::sqrt(8 * kappa / omega) ;  // std::sqrt(8*kappa/omega);

  auto acIC = std::function<double(const mfem::Vector&, double)>(
      [L, epsilon](const mfem::Vector& v, double time) {
        const double x = v[0];
        const double y = 0.;
        const auto r = std::sqrt(x * x + y * y);
        double func;
        func = 1 - (0.5 + 0.5 * std::tanh((r - 0.5 * L) / epsilon));
        return func;
      });

  auto muO = std::function<double(const mfem::Vector&, double)>(
      [L, epsilon](const mfem::Vector& v, double time) {
        const double x = v[0];
        const double y = 0.;
        const auto r = std::sqrt(x * x + y * y);
        double func;
        func = (0.5 * (0.3 + 0.4) + 0.5 * (0.3 - 0.4) * std::tanh((r - 0.5 * L) / epsilon));
        return func;
      });

  auto muu_ = std::function<double(const mfem::Vector&, double)>(
      [L, epsilon](const mfem::Vector& v, double time) {
        const double x = v[0];
        // const double y = v[1];
        const double y = 0.;
        const auto r = std::sqrt(x * x + y * y);
        double func;
        func = (0.5 * (0.6 + 0.175) + 0.5 * (0.6 - 0.175) * std::tanh((r - 0.5 * L) / epsilon));
        return func;
      });

  auto ac_ic = AnalyticalFunctions<DIM>(acIC);
  auto muO_ = AnalyticalFunctions<DIM>(muO);
  auto muU_ = AnalyticalFunctions<DIM>(muu_);
  auto mass_ic1 = AnalyticalFunctions<DIM>(muO_);
  auto mass_ic2 = AnalyticalFunctions<DIM>(muU_);

  auto vars_ac = VAR(&spatial, Calphadbcs, "phi", 2, ac_ic);
  vars_ac.set_additional_information("Temperature", "K");

  auto ac_vars = VARS(vars_ac);
  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################

  // mu
  auto muo = VAR(&spatial, Calphadbcs, "muO", 2, 0.);
  muo.set_additional_information("O", "mu");
  auto muu = VAR(&spatial, Calphadbcs, "muU", 2, 0.);
  muu.set_additional_information("U", "mu");

  // x
  // SOLID
  auto xso = VAR(&spatial, Calphadbcs, "xsO", 2, 0.);
  xso.set_additional_information("O", "SOLID", "x");
  auto xsu = VAR(&spatial, Calphadbcs, "xsU", 2, 0.);
  xsu.set_additional_information("U", "SOLID", "x");

  // LIQUID
  auto xlo = VAR(&spatial, Calphadbcs, "xlO", 2, 0.);
  xlo.set_additional_information("O", "LIQUID", "x");
  auto xlu = VAR(&spatial, Calphadbcs, "xlU", 2, 0.);
  xlu.set_additional_information("U", "LIQUID", "x");

  // g
  // SOLID
  auto gs = VAR(&spatial, Calphadbcs, "gs", 2, 0.);
  gs.set_additional_information("SOLID", "g");
  auto gl = VAR(&spatial, Calphadbcs, "gl", 2, 0.);
  gl.set_additional_information("LIQUID", "g");

  auto pres = VAR(&spatial, Calphadbcs, "pressure", 2, 50.e5);
  pres.set_additional_information("Pressure", "Pa");
  auto p_vars = VARS(pres);
  // DF if necessary
  // auto DF = VAR(&spatial, Calphadbcs, "DF", 2, 50.e5);
  // DF.set_additional_information("-", "df");
  // // auto DF_vars = VARS(DF);

  auto outputs = VARS(gs, gl, muo, xso, xlo, muu, xsu, xlu);

  auto xo = VAR(Variable<FECollection, DIM>(&spatial, Calphadbcs, "O", 2, mass_ic1));
  xo.set_additional_information("O", "X");
  auto xu = VAR(Variable<FECollection, DIM>(&spatial, Calphadbcs, "U", 2, mass_ic2));
  xu.set_additional_information("U", "X");
  auto compo_vars1 = VARS(xo);
  auto compo_vars2 = VARS(xu);
  // auto compo_vars = VARS(xo,xu);
  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");

  std::map<std::tuple<std::string, std::string>, mfem::real_t> coeff_k;
  coeff_k[std::make_tuple("SOLID", "O")] = 1.0;
  coeff_k[std::make_tuple("SOLID", "U")] = 1.0;
  coeff_k[std::make_tuple("LIQUID", "O")] = 1.0;
  coeff_k[std::make_tuple("LIQUID", "U")] = 1.0;

  auto coeff_ks = Parameter("coefficient_k", coeff_k);

  std::map<std::tuple<std::string, std::string>, mfem::real_t> c_eq;
  c_eq[std::make_tuple("SOLID", "O")] = 0.3;
  c_eq[std::make_tuple("SOLID", "U")] = 0.3;
  c_eq[std::make_tuple("LIQUID", "O")] = 0.4;
  c_eq[std::make_tuple("LIQUID", "U")] = 0.4;
  auto c_eqs = Parameter("equilibrium_composition", c_eq);
  auto calphad_parameters = Parameters(description_calphad, coeff_ks, c_eqs);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  // const std::string& main_folder_path = "Saves";
  const std::string& main_folder_path =
      "Saves_extra_fine/W=" + std::to_string(W) + "_mob=" + std::to_string(mob);
  const auto& level_of_detail = 1;
  const auto& frequency = 100;
  std::string calculation_path = "Calphad";
  auto cpst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto Calphad_pst = PST(&spatial, cpst);
  calculation_path = "Heat";
  auto hpst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto Heat_pst = PST(&spatial, hpst);
  calculation_path = "Allen-Cahn";
  auto acpst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto AC_pst = PST(&spatial, acpst);
  calculation_path = "Diffusion_1";
  auto d1pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto TD1_pst = PST(&spatial, d1pst);
  calculation_path = "Diffusion_2";
  auto d2pst =
      Parameters(Parameter("main_folder_path", main_folder_path),
                 Parameter("calculation_path", calculation_path), Parameter("frequency", frequency),
                 Parameter("level_of_detail", level_of_detail));
  auto TD2_pst = PST(&spatial, d2pst);

  // ####################
  //     Problems      //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  //---------------
  // Calphad
  //---------------
  PB Calphad_pb(calphad_parameters, outputs, Calphad_pst, convergence, heat_vars, p_vars, ac_vars,
                compo_vars1, compo_vars2);

  //---------------
  // AC
  //---------------

  OPE_AC AC_op(&spatial, params, TimeScheme::EulerImplicit);
  AC_op.overload_mobility(Parameters(Parameter("mob", mob)));
  AC_op.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("abs_tol", 1.e-10)));
  PB_AC AC_pb("AC", AC_op, ac_vars, AC_pst, convergence, outputs);

  //---------------
  // AC
  //---------------

  OPE_TD1 TD1_op(&spatial, TimeScheme::EulerImplicit);
  double zero = 0.;
  double d_stab = 0.;  // 1e-2;
  TD1_op.overload_diffusion(
      Parameters(Parameter("D", d_stab), Parameter("Dgs", zero), Parameter("D gl", zero),
                 Parameter("D muo", 1.), Parameter("D xso", zero), Parameter("D xlo", zero),
                 Parameter("D muu", 0.), Parameter("D xsu", zero), Parameter("D xlu", zero)));

  PB_TD1 TD1_pb(TD1_op, compo_vars1, TD1_pst, convergence, outputs);

  OPE_TD2 TD2_op(&spatial, TimeScheme::EulerImplicit);
  // double d_stab2 = 1e-9;
  TD2_op.overload_diffusion(
      Parameters(Parameter("D", d_stab), Parameter("Dgs", zero), Parameter("D gl", zero),
                 Parameter("D muo", 0.), Parameter("D xso", zero), Parameter("D xlo", zero),
                 Parameter("D muu", 0.8), Parameter("D xsu", zero), Parameter("D xlu", zero)));

  PB_TD2 TD2_pb(TD2_op, compo_vars2, TD2_pst, convergence, outputs);

  // ####################
  //     Coupling      //
  // ####################
  auto cc = Coupling("Calphad calculation", Calphad_pb, TD1_pb, TD2_pb, AC_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const auto& t_initial = 0.0;
  const auto& t_final = 8.e-5;
  auto time_params = Parameters(Parameter("initial_time", t_initial),
                                Parameter("final_time", t_final), Parameter("time_step", dt));
  auto time = TimeDiscretization(time_params, cc);

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
