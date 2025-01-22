/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Coupling Calphad(AnalyticalIdealSolution)/HeatTransfer : covering test
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
#include "tests/tests.hpp"
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
  /////////////////////////
  constexpr int DIM = Test<2>::dim;
  using FECollection = Test<2>::FECollection;
  using VARS = Test<2>::VARS;
  using VAR = Test<2>::VAR;
  using PSTCollection = Test<2>::PSTCollection;
  using PST = Test<2>::PST;
  using SPA = Test<2>::SPA;
  using BCS = Test<2>::BCS;
  /////////////////////////
  using PB = Calphad_Problem<AnalyticalIdealSolution<mfem::Vector>, VARS, PST>;

  // Heat
  using NLFI2 =
      HeatNLFormIntegrator<VARS, CoefficientDiscretization::Explicit, Conductivity::Constant>;
  using OPE2 = HeatOperator<FECollection, DIM, NLFI2, Density::Constant, HeatCapacity::Constant>;
  using PB2 = Problem<OPE2, VARS, PST>;

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
  auto Calphadboundaries = {Boundary("lower", 0, "Neumann", 0.),
                            Boundary("external", 2, "Neumann", 0.),
                            Boundary("upper", 1, "Neumann", 0.)};
  auto Calphadbcs = BCS(&spatial, Calphadboundaries);
  auto Tboundaries = {Boundary("lower", 0, "Neumann", 0.),
                      Boundary("external", 2, "Dirichlet", 750.),
                      Boundary("upper", 1, "Neumann", 0.)};
  auto Tbcs = BCS(&spatial, Tboundaries);
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
  auto plmax = 12.e4;
  auto pl = 2.6e4;
  auto pldt = 1.e4;
  auto TimeToIncrease = 0.;

  auto src_func = std::function<double(const mfem::Vector&, double)>(
      [pl, plmax, pldt, TimeToIncrease, pellet_radius](const mfem::Vector& vcoord, double time) {
        double puissance = pl;
        if (time > TimeToIncrease) {
          puissance += (time - TimeToIncrease) * pldt;
          puissance = std::min(plmax, puissance);
        }
        const auto func = puissance / (M_PI * 2. * pellet_radius * pellet_radius);

        return func;
      });

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  auto muo = VAR(&spatial, Calphadbcs, "muO", 2, 0.);
  muo.set_additional_information("O", "mu");
  auto xso = VAR(&spatial, Calphadbcs, "xsO", 2, 0.);
  xso.set_additional_information("O", "SOLUTION", "x");
  auto gs = VAR(&spatial, Calphadbcs, "gs", 2, 0.);
  gs.set_additional_information("SOLUTION", "g");
  auto outputs = VARS(muo, xso, gs);

  auto pres = VAR(&spatial, Calphadbcs, "pressure", 2, 50.e5);
  pres.set_additional_information("Pressure", "Pa");

  auto p_vars = VARS(pres);

  auto xo = VAR(&spatial, Calphadbcs, "O", 2, 0.66);
  xo.set_additional_information("O", "X");
  auto compo_vars = VARS(xo);
  auto description_calphad =
      Parameter("description", "Analytical thermodynamic description for an ideal solution ");
  auto calphad_parameters = Parameters(description_calphad);

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################
  const std::string& main_folder_path = "Saves";
  const auto& level_of_detail = 1;
  const auto& frequency = 1;
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

  // ####################
  //     Problems      //
  // ####################
  const auto crit_cvg_1 = 1.e-12;
  PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg_1);

  //---------------
  // Heat transfer
  //---------------
  auto source_term = AnalyticalFunctions<DIM>(src_func);
  OPE2 Heat_op(&spatial, TimeScheme::EulerImplicit, source_term);
  Heat_op.overload_density(Parameters(Parameter("rho", rho)));
  Heat_op.overload_heat_capacity(Parameters(Parameter("cp", cp)));
  Heat_op.overload_conductivity(Parameters(Parameter("lambda", cond)));
  Heat_op.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("abs_tol", 1.e-10)));
  PB2 Heat_pb("Heat", Heat_op, heat_vars, Heat_pst, convergence);

  //---------------
  // Calphad
  //---------------
  PB Calphad_pb(calphad_parameters, outputs, Calphad_pst, convergence, heat_vars, p_vars,
                compo_vars);

  // ####################
  //     Coupling      //
  // ####################
  auto cc = Coupling("Calphad calculation", Heat_pb, Calphad_pb);

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
