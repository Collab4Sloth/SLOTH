/**
 * @file main.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Allen-Cahn planar front with Dirichlet boundary conditions, AMR and mesh
 *        rebalancing (2 MPI ranks)
 * @version 0.1
 * @date 2026-10-22
 *
 * @copyright Copyright (c) 2026
 *
 * Non-regression test of the essential (Dirichlet) true dofs after a mesh modification:
 * phi = 0 on the left boundary, phi = 1 on the right boundary, planar tanh front at x = 0.5,
 * AMR driven by phi (Kelly) with rebalancing after each AMR pass. Every refinement,
 * derefinement and rebalancing changes the true dof numbering, so the list of essential
 * true dofs must be recomputed before each solve.
 *
 * At the end of the computation, the solution is checked on the Dirichlet boundaries of the
 * final mesh: the test fails (non-zero exit code) if |phi - g| > 1e-10 on any of their true
 * dofs.
 */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------
  setVerbosity(Verbosity::Normal);
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  //---------------------------------------
  // Profiling
  Profiling::getInstance().enable();
  //---------------------------------------
  using namespace Sloth2D;

  bool success = true;
  {  // scope: every SLOTH/MFEM object is destroyed before MPI_Finalize
    // ###########################################
    // ###########################################
    //         Spatial Discretization           //
    // ###########################################
    // ###########################################
    const int order_fe = 1;
    const int refinement_level = 0;
    const int nx = 8;  // coarse mesh, refined by AMR in the interface
    const int ny = 8;
    const double lx = 1.;
    const double ly = 1.;
    const auto tuple_of_dimensions = std::make_tuple(nx, ny, lx, ly);
    // Nonconforming mesh (required for derefinement and rebalancing)
    SPA spatial("InlineSquareWithQuadrangles", order_fe, refinement_level, tuple_of_dimensions,
                true);

    // ##############################
    //     Boundary conditions     //
    // ##############################
    const double phi_left = 0.;
    const double phi_right = 1.;
    auto boundaries_phi = {
        Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", phi_right),
        Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", phi_left)};
    auto bcs_phi = BCS(&spatial, boundaries_phi);

    // ###########################################
    // ###########################################
    //            Physical models               //
    // ###########################################
    // ###########################################
    const double epsilon(0.1);  // interface thickness
    const double sigma(1.);     // interfacial energy
    const double mob(1.);       // mobility
    const double lambda = 1.5 * sigma * epsilon;
    const double omega = 12. * sigma / epsilon;

    Coefficient grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda));
    Coefficient double_well(Glossary::FreeEnergy, Scheme::Explicit, W(omega));
    Coefficient capillary(Glossary::Capillary, lambda);
    Coefficient mobility(Glossary::Mobility, mob);
    Coefficients coef_pb1(double_well, capillary, mobility, grad_energy);

    // ####################
    //     variables     //
    // ####################
    // Planar front at x = 0.5 (sharp initial profile, relaxed by Allen-Cahn)
    const double center_x = 0.;
    const double center_y = 0.;
    const double a_x = 1.;
    const double a_y = 0.;
    const double thickness = 5.e-3;
    const double radius = 5.e-1;
    auto initial_condition =
        AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, center_x, center_y,
                                 a_x, a_y, thickness, radius);

    auto vars = VARS(VAR(&spatial, bcs_phi, "phi", Glossary::PhaseField, 2, initial_condition));

    // ###########################################
    // ###########################################
    //      Post-processing                     //
    // ###########################################
    // ###########################################
    const std::string main_folder_path = "Saves";
    const int level_of_detail = 1;
    const int frequency = 1;
    std::string calculation_path = "Problem1";
    auto p_pst = Parameters(Parameter("main_folder_path", main_folder_path),
                            Parameter("calculation_path", calculation_path),
                            Parameter("frequency", frequency),
                            Parameter("level_of_detail", level_of_detail));
    auto pst = PST(&spatial, p_pst);

    // ####################
    //     operators     //
    // ####################
    std::vector<SPA*> spatials{&spatial};
    TransientOPE oper(spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
    oper.overload_nl_solver(NLSolverType::NEWTON,
                            Parameters(Parameter("description", "Newton solver "),
                                       Parameter("print_level", 1), Parameter("rel_tol", 1.e-12),
                                       Parameter("abs_tol", 1.e-12), Parameter("iter_max", 100)));
    oper.overload_solver(HypreSolverType::HYPRE_GMRES);
    oper.overload_preconditioner(HyprePreconditionerType::HYPRE_ILU);

    TransientPB problem1(oper, vars, {coef_pb1}, pst);

    // AMR on phi (Kelly estimator) with rebalancing after each AMR pass
    mfem::ConstantCoefficient amr_coef{1.0};
    mfem::DiffusionIntegrator amr_integ{amr_coef};
    SlothErrorEstimators estimator(ErrorEstimatorType::KELLY, &amr_integ);

    SingleVariableAMR<VARS> amr(*spatial.get_mesh(), spatial.is_nc_simplices(), "phi");
    auto amr_params =
        Parameters(Parameter("max_elem_error", 1.e-5), Parameter("amr_max_level", 4),
                   Parameter("nc_limit", 1), Parameter("max_preref_cycles", 4),
                   Parameter("scale_down_factor", 0.25), Parameter("rebalance", true));
    amr.SetCriteria(&estimator, amr_params);
    problem1.set_amr(&amr);

    // Coupling 1
    auto cc = Coupling("AllenCahn Coupling", problem1);

    // ###########################################
    // ###########################################
    //            Time-integration              //
    // ###########################################
    // ###########################################
    const double t_initial = 0.0;
    const double t_final = 5.e-2;
    const double dt = 1.e-2;
    auto time_params = Parameters(Parameter("initial_time", t_initial),
                                  Parameter("final_time", t_final), Parameter("time_step", dt));
    auto time = TimeDiscretization(time_params, cc);

    time.solve();

    // ###########################################
    // ###########################################
    //   Check of the Dirichlet boundaries      //
    // ###########################################
    // ###########################################
    auto& phi = vars[0];
    mfem::ParFiniteElementSpace* fes = phi.get_fespace();
    const mfem::Vector u = phi.get_unknown();
    MFEM_VERIFY(u.Size() == fes->GetTrueVSize(),
                "AllenCahn2Dtest1AMR: unknown vector and finite element space sizes differ.");

    const int nb_bdr = fes->GetParMesh()->bdr_attributes.Max();
    const std::vector<std::pair<int, double>> dirichlet = {{1, phi_right}, {3, phi_left}};
    double local_err = 0.;
    for (const auto& [index, value] : dirichlet) {
      mfem::Array<int> marker(nb_bdr);
      marker = 0;
      marker[index] = 1;
      mfem::Array<int> dofs;
      fes->GetEssentialTrueDofs(marker, dofs);
      for (int k = 0; k < dofs.Size(); k++) {
        local_err = std::max(local_err, std::abs(u(dofs[k]) - value));
      }
    }
    double global_err = 0.;
    MPI_Allreduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_MAX, fes->GetComm());
    success = (global_err <= 1.e-10);
    if (mfem::Mpi::WorldRank() == 0) {
      std::cout << "AllenCahn2Dtest1AMR: max |phi - g| on the Dirichlet boundaries = " << global_err
                << (success ? " (OK)" : " (FAILED)") << std::endl;
    }
    //---------------------------------------
    // Profiling stop
    //---------------------------------------
    Profiling::getInstance().print();
  }
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  MPI_Finalize();
  //---------------------------------------
  return success ? 0 : 1;
}
