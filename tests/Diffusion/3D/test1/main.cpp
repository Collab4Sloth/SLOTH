/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
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
//   std::vector<int> vecnbTemp = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
//   std::vector<int> vecnbc = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    std::vector<int> vecnbTemp = {257};
  std::vector<int> vecnbc = {257};
  std::vector<double> vec_duration;
  std::vector<double> vec_temp;
  std::vector<double> vec_c;
  for (int kk = 0; kk < vecnbTemp.size(); kk++) {
    for (int kkk = 0; kkk < vecnbc.size(); kkk++) {
      auto start = std::chrono::high_resolution_clock::now();
      int NN = 200;
      //   SPA spatial("GMSH", fe_order, refinement_level, "camembert2D.msh", false);

      SPA spatial("InlineLineWithSegments", 1, refinement_level, std::make_tuple(NN, 4.65e-3));
      int nbTemp = vecnbTemp[kk];
      int nbc = vecnbc[kkk];

      // ##############################
      //     Boundary conditions     //
      // ##############################

      //   auto interdiffu_bcs =
      //       BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Neumann",
      //       0.),
      //           Boundary("upper", 1, "Neumann", 0.));
      //   auto thermal_bcs =
      //       BCS(&spatial, Boundary("lower", 0, "Neumann", 0.),
      //           Boundary("external", 2, "Dirichlet", 835.), Boundary("upper", 1, "Neumann", 0.));
      //   auto calphad_bcs =
      //       BCS(&spatial, Boundary("lower", 0, "Neumann", 0.), Boundary("external", 2, "Neumann",
      //       0.),
      //           Boundary("upper", 1, "Neumann", 0.));
      //   auto pressure_bcs =
      //       BCS(&spatial, Boundary("lower", 0, "Dirichlet", 5.e6),
      //           Boundary("external", 2, "Dirichlet", 5.e6), Boundary("upper", 1,
      //           "Dirichlet", 5.e6));
      auto interdiffu_bcs =
          BCS(&spatial, Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.));
      auto thermal_bcs = BCS(&spatial, Boundary("left", 0, "Neumann", 0.),
                             Boundary("right", 1, "Dirichlet", 835.));
      auto calphad_bcs =
          BCS(&spatial, Boundary("left", 0, "Neumann", 0.), Boundary("right", 1, "Neumann", 0.));
      auto pressure_bcs = BCS(&spatial, Boundary("left", 0, "Dirichlet", 5.e6),
                              Boundary("right", 1, "Dirichlet", 5.e6));

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
            // const double rr = vcoord[0] * vcoord[0] + vcoord[1] * vcoord[1];
            const double rr = vcoord[0] * vcoord[0];

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
      const double initial_compo = 1.995 / Nmol;
      auto xo = VAR(&spatial, interdiffu_bcs, "O", level_of_storage, initial_compo);
      xo.set_additional_information("O", "x");
      auto xo_vars = VARS(xo);

      auto nu = VAR(&spatial, calphad_bcs, "U", level_of_storage, Nmol);
      nu.set_additional_information("U", "N");
      auto nu_vars = VARS(nu);

      // Chemical potential
      auto muo = VAR(&spatial, calphad_bcs, "muO", level_of_storage, 0.);
      muo.set_additional_information("O", "mu");
      auto mu_var = VARS(muo);
      auto muu = VAR(&spatial, calphad_bcs, "muU", level_of_storage, 0.);
      muu.set_additional_information("U", "mu");
      auto muu_var = VARS(muu);
      // Mobilities
      auto mobO = VAR(&spatial, calphad_bcs, "Mo", level_of_storage, 0.);
      mobO.set_additional_information("C1_MO2", "O", "mob");
      auto mobU = VAR(&spatial, calphad_bcs, "Mu", level_of_storage, 0.);
      mobU.set_additional_information("C1_MO2", "U", "mob");
      auto calphad_outputs = VARS(muo, muu, mobO, mobU);

      // TDB file
      auto description_calphad =
          Parameter("description", "Calphad description for a U-O binary system");
      auto element_removed_from_ic = Parameter("element_removed_from_ic", "U");
      std::map<std::string, double> map_unsuspended_phases = {{"C1_MO2", -1}};
      auto unsuspended_phases = Parameter("unsuspended_phases", map_unsuspended_phases);
      //   std::map<std::tuple<std::string, std::string>, std::tuple<double, double>>
      //   map_reference_states =
      //       {{std::make_tuple("O", "GAS"), std::make_tuple(-1, 1.e5)}};
      //   auto reference_states = Parameter("reference_states", map_reference_states);
      std::string hffilename =
          "database/NT_" + std::to_string(nbTemp) + "_Nc_" + std::to_string(nbc) + ".h5";
      auto paramh5file = Parameter("data filename", hffilename);
      auto calphad_parameters = Parameters(description_calphad, element_removed_from_ic,
                                           unsuspended_phases, paramh5file);

      //==========================================
      //======      Inter-diffusion         ======
      //==========================================
      //--- Variables
      const double& stabCoeff(1.e-7);

      auto td_parameters = Parameters(Parameter("last_component", "U"),
                                      Parameter("InterdiffusionScalingByTemperature", true));
      //--- Integrator : alias definition for the sake of clarity
      using InterDiffusionIntegrator =
          BinaryInterDiffusionNLFormIntegrator<VARS, CoefficientDiscretization::Explicit,
                                               Diffusion::Constant>;
      //--- Operator definition
      DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>
          interdiffu_oper(&spatial, td_parameters, TimeScheme::EulerImplicit);
      interdiffu_oper.overload_diffusion(Parameters(Parameter("D", stabCoeff)));

      //==========================================
      //==========================================
      //--- Post-Processing
      const std::string& main_folder_path = "Saves";
      std::string Savefiles = "SavesNT_" + std::to_string(nbTemp) + "_Nc_" + std::to_string(nbc);
      std::string calculation_path = "InterDiffusion";
      const auto& frequency = 10;
      auto pst_parameters = Parameters(Parameter("main_folder_path", Savefiles),
                                       Parameter("calculation_path", calculation_path),
                                       Parameter("frequency", frequency));
      auto interdiffu_pst = PST(&spatial, pst_parameters);
      calculation_path = "Calphad";
      auto cc_pst_parameters = Parameters(Parameter("main_folder_path", Savefiles),
                                          Parameter("calculation_path", calculation_path),
                                          Parameter("frequency", frequency));
      auto cc_pst = PST(&spatial, cc_pst_parameters);

      //--- Physical Convergence
      const double crit_cvg = 1.e-12;
      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, crit_cvg);

      //-----------------------
      // Problems
      //-----------------------
  Calphad_Problem<MultiParamsTabulation<mfem::Vector>, VARS, PST> cc_problem(
      calphad_parameters, calphad_outputs, cc_pst, convergence, heat_vars, p_vars, xo_vars,
      nu_vars);

      Problem<DiffusionOperator<FECollection, DIM, InterDiffusionIntegrator, Density::Constant>,
              VARS, PST>
          interdiffu_problem(interdiffu_oper, xo_vars, interdiffu_pst, convergence, calphad_outputs,
                             heat_vars);

      //-----------------------
      // Coupling
      //-----------------------
      auto main_coupling = Coupling("Main coupling", cc_problem, interdiffu_problem);

      //---------------------------------------
      // Time discretization
      //---------------------------------------
      const auto& dt = 0.5;
      const auto& t_initial = 0.0;
      const auto& t_final = 500.0;
      auto time_parameters =
          Parameters(Parameter("initial_time", t_initial), Parameter("final_time", t_final),
                     Parameter("time_step", dt));
      auto time = TimeDiscretization(time_parameters, main_coupling);

      time.solve();
      //---------------------------------------
      // Profiling stop
      //---------------------------------------
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - start;
      vec_duration.emplace_back(std::move(duration.count()));
      vec_temp.emplace_back(std::move(nbTemp));
      vec_c.emplace_back(std::move(nbc));
      std::cout << "Temps écoulé : " << duration.count() << " secondes" << std::endl;
      Profiling::getInstance().print();
    }
  }
  std::ofstream fichierCSV("perf_time.csv");
  if (!fichierCSV.is_open()) {
    std::cerr << "error" << std::endl;
    return 1;
  }

  // Écrire les en-têtes de colonnes
  fichierCSV << "T,x_O,time\n";

  // Écrire les données des vecteurs dans le fichier CSV
  for (size_t i = 0; i < vec_duration.size(); ++i) {
    fichierCSV << vec_temp[i] << "," << vec_c[i] << "," << vec_duration[i] << "\n";
  }

  // Fermer le fichier
  fichierCSV.close();

  std::cout << "res write -> perf_time.csv" << std::endl;
  //---------------------------------------
  // Finalize MPI
  //---------------------------------------
  mfem::Mpi::Finalize();
}
