/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Sintering test taken from "Programming Phase-field modeling - Biner 2017"
 * @version 0.1
 * @date 2026-05-07
 *
 * Copyright CEA (c) 2026
 *
 */

#include "./SinteringCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  setVerbosity(Verbosity::Verbose);   // output level
  mfem::Mpi::Init(argc, argv);        // MPI
  mfem::Hypre::Init();                // Hypre solvers
  Profiling::getInstance().enable();  // Profiling
  using namespace Sloth2D;            // dim., aliases

  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  const int nx = 200, ny = 200;
  const double lx = 50, ly = 50;
  const int fe_order = 1;
  const int level_uniform_refinement = 0;
  SPA spatial("InlineSquareWithQuadrangles", fe_order, level_uniform_refinement,
              std::make_tuple(nx, ny, lx, ly));
  auto boundaries = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                     Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};

  // CahnHilliard
  auto ch_bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);
  // allenCahn
  auto ac_bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ####################
  //     Coefficients    //
  // ####################

  // CahnHilliard
  Coefficient ch_grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradC());
  Coefficient ch_energy(Glossary::FreeEnergy, Scheme::Implicit, Gc());
  Coefficient ch_capillary(Glossary::Capillary, 5.);
  Coefficient ch_mobility(Glossary::Mobility, Scheme::Explicit, D());
  Coefficients ch_coeffs(ch_energy, ch_capillary, ch_mobility, ch_grad_energy);

  // AllenCahn
  Coefficient ac_grad_energy(Glossary::GradEnergy, Scheme::Implicit, Gradeta());
  Coefficient ac_energy(Glossary::FreeEnergy, Scheme::Implicit, Geta());
  Coefficient ac_capillary(Glossary::Capillary, 2.0);
  Coefficient ac_mobility(Glossary::Mobility, 10.0);
  Coefficients ac_coeffs(ac_energy, ac_capillary, ac_mobility, ac_grad_energy);

  // ####################
  //     variables     //
  // ####################

  // CahnHilliard
  auto ch_initial_condition =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;

        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;

        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;

        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;

        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;

        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 5;
        if (rr1 < R || rr2 < R || rr3 < R || rr4 < R || rr5 < R || rr6 < 0.5 * R || rr7 < 0.5 * R ||
            rr8 < 0.5 * R || rr9 < 0.5 * R) {
          sol = 1;
        }

        return sol;
      });
  auto ch_initial_condition_mu =
      std::function<double(const mfem::Vector&, double)>([&](const mfem::Vector& x, double time) {
        double phi = ch_initial_condition(x, time);
        double sol = 2 * phi * (1 - phi) * (1 - 2 * phi);

        return sol;
      });

  auto ch_initial_condition_1 = AnalyticalFunctions<DIM>(ch_initial_condition);
  double mu_initial_condition = 0.0;
  auto ch_v1 = VAR(&spatial, ch_bcs, "c", Glossary::PhaseField, 2, ch_initial_condition);
  ch_v1.set_additional_information("c");
  auto ch_v2 = VAR(&spatial, ch_bcs, "mu", Glossary::ChemicalPotential, 2, ch_initial_condition_mu);
  ch_v2.set_additional_information("mu");
  auto ch_vars = VARS(ch_v1, ch_v2);

  // AllenCahn

  auto ac_initial_condition_1 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;

        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;

        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 5;
        if (rr1 < R) {
          sol = 1;
        }

        return sol;
      });
  auto ac_initial_condition_2 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 5;
        if (rr2 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_3 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;

        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 5;
        if (rr3 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_4 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 5;
        if (rr4 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_5 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 5;
        if (rr5 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_6 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 2.5;
        if (rr6 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_7 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 2.5;
        if (rr7 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_8 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 2.5;
        if (rr8 < R) {
          sol = 1;
        }
        return sol;
      });
  auto ac_initial_condition_9 =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double xx2 = x[0] - 35.5;
        double yy2 = x[1] - 25;
        double xx3 = x[0] - 25;
        double yy3 = x[1] - 14.5;
        double xx4 = x[0] - 25;
        double yy4 = x[1] - 35.5;

        double xx5 = x[0] - 25;
        double yy5 = x[1] - 25;

        double xx6 = x[0] - 19.5;
        double yy6 = x[1] - 19.5;
        double xx7 = x[0] - 30.5;
        double yy7 = x[1] - 19.5;
        double xx8 = x[0] - 19.5;
        double yy8 = x[1] - 30.5;
        double xx9 = x[0] - 30.5;
        double yy9 = x[1] - 30.5;

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9);

        double sol = 0.0;
        double R = 2.5;
        if (rr9 < R) {
          sol = 1;
        }
        return sol;
      });

  auto ac_v1 = VAR(&spatial, ac_bcs, "eta_1", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_1));
  ac_v1.set_additional_information("eta");
  auto ac_v2 = VAR(&spatial, ac_bcs, "eta_2", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_2));
  ac_v2.set_additional_information("eta");
  auto ac_v3 = VAR(&spatial, ac_bcs, "eta_3", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_3));
  ac_v3.set_additional_information("eta");
  auto ac_v4 = VAR(&spatial, ac_bcs, "eta_4", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_4));
  ac_v4.set_additional_information("eta");
  auto ac_v5 = VAR(&spatial, ac_bcs, "eta_5", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_5));
  ac_v5.set_additional_information("eta");
  auto ac_v6 = VAR(&spatial, ac_bcs, "eta_6", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_6));
  ac_v6.set_additional_information("eta");
  auto ac_v7 = VAR(&spatial, ac_bcs, "eta_7", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_7));
  ac_v7.set_additional_information("eta");
  auto ac_v8 = VAR(&spatial, ac_bcs, "eta_8", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_8));
  ac_v8.set_additional_information("eta");
  auto ac_v9 = VAR(&spatial, ac_bcs, "eta_9", Glossary::PhaseField, 2,
                   AnalyticalFunctions<DIM>(ac_initial_condition_9));
  ac_v9.set_additional_information("eta");
  auto ac_vars = VARS(ac_v1, ac_v2, ac_v3, ac_v4, ac_v5, ac_v6, ac_v7, ac_v8, ac_v9);


  // ###################################################
  //   Two physics = two operators + integrators      //
  //   TransientOPE + time scheme (eg. EulerImplicit) //
  //   F(u) : SplitTimeDerivative  , TimeDerivative   //
  //   G(U) : CahnHilliard         ; AllenCahn        //
  // ###################################################

  // CahnHilliard
  std::vector<SPA*> spatials{&spatial, &spatial};  // One FE space per variable (here, the same)
  TransientOPE ch_oper(spatials, {"CahnHilliard"}, TimeScheme::EulerImplicit,
                       "SplitTimeDerivative");
  ch_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));
  // AllenCahn
  std::vector<SPA*> ac_spatials{
      &spatial, &spatial, &spatial, &spatial, &spatial,
      &spatial, &spatial, &spatial, &spatial};  // One FE space per variable (here, the same)
  TransientOPE ac_oper(ac_spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  ac_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", 1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));

  // ###################################################
  // Post-Processing
  // VTK -> paraview
  // CSV -> python script
  // ###################################################
  // location
  const std::string& main_folder_path = "Saves";
  std::string calculation_path = "CahnHilliard_AllenCahn";
  // VTK
  const int vtk_frequency = 1;

  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {{"phi", {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  int specialized_frequency = 1;
  auto param_unified_pst = Parameters(
      Parameter("main_folder_path", main_folder_path),
      Parameter("calculation_path", calculation_path), Parameter("frequency", vtk_frequency),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter),
      Parameter("specialized_frequency", specialized_frequency),
      Parameter("enable_compute_energies", true),
      Parameter("integral_to_compute", map_threshold_integral));
  bool vtk_unified = false;

  // Post-processing object (unified) associated with a mesh
  auto unified_pst = PST(&spatial, param_unified_pst);

  // ###################################################
  // Problems
  // ###################################################
  // A vector of coefficients must be defined : one per variable, here the same
  // ch_vars = primary variables
  // ac_vars = auxiliary variables
  TransientPB ch_pb("Cahn-Hilliard", ch_oper, ch_vars, {ch_coeffs, ch_coeffs}, unified_pst,
                    ac_vars);

  // ac_vars = primary variables
  // ch_vars = auxiliary variables
  TransientPB ac_pb("Allen-Cahn", ac_oper, ac_vars,
                    {ac_coeffs, ac_coeffs, ac_coeffs, ac_coeffs, ac_coeffs, ac_coeffs, ac_coeffs,
                     ac_coeffs, ac_coeffs},
                    unified_pst, ch_vars);

  // ###################################################
  // Coupling : here  composed of two problem
  // ###################################################
  auto cc = Coupling("Main coupling ", ac_pb, ch_pb);

  // ###########################################
  // ###########################################
  //            Resolution                    //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 2.e-4; 2. ; //100.;
  const double dt = 1.e-4;
  auto time_params = Parameters(
      Parameter("initial_time", t_initial), Parameter("final_time", t_final),
      Parameter("time_step", dt), Parameter("adaptive_time_step", false),
      Parameter("time_step_min", 1.e-4), Parameter("time_step_max", 20.),
      Parameter("time_step_growth_factor", 1.5), Parameter("time_step_cutback_factor", 0.5),
      Parameter("newton_iterations_low", 4), Parameter("newton_iterations_high", 10),
      Parameter("max_time_step_rejections", 20));
  auto time = TimeDiscretization(time_params, cc);

  time.solve();
  //---------------------------------------
  // footer
  //---------------------------------------
  Profiling::getInstance().print();
  MPI_Finalize();
  return 0;
}
