/**
 * @file main.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Sintering test taken from Biner book (3D extension)
 * @version 0.1
 * @date 2026-05-07
 *
 * Copyright CEA (c) 2026
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "./SinteringCoefficients.hpp"
#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"

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
  const int DIM = 3;
  using FECollection = Test<DIM>::FECollection;
  using VARS = Test<DIM>::VARS;
  using VAR = Test<DIM>::VAR;
  using PST = Test<DIM>::PST;
  using SPA = Test<DIM>::SPA;
  /////////////////////////
  using OPE = TransientOperator<FECollection, DIM>;
  using PB = Problem<OPE, VARS, PST>;
  // ###########################################
  // ###########################################
  //         Spatial Discretization           //
  // ###########################################
  // ###########################################
  // ##############################
  //           Meshing           //
  // ##############################
  const std::string mesh_type =
      "InlineSquareWithTetraedres";  // type of mesh // "InlineSquareWithTriangles"
  const int order_fe = 1;            // finite element order

  const int refinement_level = 2;  // number of levels of uniform refinement
  const int nx = 50;
  const int ny = 50;
  const int nz = 30;
  const double lx = 50;
  const double ly = 50;
  const double lz = 30;
  const std::tuple<int, int, int, double, double, double>& tuple_of_dimensions = std::make_tuple(
      nx, ny, nz, lx, ly, lz);  // Number of elements and maximum length in each direction

  SPA spatial(mesh_type, order_fe, refinement_level, tuple_of_dimensions);
  // ##############################
  //     Boundary conditions     //
  // ##############################

  auto boundaries = {Boundary("rear", 0, "Neumann"),  Boundary("lower", 1, "Neumann"),
                     Boundary("right", 2, "Neumann"), Boundary("upper", 3, "Neumann"),
                     Boundary("left", 4, "Neumann"),  Boundary("front", 5, "Neumann")};
  // CahnHilliard
  auto ch_bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);
  // allenCahn
  auto ac_bcs = BoundaryConditions<FECollection, DIM>(&spatial, boundaries);

  // ###########################################
  // ###########################################
  //            Physical models               //
  // ###########################################
  // ###########################################
  // ####################
  //     Coefficients    //
  // ####################

  // CahnHilliard
  const double ch_lambda = 5.;
  Coefficient ch_grad_energy(Glossary::GradEnergy, Scheme::Implicit, GradC());
  Coefficient ch_energy(Glossary::FreeEnergy, Scheme::Implicit, Gc());
  Coefficient ch_capillary(Glossary::Capillary, ch_lambda);
  Coefficient ch_mobility(Glossary::Mobility, Scheme::Explicit, D());
  Coefficients ch_coef(ch_energy, ch_capillary, ch_mobility, ch_grad_energy);

  // AllenCahn
  const double ac_mob(10.);
  const double ac_lambda = 2.;
  Coefficient ac_grad_energy(Glossary::GradEnergy, Scheme::Implicit, Gradeta());
  Coefficient ac_energy(Glossary::FreeEnergy, Scheme::Implicit, Geta());
  Coefficient ac_capillary(Glossary::Capillary, ac_lambda);
  Coefficient ac_mobility(Glossary::Mobility, ac_mob);
  Coefficients ac_coef(ac_energy, ac_capillary, ac_mobility, ac_grad_energy);

  // ####################
  //     variables     //
  // ####################

  // CahnHilliard
  auto ch_initial_condition =
      std::function<double(const mfem::Vector&, double)>([](const mfem::Vector& x, double time) {
        double xx1 = x[0] - 14.5;
        double yy1 = x[1] - 25;
        double zz = x[2] - 15;

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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;

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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;

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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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
        double zz = x[2] - 15;
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

        double rr1 = std::sqrt(xx1 * xx1 + yy1 * yy1 + zz * zz);
        double rr2 = std::sqrt(xx2 * xx2 + yy2 * yy2 + zz * zz);
        double rr3 = std::sqrt(xx3 * xx3 + yy3 * yy3 + zz * zz);
        double rr4 = std::sqrt(xx4 * xx4 + yy4 * yy4 + zz * zz);
        double rr5 = std::sqrt(xx5 * xx5 + yy5 * yy5 + zz * zz);
        double rr6 = std::sqrt(xx6 * xx6 + yy6 * yy6 + zz * zz);
        double rr7 = std::sqrt(xx7 * xx7 + yy7 * yy7 + zz * zz);
        double rr8 = std::sqrt(xx8 * xx8 + yy8 * yy8 + zz * zz);
        double rr9 = std::sqrt(xx9 * xx9 + yy9 * yy9 + zz * zz);

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

  // ###########################################
  // ###########################################
  //      Post-processing                     //
  // ###########################################
  // ###########################################

  const std::string& gf_folder_path = "Restart";
  std::vector<int> iterations_list_save_gf = {500, 1000, 1250, 1500, 2500};

  const std::string& main_folder_path = "Saves";
  const int level_of_detail = 1;
  const auto& frequency = 100;

  std::string calculation_path = "CahnHilliard";
  const double threshold = 10.;
  std::map<std::string, std::tuple<double, double>> map_threshold_integral = {{"phi", {-1.1, 1.1}}};
  bool enable_save_specialized_at_iter = true;
  auto ch_p_pst = Parameters(
      Parameter("main_folder_path", main_folder_path), Parameter("gf_folder_path", gf_folder_path),
      Parameter("iterations_list_save_gf", iterations_list_save_gf),
      Parameter("calculation_path", "CahnHilliard"),

      Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
      Parameter("integral_to_compute", map_threshold_integral),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));

  auto ac_p_pst = Parameters(
      Parameter("main_folder_path", main_folder_path), Parameter("calculation_path", "AllenCahn"),
      Parameter("gf_folder_path", gf_folder_path),
      Parameter("iterations_list_save_gf", iterations_list_save_gf),

      Parameter("frequency", frequency), Parameter("level_of_detail", level_of_detail),
      Parameter("enable_save_specialized_at_iter", enable_save_specialized_at_iter));
  // ####################
  //     operators     //
  // ####################

  // CahnHilliard
  std::vector<SPA*> spatials{&spatial, &spatial};
  OPE ch_oper(spatials, {"CahnHilliard"}, TimeScheme::EulerImplicit, "SplitTimeDerivative");
  ch_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));
  const auto& ch_solver = HypreSolverType::HYPRE_GMRES;
  const auto& ch_precond = HyprePreconditionerType::HYPRE_ILU;
  ch_oper.overload_solver(ch_solver);
  ch_oper.overload_preconditioner(ch_precond);

  auto ch_pst = PST(&spatial, ch_p_pst);
  PB ch_pb("CahnHilliard", ch_oper, ch_vars, {ch_coef, ch_coef}, ch_pst, ac_vars);
  //
  //
  //
  std::vector<SPA*> ac_spatials{&spatial, &spatial, &spatial, &spatial, &spatial,
                                &spatial, &spatial, &spatial, &spatial};
  OPE ac_oper(ac_spatials, {"AllenCahn"}, TimeScheme::EulerImplicit, "TimeDerivative");
  ac_oper.overload_nl_solver(
      NLSolverType::NEWTON,
      Parameters(Parameter("description", "Newton solver "), Parameter("print_level", -1),
                 Parameter("rel_tol", 1.e-10), Parameter("abs_tol", 1.e-14)));
  const auto& solver = HypreSolverType::HYPRE_GMRES;
  const auto& precond = HyprePreconditionerType::HYPRE_ILU;
  ac_oper.overload_solver(solver);
  ac_oper.overload_preconditioner(precond);
  auto ac_pst = PST(&spatial, ac_p_pst);
  PB ac_pb("AllenCahn", ac_oper, ac_vars,
           {ac_coef, ac_coef, ac_coef, ac_coef, ac_coef, ac_coef, ac_coef, ac_coef, ac_coef},
           ac_pst, ch_vars);

  // Coupling 1
  auto cc = Coupling("CahnHilliard/AllenCahn Coupling", ac_pb, ch_pb);

  // ###########################################
  // ###########################################
  //            Time-integration              //
  // ###########################################
  // ###########################################
  const double t_initial = 0.0;
  const double t_final = 0.25;
  const double dt = 1.e-4;
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
