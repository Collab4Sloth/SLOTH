/**
 * @file test_Couplings.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for the coupling of several problems (Coupling, TimeDiscretization)
 * @version 0.1
 * @date 2026-10-22
 *
 * Copyright CEA (C) 2026
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <filesystem>
#include <functional>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

namespace {
using namespace Sloth2D;

/// @brief Coefficient equal to the auxiliary field number k (no dependence on u)
class AuxValue final : public FunctionCoefficient {
 public:
  explicit AuxValue(std::size_t k) : k_(k) {}

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [this](const std::span<const double>&, const std::span<const double>&,
                  const std::span<const double>& aux, const unsigned int) { return aux[k_]; };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  GradientF() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&, const unsigned int) {
      return std::vector<double>(u.size(), 0.);
    };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  HessianF() override {
    return [](const std::span<const double>& u, const std::span<const double>&,
              const std::span<const double>&, const unsigned int) {
      return std::vector<double>(u.size() * u.size(), 0.);
    };
  }

 private:
  std::size_t k_;
};

Parameters newton() {
  return Parameters(Parameter("description", std::string("Newton")), Parameter("rel_tol", 1.e-12),
                    Parameter("abs_tol", 1.e-12));
}

/**
 * @brief Two problems on the unit square, solved one after the other at each
 * time step:
 *   P1 (steady):    - div(grad T) = 0, T = 1 on x = 0, T = 2 on x = 1  ->  T = 1 + x
 *   P2 (transient): du/dt = T, u(0) = 0, T auxiliary variable of P2
 * If P1 is solved before P2 at each step, P2 sees T = 1 + x from the first
 * step on, and u(t) = t (1 + x) exactly (implicit Euler, Q1). If P2 were solved
 * first, it would see the initial guess T = 1.5 at the first step.
 */
struct TwoProblems {
  SPA spatial{"InlineSquareWithQuadrangles", 1, 0, std::make_tuple(8, 8, 1., 1.)};
  BCS bcs_T{&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Dirichlet", 2.),
            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Dirichlet", 1.)};
  BCS bcs_u{&spatial, Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
            Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
  VARS T_vars;
  VARS u_vars;
  std::vector<SPA*> spatials{&spatial};
  Coefficients coefs_T{Coefficient(Glossary::Conductivity, 1.)};
  Coefficients coefs_u{Coefficient(Glossary::SourceTerm, Scheme::Implicit, AuxValue(0))};

  static VAR make_T(SPA* s, BCS& b) {
    auto v = VAR(s, b, "T", Glossary::Temperature, 2, 1.5);  // initial guess 1.5
    v.set_additional_information("T");
    return v;
  }

  TwoProblems()
      : T_vars(make_T(&spatial, bcs_T)),
        u_vars(VAR(&spatial, bcs_u, "u", Glossary::Concentration, 2, 0.)) {
    Profiling::getInstance().disable();
  }

  /// @brief L2 distance between u and t (1 + x) (all ranks)
  double u_error(double t) {
    mfem::FunctionCoefficient exact([t](const mfem::Vector& p) { return t * (1. + p[0]); });
    return u_vars[0].get_ref_gf().ComputeL2Error(exact);
  }
};

Parameters time_parameters() {
  return Parameters(Parameter("initial_time", 0.), Parameter("final_time", 0.3),
                    Parameter("time_step", 0.1));
}
}  // namespace

TEST_CASE_METHOD(TwoProblems, "Coupling: two problems in one coupling, solved in order",
                 "[Couplings][Problems]") {
  SteadyOPE oper_T(spatials, {"Fourier"});
  oper_T.overload_nl_solver(NLSolverType::NEWTON, newton());
  TransientOPE oper_u(spatials, {"SourceTerm"}, TimeScheme::EulerImplicit, "TimeDerivative");
  oper_u.overload_nl_solver(NLSolverType::NEWTON, newton());

  SteadyPB p1("T", oper_T, T_vars, {coefs_T});
  TransientPB p2("u", oper_u, u_vars, {coefs_u}, T_vars);  // T: auxiliary variable of P2

  auto coupling = Coupling("T then u", p1, p2);
  auto time = TimeDiscretization(time_parameters(), coupling);
  time.solve();

  CHECK(u_error(0.3) < 1.e-10);
}

TEST_CASE_METHOD(TwoProblems, "Coupling: two couplings, solved in order", "[Couplings]") {
  SteadyOPE oper_T(spatials, {"Fourier"});
  oper_T.overload_nl_solver(NLSolverType::NEWTON, newton());
  TransientOPE oper_u(spatials, {"SourceTerm"}, TimeScheme::EulerImplicit, "TimeDerivative");
  oper_u.overload_nl_solver(NLSolverType::NEWTON, newton());

  SteadyPB p1("T", oper_T, T_vars, {coefs_T});
  TransientPB p2("u", oper_u, u_vars, {coefs_u}, T_vars);
  auto c1 = Coupling("T", p1);
  auto c2 = Coupling("u", p2);
  auto time = TimeDiscretization(time_parameters(), c1, c2);
  time.solve();

  CHECK(u_error(0.3) < 1.e-10);
}

TEST_CASE_METHOD(TwoProblems, "Coupling: unified ParaView output of all the problems",
                 "[Couplings][PostProcessing]") {
  sloth_unit::TemporaryDirectory dir("coupling");
  const Parameters p_T(Parameter("main_folder_path", dir.path()),
                       Parameter("calculation_path", std::string("T")), Parameter("frequency", 1),
                       Parameter("enable_compute_energies", false));
  const Parameters p_u(Parameter("main_folder_path", dir.path()),
                       Parameter("calculation_path", std::string("u")), Parameter("frequency", 1),
                       Parameter("enable_compute_energies", false));
  PST pst_T(&spatial, p_T);
  PST pst_u(&spatial, p_u);

  SteadyOPE oper_T(spatials, {"Fourier"});
  oper_T.overload_nl_solver(NLSolverType::NEWTON, newton());
  TransientOPE oper_u(spatials, {"SourceTerm"}, TimeScheme::EulerImplicit, "TimeDerivative");
  oper_u.overload_nl_solver(NLSolverType::NEWTON, newton());
  SteadyPB p1("T", oper_T, T_vars, {coefs_T}, pst_T);
  TransientPB p2("u", oper_u, u_vars, {coefs_u}, pst_u, T_vars);

  auto coupling = Coupling("T then u", p1, p2);
  auto time = TimeDiscretization(time_parameters() + Parameter("vtk_unified", true), coupling);
  time.solve();

  CHECK(u_error(0.3) < 1.e-10);
  MPI_Barrier(MPI_COMM_WORLD);
  if (mfem::Mpi::Root()) {
    bool pvd_found = false;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(dir.path())) {
      pvd_found = pvd_found || entry.path().extension() == ".pvd";
    }
    CHECK(pvd_found);
  }
}

// -----------------------------------------------------------------------------
// Stationarity criteria: the time loop stops as soon as every criterion given
// is met. Problems without criterion are not considered. Here the criterion is
// on T only: T is the same at steps 1 and 2, so the loop stops at t = 0.2 while
// u still evolves (u = 0.2 (1 + x) instead of 0.3 (1 + x)).
// -----------------------------------------------------------------------------
TEST_CASE_METHOD(TwoProblems, "Coupling: the time loop stops when the criteria given are met",
                 "[Couplings][Convergence]") {
  SteadyOPE oper_T(spatials, {"Fourier"});
  oper_T.overload_nl_solver(NLSolverType::NEWTON, newton());
  TransientOPE oper_u(spatials, {"SourceTerm"}, TimeScheme::EulerImplicit, "TimeDerivative");
  oper_u.overload_nl_solver(NLSolverType::NEWTON, newton());

  auto cvg = Convergence(PhysicalConvergence(ConvergenceType::ABSOLUTE_MAX, 1.e-12));
  SteadyPB p1("T", oper_T, T_vars, {coefs_T}, cvg);
  TransientPB p2("u", oper_u, u_vars, {coefs_u}, T_vars);
  auto coupling = Coupling("T then u", p1, p2);
  auto time = TimeDiscretization(time_parameters(), coupling);
  time.solve();

  CHECK(u_error(0.2) < 1.e-10);
}
