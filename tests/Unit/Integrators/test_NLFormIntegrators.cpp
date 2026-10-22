/**
 * @file test_NLFormIntegrators.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests of the SLOTH nonlinear integrators: Jacobian (AssembleElementGrad)
 *        against finite differences of the residual (AssembleElementVector)
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
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;

namespace {
using VARS2D = Sloth2D::VARS;

/// @brief What an integrator needs to be built (the arguments of every SLOTH integrator)
struct Context {
  double dt = 0.1;
  std::vector<mfem::ParGridFunction> u_old;  // one per block
  std::vector<mfem::ParGridFunction> aux_old;
  Parameters params{std::vector<Parameter>{}};
  std::vector<VARS2D*> auxvars;
};

/// @brief One integrator to test: name, number of blocks, how to build it
struct Case {
  std::string name;
  int nb_blocks;
  std::function<mfem::BlockNonlinearFormIntegrator*(const Context&)> make;
  int boundary = -1;  ///< boundary index (attribute - 1) for a boundary integrator, -1 otherwise
};

// Coefficients depending on the unknown (implicit): the Jacobian then also
// contains the derivatives of the coefficients, which are checked as well.
Coefficients phase_field_coefficients() {
  const double lambda = 0.15, omega = 120.;
  return Coefficients(Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(omega)),
                      Coefficient(Glossary::Capillary, lambda),
                      Coefficient(Glossary::Mobility, 1.),
                      Coefficient(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda)));
}

template <template <class> class INTEGRATOR>
std::function<mfem::BlockNonlinearFormIntegrator*(const Context&)> maker(
    std::vector<Coefficients> coefficients) {
  return [coefficients](const Context& c) {
    // Through the base class, as OperatorBase::set_nlfi_ptr: init() is public
    // in SlothNLFormIntegrator but redeclared protected in some derived classes
    SlothNLFormIntegrator<VARS2D>* integ = new INTEGRATOR<VARS2D>(
        Geometry::Cartesian, c.dt, c.u_old, c.aux_old, c.params, c.auxvars, coefficients);
    integ->init();
    return integ;
  };
}

/// @brief Builds a boundary integrator (extra arguments: block, boundary index)
template <template <class> class INTEGRATOR>
std::function<mfem::BlockNonlinearFormIntegrator*(const Context&)> bdr_maker(
    std::vector<Coefficients> coefficients, unsigned int block, unsigned int bdr) {
  return [coefficients, block, bdr](const Context& c) {
    SlothNLFormIntegrator<VARS2D>* integ =
        new INTEGRATOR<VARS2D>(Geometry::Cartesian, c.dt, c.u_old, c.aux_old, c.params, c.auxvars,
                               coefficients, block, bdr);
    integ->init();
    return integ;
  };
}

/// @brief Coefficient attached to one boundary (as given by the user in SLOTH)
Coefficient on_boundary(Coefficient c, int bdr) {
  c.set_bdr_index_coef({bdr});
  return c;
}

// Boundary used for the boundary integrators: index 1 = attribute 2 = right
// side of the unit square (x = 1), of length 1
constexpr int kRight = 1;

std::vector<Case> all_cases() {
  const Coefficients none;
  return {
      {"TimeDerivative", 1, maker<TimeNLFormIntegrator>({none})},
      {"SplitTimeDerivative", 2, maker<TimeCHNLFormIntegrator>({none, none})},
      {"AllenCahn", 1, maker<AllenCahnNLFormIntegrator>({phase_field_coefficients()})},
      {"CahnHilliard", 2,
       maker<CahnHilliardNLFormIntegrator>(
           {phase_field_coefficients(), phase_field_coefficients()})},
      {"Fick", 1,
       maker<FickNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::Diffusivity, Scheme::Implicit, H(2.)))})},
      {"Fourier", 1,
       maker<FourierNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::Conductivity, Scheme::Implicit, H(3.)))})},
      // Time-derivative integrators receive k = du/dt (left-hand side of the
      // reduced operator): an implicit coefficient would be evaluated on k,
      // not on u. Cp(T) must be explicit (evaluated on T at the previous step).
      {"HeatTimeDerivative", 1,
       maker<HeatTimeNLFormIntegrator>({Coefficients(
           Coefficient(Glossary::Concentration, 2.),
           Coefficient(Glossary::Cp, Scheme::Explicit, W(5.)))})},
      {"SourceTerm", 1,
       maker<SourceTermNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::SourceTerm, Scheme::Implicit, W(4.)))})},
      {"SplitAllenCahn", 2,
       maker<BlockAllenCahnNLFormIntegrator>(
           {phase_field_coefficients(), phase_field_coefficients()})},
      // Boundary integrators, on the right side
      {"Robin (a(u), b constant)", 1,
       bdr_maker<RobinNLFormIntegrator>(
           {Coefficients(
               on_boundary(Coefficient(Glossary::Robin_a, Scheme::Implicit, H(2.)), kRight),
               on_boundary(Coefficient(Glossary::Robin_b, 3.), kRight))},
           0, kRight),
       kRight},
      {"Neumann (constant flux)", 1,
       bdr_maker<NeumannNLFormIntegrator>(
           {Coefficients(on_boundary(Coefficient(Glossary::Neumann, 2.), kRight))}, 0, kRight),
       kRight},
  };
}

/**
 * @brief Checks the Jacobian of a nonlinear form at x against centred finite
 * differences of its residual, in a random direction v:
 *   J(x) v  ~  (F(x + h v) - F(x - h v)) / 2h.
 * The error of the finite difference is O(h^2): with h = 1e-6 and a tolerance
 * of 1e-5 relative to max |J v|, a wrong term of the Jacobian is detected.
 */
void check_gradient(mfem::ParBlockNonlinearForm& form, const mfem::Vector& x) {
  const double h = 1.e-6;
  mfem::Vector v(x.Size());
  v.Randomize(11);
  mfem::Vector xp(x), xm(x), fp(x.Size()), fm(x.Size()), fd(x.Size()), jv(x.Size());
  xp.Add(h, v);
  xm.Add(-h, v);
  form.Mult(xp, fp);
  form.Mult(xm, fm);
  subtract(0.5 / h, fp, fm, fd);
  form.GetGradient(x).Mult(v, jv);

  // Scale of the problem (all ranks)
  double local_max = jv.Normlinf(), scale = 0.;
  MPI_Allreduce(&local_max, &scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  // Zero Jacobian (e.g. Neumann with a constant flux): the residual must not
  // depend on x either
  const double tol = scale > 0. ? 1.e-5 * scale : 1.e-8;
  for (int r = 0; r < x.Size(); ++r) {
    CAPTURE(r);
    CHECK_THAT(jv(r), WithinAbs(fd(r), tol));
  }
}

double phi_old(const mfem::Vector& p) {
  return 0.5 + 0.2 * std::sin(3. * p[0]) * std::cos(2. * p[1]);
}
}  // namespace

TEST_CASE_METHOD(sloth_unit::UnitSquare, "NLFormIntegrators: Jacobian matches the residual",
                 "[Integrators][derivatives]") {
  // The kernel timers need the MATools timers (see Problems/test_CahnHilliardDirichlet.cpp)
  Profiling::getInstance().disable();

  const Case c = GENERATE(Catch::Generators::from_range(all_cases()));
  INFO(c.name);

  Context ctx;
  for (int i = 0; i < c.nb_blocks; ++i) ctx.u_old.push_back(project(phi_old));

  mfem::Array<mfem::ParFiniteElementSpace*> spaces(c.nb_blocks);
  for (int i = 0; i < c.nb_blocks; ++i) spaces[i] = fes.get();
  mfem::ParBlockNonlinearForm form(spaces);
  mfem::Array<int> marker(pmesh->bdr_attributes.Max());  // must outlive the form
  marker = 0;
  if (c.boundary < 0) {
    form.AddDomainIntegrator(c.make(ctx));
  } else {
    marker[c.boundary] = 1;
    form.AddBoundaryIntegrator(c.make(ctx), marker);
  }

  // Evaluation point: values in [0.2, 0.8] (inside the wells of W, away from 0 and 1)
  mfem::Vector x(c.nb_blocks * fes->GetTrueVSize());
  x.Randomize(5);
  x *= 0.6;
  x += 0.2;
  check_gradient(form, x);
}

/// @brief Sum of the entries of a true-dof vector over all the ranks
double global_sum(const mfem::Vector& v) {
  double local = v.Sum(), global = 0.;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

TEST_CASE_METHOD(sloth_unit::UnitSquare, "NeumannNLFormIntegrator: residual = - integral of g",
                 "[Integrators]") {
  Profiling::getInstance().disable();
  // The shape functions sum to 1: sum_i F_i = - int_right g = - g (length 1)
  const double g = 2.;
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes.get();
  mfem::ParBlockNonlinearForm form(spaces);
  mfem::Array<int> marker(pmesh->bdr_attributes.Max());
  marker = 0;
  marker[kRight] = 1;
  form.AddBoundaryIntegrator(
      bdr_maker<NeumannNLFormIntegrator>(
          {Coefficients(on_boundary(Coefficient(Glossary::Neumann, g), kRight))}, 0, kRight)(ctx),
      marker);

  mfem::Vector x(fes->GetTrueVSize()), y(fes->GetTrueVSize());
  x = 0.5;
  form.Mult(x, y);
  CHECK_THAT(global_sum(y), Catch::Matchers::WithinRel(-g, 1.e-12));
}

// Non-regression: residual (a(u) u - b(u)), Jacobian a + a'(u) u - b'(u).
// The term -b'(u) was missing (wrong Jacobian when b depends on the unknown).
TEST_CASE_METHOD(sloth_unit::UnitSquare, "RobinNLFormIntegrator: Jacobian with b(u)",
                 "[Integrators][derivatives]") {
  Profiling::getInstance().disable();
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes.get();
  mfem::ParBlockNonlinearForm form(spaces);
  mfem::Array<int> marker(pmesh->bdr_attributes.Max());
  marker = 0;
  marker[kRight] = 1;
  form.AddBoundaryIntegrator(
      bdr_maker<RobinNLFormIntegrator>(
          {Coefficients(on_boundary(Coefficient(Glossary::Robin_a, 1.), kRight),
                        on_boundary(Coefficient(Glossary::Robin_b, Scheme::Implicit, W(3.)),
                                    kRight))},
          0, kRight)(ctx),
      marker);

  mfem::Vector x(fes->GetTrueVSize());
  x.Randomize(5);
  x *= 0.6;
  x += 0.2;
  check_gradient(form, x);
}
