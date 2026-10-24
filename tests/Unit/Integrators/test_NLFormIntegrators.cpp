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
#include <span>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

// Mobility laws generated at configure time from Integrators/MobilityCoefficients.json:
//   MobilityPhi   = m0 (1 + phi^2)          (one variable)
//   MobilityCross = m0 phi0 (1 - phi1)      (two phases, monolithic coupling)
//   MobilityDegenerate = m0 phi0 (1 - phi0) (two phases, depends on phi0 only)
#include "MobilityCross.hpp"
#include "MobilityDegenerate.hpp"
#include "MobilityPhi.hpp"

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
Coefficients phase_field_coefficients(Coefficient mobility = Coefficient(Glossary::Mobility, 1.)) {
  const double lambda = 0.15, omega = 120.;
  return Coefficients(Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(omega)),
                      Coefficient(Glossary::Capillary, lambda), mobility,
                      Coefficient(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda)));
}

/// @brief Allen-Cahn coefficients with a mobility defined on the two phases
template <class MOBILITY>
Coefficients two_phase_mobility_coefficients() {
  const double lambda = 0.15;
  return Coefficients(Coefficient(Glossary::FreeEnergy, 0.),
                      Coefficient(Glossary::Capillary, lambda),
                      Coefficient(Glossary::Mobility, Scheme::Implicit, MOBILITY()),
                      Coefficient(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda)));
}

template <template <class> class INTEGRATOR>
std::function<mfem::BlockNonlinearFormIntegrator*(const Context&)> maker(
    std::vector<Coefficients> coefficients) {
  return [coefficients](const Context& c) {
    // Through the base class, as OperatorBase::set_nlfi_ptr
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

// Boundaries used for the boundary integrators (index = attribute - 1) of the
// unit square: right side x = 1, left side x = 0, both of length 1
constexpr int kRight = 1;
constexpr int kLeft = 3;

/// @brief Constant advection velocity v = (vx, vy)
Coefficients velocity(double vx, double vy) {
  return Coefficients(Coefficient(Glossary::Velocity_X, vx), Coefficient(Glossary::Velocity_Y, vy));
}

/// @brief Constant velocity and an inflow value on one boundary
Coefficients velocity_with_inflow(double vx, double vy, double inflow, int bdr) {
  Coefficients c = velocity(vx, vy);
  c.add(on_boundary(Coefficient(Glossary::Inflow, inflow), bdr));
  return c;
}

std::vector<Case> all_cases() {
  const Coefficients none;
  return {
      {"TimeDerivative", 1, maker<TimeNLFormIntegrator>({none})},
      {"SplitTimeDerivative", 2, maker<TimeCHNLFormIntegrator>({none, none})},
      {"AllenCahn", 1, maker<AllenCahnNLFormIntegrator>({phase_field_coefficients()})},
      {"AllenCahn, mobility M(phi)", 1,
       maker<AllenCahnNLFormIntegrator>({phase_field_coefficients(
           Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityPhi()))})},
      // Two phases solved together, mobility depending on both: the Jacobian
      // gets dM/dphi_j in every block (b, j). Constant free energy, so that
      // only the mobility couples the blocks.
      {"AllenCahn, 2 phases, mobility M(phi0, phi1)", 2,
       maker<AllenCahnNLFormIntegrator>({two_phase_mobility_coefficients<MobilityCross>(),
                                         two_phase_mobility_coefficients<MobilityCross>()})},
      // Defined on the two phases but depending on phi0 only: dM/dphi1 = 0
      {"AllenCahn, 2 phases, mobility M(phi0) = phi0 (1 - phi0)", 2,
       maker<AllenCahnNLFormIntegrator>({two_phase_mobility_coefficients<MobilityDegenerate>(),
                                         two_phase_mobility_coefficients<MobilityDegenerate>()})},
      {"CahnHilliard", 2,
       maker<CahnHilliardNLFormIntegrator>(
           {phase_field_coefficients(), phase_field_coefficients()})},
      // Mobility depending on phi (block 1 = equation of mu), e.g. a degenerate mobility
      {"CahnHilliard, mobility M(phi)", 2,
       maker<CahnHilliardNLFormIntegrator>(
           {phase_field_coefficients(),
            phase_field_coefficients(
                Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityPhi()))})},
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
      // Advection with a constant velocity v = (1, 0.5)
      {"Advection", 1, maker<AdvectionNLFormIntegrator>({velocity(1., 0.5)})},
      {"SplitAdvection", 2, maker<SplitAdvectionNLFormIntegrator>({velocity(1., 0.5), none})},
      {"AdvectionBoundary, outflow (v.n > 0)", 1,
       bdr_maker<AdvectionBoundaryNLFormIntegrator>(
           {velocity_with_inflow(1., 0.5, 3., kRight)}, 0, kRight),
       kRight},
      {"AdvectionBoundary, inflow (v.n < 0)", 1,
       bdr_maker<AdvectionBoundaryNLFormIntegrator>({velocity_with_inflow(1., 0.5, 3., kLeft)},
                                                    0, kLeft),
       kLeft},
      {"AdvectionBoundaryCH, outflow", 2,
       bdr_maker<AdvectionBoundaryCHNLFormIntegrator>(
           {velocity_with_inflow(1., 0.5, 3., kRight), none}, 0, kRight),
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

// Residual (a(u) u - b(u)), Jacobian a + a'(u) u - b'(u), with a and b
// depending on the unknown.
TEST_CASE_METHOD(sloth_unit::UnitSquare, "RobinNLFormIntegrator: Jacobian with b(u)",
                 "[Integrators][derivatives]") {
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

// =============================================================================
// Integrators with auxiliary variables (melting, latent heat)
// =============================================================================
// The auxiliary variables are SLOTH Variables, recognised by their additional
// information ("T", "phi"). Their grid functions must live on the same mesh as
// the tested form: everything is built on the mesh of a SpatialDiscretization.
namespace {
using SPA2D = Sloth2D::SPA;
using BCS2D = Sloth2D::BCS;
using VAR2D = Sloth2D::VAR;

struct SlothSquare {
  std::unique_ptr<SPA2D> spatial;
  std::unique_ptr<BCS2D> bcs;  // Neumann everywhere

  SlothSquare() {
    spatial = std::make_unique<SPA2D>("InlineSquareWithQuadrangles", 1, 0,
                                      std::make_tuple(4, 4, 1., 1.));
    auto neumann = {Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
                    Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann")};
    bcs = std::make_unique<BCS2D>(spatial.get(), neumann);
  }
  mfem::ParFiniteElementSpace* fes() const { return spatial->get_finite_element_space(); }

  mfem::ParGridFunction project(const std::function<double(const mfem::Vector&)>& f) const {
    mfem::ParGridFunction gf(fes());
    mfem::FunctionCoefficient c(f);
    gf.ProjectCoefficient(c);
    return gf;
  }

  /// @brief Auxiliary variable with its additional information (e.g. "T")
  VARS2D aux_variable(const std::string& name, const std::string& info, GlossaryQuantity qty,
                      std::function<double(const mfem::Vector&, double)> f) const {
    auto v = VAR2D(spatial.get(), *bcs, name, qty, 2, AnalyticalFunctions<DIM2>(f));
    v.set_additional_information(info);
    return VARS2D(v);
  }
  static constexpr int DIM2 = 2;

  /// @brief Auxiliary variable with several additional informations (e.g. {"U", "mu"})
  VARS2D aux_variable(const std::string& name, const std::vector<std::string>& infos,
                      GlossaryQuantity qty,
                      std::function<double(const mfem::Vector&, double)> f) const {
    auto v = VAR2D(spatial.get(), *bcs, name, qty, 2, AnalyticalFunctions<DIM2>(f));
    for (const auto& info : infos) v.set_additional_information(info);
    return VARS2D(v);
  }

  /// @brief Random point in [0.2, 0.8] (one block)
  mfem::Vector random_point() const {
    mfem::Vector x(fes()->GetTrueVSize());
    x.Randomize(5);
    x *= 0.6;
    x += 0.2;
    return x;
  }
};

/// @brief One-block form with one domain integrator, on the given space
std::unique_ptr<mfem::ParBlockNonlinearForm> one_block_form(
    mfem::ParFiniteElementSpace* fes, mfem::BlockNonlinearFormIntegrator* i) {
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes;
  auto form = std::make_unique<mfem::ParBlockNonlinearForm>(spaces);
  form->AddDomainIntegrator(i);
  return form;
}

Coefficients melting_coefficients(Coefficient mobility) {
  return Coefficients(mobility,
                      Coefficient(Glossary::InterpolationFunction, Scheme::Implicit, H()));
}
}  // namespace

TEST_CASE_METHOD(SlothSquare, "MeltingConstantNLFormIntegrator: Jacobian",
                 "[Integrators][derivatives]") {
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  ctx.params = Parameters(Parameter("melting_factor", 2.));
  auto form = one_block_form(
      fes(), maker<MeltingConstantNLFormIntegrator>(
                 {melting_coefficients(Coefficient(Glossary::Mobility, 3.))})(ctx));
  check_gradient(*form, random_point());
}

TEST_CASE_METHOD(SlothSquare, "MeltingTemperatureNLFormIntegrator: Jacobian",
                 "[Integrators][derivatives]") {
  // T = 1000 + 200 x: melting (T > 1100) on half of the domain
  auto temperature =
      aux_variable("T", "T", Glossary::Temperature,
                   [](const mfem::Vector& p, double) { return 1000. + 200. * p[0]; });
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  ctx.auxvars = {&temperature};
  ctx.aux_old = {temperature.getVariables()[0].get_gf()};
  ctx.params = Parameters(Parameter("melting_temperature", 1100.),
                          Parameter("melting_enthalpy", 5.));
  auto form = one_block_form(
      fes(), maker<MeltingTemperatureNLFormIntegrator>(
                 {melting_coefficients(Coefficient(Glossary::Mobility, 3.))})(ctx));
  check_gradient(*form, random_point());
}

TEST_CASE_METHOD(SlothSquare, "LatentHeatNLFormIntegrator: residual and Jacobian",
                 "[Integrators][derivatives]") {
  // phi - phi_old = c everywhere: residual = - ((c / dt)^2 / M) psi
  const double c = 0.3, dt = 0.1, M = 2.;
  auto phi = aux_variable("phi", "phi", Glossary::PhaseField,
                          [](const mfem::Vector& p, double) { return 0.4 + 0.1 * p[0]; });
  Context ctx;
  ctx.dt = dt;
  ctx.u_old.push_back(project(phi_old));
  ctx.auxvars = {&phi};
  ctx.aux_old = {project([c](const mfem::Vector& p) { return 0.4 + 0.1 * p[0] - c; })};
  auto form = one_block_form(
      fes(), maker<LatentHeatNLFormIntegrator>(
                 {Coefficients(Coefficient(Glossary::Mobility, M))})(ctx));

  const mfem::Vector x = random_point();
  SECTION("residual: sum = - (c / dt)^2 / M over the unit square") {
    mfem::Vector y(x.Size());
    form->Mult(x, y);
    CHECK_THAT(global_sum(y), Catch::Matchers::WithinRel(-(c / dt) * (c / dt) / M, 1.e-10));
  }
  SECTION("Jacobian: zero (the residual does not depend on the unknown)") {
    check_gradient(*form, x);
  }
}

// Residual M(u) alpha h'(phi), Jacobian M alpha h''(phi) + M'(u) alpha h'(phi).
// The second term comes from the mobility depending on u.
TEST_CASE_METHOD(SlothSquare, "MeltingConstantNLFormIntegrator: Jacobian with M(u)",
                 "[Integrators][derivatives]") {
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  auto form = one_block_form(
      fes(), maker<MeltingConstantNLFormIntegrator>({melting_coefficients(
                 Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityPhi()))})(ctx));
  check_gradient(*form, random_point());
}

// =============================================================================
// Mass diffusion flux
//   residual: (sum_i M_i grad mu_i + D (grad u - grad u_old), grad psi)
// Physical flux M_i grad mu_i (mu_i, M_i: auxiliary variables [element, "mu"],
// [element, "inter_mob"]) and a stabilisation term D grad(u - u_old), which
// vanishes when u no longer changes over the time step.
// =============================================================================
namespace {
struct MassFluxSetup : SlothSquare {
  VARS2D mu = aux_variable("mu_U", std::vector<std::string>{"U", "mu"},
                           Glossary::ChemicalPotential,
                           [](const mfem::Vector& p, double) { return p[0]; });  // mu = x
  VARS2D mob = aux_variable("M_U", std::vector<std::string>{"U", "inter_mob"},
                            Glossary::InterDiffusionMobility,
                            [](const mfem::Vector&, double) { return 2.; });  // M = 2

  std::unique_ptr<mfem::ParBlockNonlinearForm> form(Coefficient diffusivity, Context& ctx) {
    ctx.u_old.push_back(project(phi_old));
    ctx.auxvars = {&mu, &mob};
    ctx.aux_old = {mu.getVariables()[0].get_gf(), mob.getVariables()[0].get_gf()};
    return one_block_form(fes(), maker<MassDiffusionFluxNLFormIntegrator>(
                                     {Coefficients(diffusivity)})(ctx));
  }
};
}  // namespace

TEST_CASE_METHOD(MassFluxSetup, "MassDiffusionFluxNLFormIntegrator: residual and Jacobian",
                 "[Integrators][derivatives]") {
  Context ctx;
  auto f = form(Coefficient(Glossary::Diffusivity, 1.), ctx);

  SECTION("Jacobian") { check_gradient(*f, random_point()); }

  SECTION("residual at u = u_old: no stabilisation, physical flux M grad(mu) = (2, 0)") {
    // sum_i X_i F_i = (flux, grad x) = 2 over the unit square, X = interpolant of x
    mfem::Vector u_old, y(fes()->GetTrueVSize()), X;
    ctx.u_old[0].GetTrueDofs(u_old);
    f->Mult(u_old, y);
    project([](const mfem::Vector& p) { return p[0]; }).GetTrueDofs(X);
    CHECK_THAT(mfem::InnerProduct(MPI_COMM_WORLD, X, y), Catch::Matchers::WithinRel(2., 1.e-10));
  }
}

TEST_CASE_METHOD(SlothSquare, "MassDiffusionFluxNLFormIntegrator: invalid data",
                 "[Integrators][errors]") {
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  SECTION("no chemical potential") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        maker<MassDiffusionFluxNLFormIntegrator>(
            {Coefficients(Coefficient(Glossary::Diffusivity, 1.))})(ctx),
        "Neither chemical potentials");
  }
  SECTION("chemical potential without element") {
    auto mu = aux_variable("mu", std::vector<std::string>{"mu"}, Glossary::ChemicalPotential,
                           [](const mfem::Vector&, double) { return 0.; });
    ctx.auxvars = {&mu};
    ctx.aux_old = {mu.getVariables()[0].get_gf()};
    SLOTH_CHECK_MFEM_ERROR_WITH(
        maker<MassDiffusionFluxNLFormIntegrator>(
            {Coefficients(Coefficient(Glossary::Diffusivity, 1.))})(ctx),
        "Expected [element, 'mu']");
  }
}

// Physical flux scaled by the temperature T (auxiliary variable ["T"]):
//   ScaleVariablesByTemperature (diffusion chemical potentials dmu, [element, "dmu"]):
//     M grad(dmu / T) = M (grad dmu / T - dmu grad T / T^2)
//   ScaleCoefficientsByTemperature: (M / (R T)) grad mu
// Moments (F, x) at u = u_old (no stabilisation), fields exact in Q1.
TEST_CASE_METHOD(SlothSquare, "MassDiffusionFluxNLFormIntegrator: scaling by the temperature",
                 "[Integrators]") {
  using Field = std::function<double(const mfem::Vector&, double)>;
  auto make = [this](const std::string& name, const std::vector<std::string>& infos,
                     GlossaryQuantity qty, Field f) { return aux_variable(name, infos, qty, f); };
  auto mob = make("M_U", {"U", "inter_mob"}, Glossary::InterDiffusionMobility,
                  [](const mfem::Vector&, double) { return 2.; });
  auto T_uniform = make("T", {"T"}, Glossary::Temperature,
                        [](const mfem::Vector&, double) { return 2.; });
  auto T_affine = make("T", {"T"}, Glossary::Temperature,
                       [](const mfem::Vector& p, double) { return 1. + p[0]; });
  auto dmu_x = make("dmu_U", {"U", "dmu"}, Glossary::ChemicalPotential,
                    [](const mfem::Vector& p, double) { return p[0]; });
  auto dmu_3T = make("dmu_U", {"U", "dmu"}, Glossary::ChemicalPotential,
                     [](const mfem::Vector& p, double) { return 3. * (1. + p[0]); });
  auto mu_x = make("mu_U", {"U", "mu"}, Glossary::ChemicalPotential,
                   [](const mfem::Vector& p, double) { return p[0]; });

  const Parameters scale_variables(Parameter("EnableDiffusionChemicalPotentials", true),
                                   Parameter("ScaleVariablesByTemperature", true));

  // (F, x) at u = u_old, X = interpolant of x
  auto moment = [this](const std::vector<VARS2D*>& aux, const Parameters& params) {
    Context ctx;
    ctx.u_old.push_back(project(phi_old));
    ctx.auxvars = aux;
    for (auto* v : aux) ctx.aux_old.push_back(v->getVariables()[0].get_gf());
    ctx.params = params;
    auto f = one_block_form(fes(), maker<MassDiffusionFluxNLFormIntegrator>(
                                       {Coefficients(Coefficient(Glossary::Diffusivity, 1.))})(ctx));
    mfem::Vector u_old, y(fes()->GetTrueVSize()), X;
    ctx.u_old[0].GetTrueDofs(u_old);
    f->Mult(u_old, y);
    project([](const mfem::Vector& p) { return p[0]; }).GetTrueDofs(X);
    return mfem::InnerProduct(MPI_COMM_WORLD, X, y);
  };

  SECTION("dmu = x, T = 2: flux M grad(x) / 2 = (1, 0)") {
    const std::vector<VARS2D*> aux{&dmu_x, &mob, &T_uniform};
    CHECK_THAT(moment(aux, scale_variables),
               Catch::Matchers::WithinRel(1., 1.e-10));
  }
  SECTION("dmu = 3 T, T = 1 + x: dmu / T is uniform, no flux") {
    const std::vector<VARS2D*> aux{&dmu_3T, &mob, &T_affine};
    CHECK_THAT(moment(aux, scale_variables),
               Catch::Matchers::WithinAbs(0., 1.e-12));
  }
  SECTION("mu = x, T = 2, mobility scaled: flux (M / (R T)) grad(x) = (1 / R, 0)") {
    const Parameters scale_coefficients(Parameter("ScaleCoefficientsByTemperature", true));
    const std::vector<VARS2D*> aux{&mu_x, &mob, &T_uniform};
    CHECK_THAT(moment(aux, scale_coefficients),
               Catch::Matchers::WithinRel(1. / Physical::R, 1.e-10));
  }
  SECTION("the temperature is required to scale") {
    const std::vector<VARS2D*> no_temperature{&dmu_x, &mob};
    SLOTH_CHECK_MFEM_ERROR_WITH(moment(no_temperature, scale_variables), "No temperature found");
  }
}

// Stabilisation coefficient depending on u (not the usual case, D is generally
// constant): Jacobian D grad psi_j . grad psi_i + D'(u) psi_j grad(u - u_old) . grad psi_i
TEST_CASE_METHOD(MassFluxSetup, "MassDiffusionFluxNLFormIntegrator: Jacobian with D(u)",
                 "[Integrators][derivatives]") {
  Context ctx;
  auto f = form(Coefficient(Glossary::Diffusivity, Scheme::Implicit, H(3.)), ctx);
  check_gradient(*f, random_point());
}

// =============================================================================
// Advection in divergence form: div(v u), v not necessarily solenoidal
//   domain:   - (u v, grad psi)
//   boundary: + (u v.n, psi) where v.n > 0, (g v.n, psi) where v.n < 0 (inflow g)
// =============================================================================
namespace {
/// @brief Sum_i X_i F_i over all the ranks, X = interpolant of x: (F, x)
double moment_x(const SlothSquare& sq, const mfem::Vector& F) {
  mfem::Vector X;
  sq.project([](const mfem::Vector& p) { return p[0]; }).GetTrueDofs(X);
  return mfem::InnerProduct(MPI_COMM_WORLD, X, F);
}
}  // namespace

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

TEST_CASE_METHOD(SlothSquare,
                 "AdvectionNLFormIntegrator: divergence form with a non-solenoidal velocity",
                 "[Integrators]") {
  // v = (b, 0) with the auxiliary field b = x (exact in Q1): v = (x, 0), div v = 1.
  // With u = 1: div(v u) = 1, so (div(v u), x) = int x = 1/2 over the unit square.
  auto b = aux_variable("b", "b", Glossary::Default,
                        [](const mfem::Vector& p, double) { return p[0]; });
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  ctx.auxvars = {&b};
  ctx.aux_old = {b.getVariables()[0].get_gf()};
  auto v = [](int bdr) {
    Coefficients c(Coefficient(Glossary::Velocity_X, Scheme::Implicit, AuxValue(0)),
                   Coefficient(Glossary::Velocity_Y, 0.));
    if (bdr >= 0) c.add(on_boundary(Coefficient(Glossary::Inflow, 0.), bdr));
    return c;
  };

  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes();
  mfem::ParBlockNonlinearForm form(spaces);
  form.AddDomainIntegrator(maker<AdvectionNLFormIntegrator>({v(-1)})(ctx));

  mfem::Vector u(fes()->GetTrueVSize()), F(fes()->GetTrueVSize());
  u = 1.;

  SECTION("domain term alone: - (v, grad x) = - int x = -1/2") {
    form.Mult(u, F);
    CHECK_THAT(moment_x(*this, F), Catch::Matchers::WithinAbs(-0.5, 1.e-10));
  }
  SECTION("with the outflow boundary x = 1: (div(v u), x) = 1/2") {
    // Outflow on x = 1 (v.n = 1); v = 0 on x = 0, v.n = 0 on y = 0 and y = 1
    mfem::Array<int> marker(spatial->get_mesh()->bdr_attributes.Max());
    marker = 0;
    marker[kRight] = 1;
    form.AddBoundaryIntegrator(bdr_maker<AdvectionBoundaryNLFormIntegrator>({v(kRight)}, 0,
                                                                             kRight)(ctx),
                               marker);
    form.Mult(u, F);
    CHECK_THAT(moment_x(*this, F), Catch::Matchers::WithinAbs(0.5, 1.e-10));
  }
}

TEST_CASE_METHOD(SlothSquare, "AdvectionNLFormIntegrator: gradient velocity v = alpha grad(a)",
                 "[Integrators]") {
  // a = x^2 / 2: in Q1, grad(a_h) is constant per element (mean slope), and
  // - (grad(a_h), grad x) = - int d(a_h)/dx = -(a(1) - a(0)) = -1/2 exactly (u = 1)
  auto a = aux_variable("a", "a", Glossary::Default,
                        [](const mfem::Vector& p, double) { return 0.5 * p[0] * p[0]; });
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  ctx.auxvars = {&a};
  ctx.aux_old = {a.getVariables()[0].get_gf()};
  ctx.params = Parameters(Parameter("velocity_gradient", std::string("a")));
  auto form = one_block_form(fes(), maker<AdvectionNLFormIntegrator>(
                                        {Coefficients(Coefficient(Glossary::VelocityPrefactor,
                                                                  1.))})(ctx));
  mfem::Vector u(fes()->GetTrueVSize()), F(fes()->GetTrueVSize());
  u = 1.;
  form->Mult(u, F);
  CHECK_THAT(moment_x(*this, F), Catch::Matchers::WithinAbs(-0.5, 1.e-10));
}

TEST_CASE_METHOD(SlothSquare, "AdvectionBoundaryNLFormIntegrator: inflow value",
                 "[Integrators]") {
  // v = (1, 0), inflow g = 3 on x = 0 (v.n = -1): residual (g v.n, psi), sum = -3
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes();
  mfem::ParBlockNonlinearForm form(spaces);
  mfem::Array<int> marker(spatial->get_mesh()->bdr_attributes.Max());
  marker = 0;
  marker[kLeft] = 1;
  form.AddBoundaryIntegrator(bdr_maker<AdvectionBoundaryNLFormIntegrator>(
                                 {velocity_with_inflow(1., 0., 3., kLeft)}, 0, kLeft)(ctx),
                             marker);
  mfem::Vector u = random_point(), F(u.Size());
  form.Mult(u, F);  // independent of u on an inflow boundary
  CHECK_THAT(global_sum(F), Catch::Matchers::WithinRel(-3., 1.e-12));
}

TEST_CASE_METHOD(SlothSquare, "Advection integrators: invalid data", "[Integrators][errors]") {
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  SECTION("no velocity") {
    SLOTH_CHECK_MFEM_ERROR_WITH(maker<AdvectionNLFormIntegrator>({Coefficients()})(ctx),
                                "no advected variable");
  }
  SECTION("only one velocity component in 2D") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        maker<AdvectionNLFormIntegrator>({Coefficients(Coefficient(Glossary::Velocity_X, 1.))})(
            ctx),
        "only some of the velocity components");
  }
  SECTION("velocity prefactor without velocity_gradient") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        maker<AdvectionNLFormIntegrator>(
            {Coefficients(Coefficient(Glossary::VelocityPrefactor, 1.))})(ctx),
        "velocity_gradient");
  }
  SECTION("boundary without inflow coefficient") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        bdr_maker<AdvectionBoundaryNLFormIntegrator>({velocity(1., 0.)}, 0, kLeft)(ctx),
        "Inflow coefficient not found");
  }
}

// =============================================================================
// Residual of the phase-field integrators against their weak form
// (each term is isolated on fields for which it is computed exactly)
//   Allen-Cahn:    R(u)   = M [ (W'(u), psi) + lambda (grad u, grad psi) ]
//   Cahn-Hilliard: R_phi  = (mu - W'(phi), psi) - lambda (grad phi, grad psi)
//                  R_mu   = (M grad mu, grad psi)
// with W(u) = omega u^2 (1 - u)^2, W'(u) = 2 omega u (1 - u) (1 - 2 u).
// =============================================================================
namespace {
constexpr double kM = 2., kLambda = 0.15, kOmega = 120.;

double dW(double u) { return 2. * kOmega * u * (1. - u) * (1. - 2. * u); }

/// @brief Coefficients of a phase-field integrator; omega = 0 removes the double well
Coefficients pf_coefficients(double omega) {
  return Coefficients(Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(omega)),
                      Coefficient(Glossary::Capillary, kLambda),
                      Coefficient(Glossary::Mobility, kM),
                      Coefficient(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(kLambda)));
}

/// @brief Residual of a block form on the unit square, for given fields (one per block)
struct BlockResidual {
  mfem::Vector F;  // true dofs, blocks one after the other
  int tsize;

  BlockResidual(const sloth_unit::UnitSquare& sq, mfem::BlockNonlinearFormIntegrator* integ,
                const std::vector<double (*)(const mfem::Vector&)>& fields) {
    const int nb = static_cast<int>(fields.size());
    mfem::Array<mfem::ParFiniteElementSpace*> spaces(nb);
    for (int i = 0; i < nb; ++i) spaces[i] = sq.fes.get();
    mfem::ParBlockNonlinearForm form(spaces);
    form.AddDomainIntegrator(integ);
    tsize = sq.fes->GetTrueVSize();
    mfem::Vector x(nb * tsize), xi;
    for (int i = 0; i < nb; ++i) {
      sq.project(fields[static_cast<size_t>(i)]).GetTrueDofs(xi);
      for (int k = 0; k < tsize; ++k) x(i * tsize + k) = xi(k);
    }
    F.SetSize(x.Size());
    form.Mult(x, F);
  }

  mfem::Vector block(int i) const { return mfem::Vector(F.GetData() + i * tsize, tsize); }

  /// @brief (R_i, 1): sum over all the dofs and all the ranks
  double sum(int i) const { return global_sum(block(i)); }

  /// @brief (R_i, x): moment against the interpolant of x
  double moment_x(const sloth_unit::UnitSquare& sq, int i) const {
    mfem::Vector X;
    sq.project([](const mfem::Vector& p) { return p[0]; }).GetTrueDofs(X);
    mfem::Vector bi = block(i);
    return mfem::InnerProduct(MPI_COMM_WORLD, X, bi);
  }
};

double c03(const mfem::Vector&) { return 0.3; }
double c07(const mfem::Vector&) { return 0.7; }
double x_of(const mfem::Vector& p) { return p[0]; }
}  // namespace

TEST_CASE_METHOD(sloth_unit::UnitSquare, "AllenCahnNLFormIntegrator: residual",
                 "[Integrators][residual]") {
  Context ctx;
  ctx.u_old.push_back(project(phi_old));

  SECTION("u = 0.3: (R, 1) = M W'(0.3)") {
    const BlockResidual r(*this, maker<AllenCahnNLFormIntegrator>({pf_coefficients(kOmega)})(ctx),
                          {c03});
    CHECK_THAT(r.sum(0), Catch::Matchers::WithinRel(kM * dW(0.3), 1.e-12));
  }
  SECTION("u = x, no double well: (R, x) = M lambda (grad x, grad x) = M lambda") {
    const BlockResidual r(*this, maker<AllenCahnNLFormIntegrator>({pf_coefficients(0.)})(ctx),
                          {x_of});
    CHECK_THAT(r.moment_x(*this, 0), Catch::Matchers::WithinRel(kM * kLambda, 1.e-12));
  }
}

TEST_CASE_METHOD(sloth_unit::UnitSquare, "CahnHilliardNLFormIntegrator: residual",
                 "[Integrators][residual]") {
  Context ctx;
  ctx.u_old = {project(phi_old), project(phi_old)};

  SECTION("phi = 0.3, mu = 0.7: (R_phi, 1) = mu - W'(phi), R_mu = 0") {
    const BlockResidual r(*this,
                          maker<CahnHilliardNLFormIntegrator>(
                              {pf_coefficients(kOmega), pf_coefficients(kOmega)})(ctx),
                          {c03, c07});
    CHECK_THAT(r.sum(0), Catch::Matchers::WithinRel(0.7 - dW(0.3), 1.e-12));
    CHECK(r.block(1).Normlinf() < 1.e-12);
  }
  SECTION("phi = x, mu = 0, no double well: (R_phi, x) = - lambda") {
    const BlockResidual r(*this,
                          maker<CahnHilliardNLFormIntegrator>(
                              {pf_coefficients(0.), pf_coefficients(0.)})(ctx),
                          {x_of, [](const mfem::Vector&) { return 0.; }});
    CHECK_THAT(r.moment_x(*this, 0), Catch::Matchers::WithinRel(-kLambda, 1.e-12));
  }
  SECTION("phi = 0.3, mu = x: (R_mu, x) = M (grad x, grad x) = M") {
    const BlockResidual r(*this,
                          maker<CahnHilliardNLFormIntegrator>(
                              {pf_coefficients(kOmega), pf_coefficients(kOmega)})(ctx),
                          {c03, x_of});
    CHECK_THAT(r.moment_x(*this, 1), Catch::Matchers::WithinRel(kM, 1.e-12));
  }
}

// -----------------------------------------------------------------------------
// Coefficient of the spatial coordinates in an integrator: the integrator
// gives it the physical point of each integration point.
// Source f = 1 + x + y: |sum_i F_i| = |(f, 1)| = 2 over the unit square.
// -----------------------------------------------------------------------------
namespace {
class OnePlusXPlusY final : public FunctionCoefficient {
 public:
  bool uses_coordinates() const override { return true; }

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [this](const std::span<const double>&, const std::span<const double>&,
                  const std::span<const double>&, const unsigned int) {
      return 1. + this->coordinates_[0] + this->coordinates_[1];
    };
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
};
}  // namespace

TEST_CASE_METHOD(sloth_unit::UnitSquare, "SourceTermNLFormIntegrator: source of the coordinates",
                 "[Integrators]") {
  Context ctx;
  ctx.u_old.push_back(project(phi_old));
  auto form = one_block_form(
      fes.get(), maker<SourceTermNLFormIntegrator>({Coefficients(
                     Coefficient(Glossary::SourceTerm, Scheme::Implicit, OnePlusXPlusY()))})(ctx));
  mfem::Vector x(fes->GetTrueVSize()), y(fes->GetTrueVSize());
  x = 0.5;
  form->Mult(x, y);
  CHECK_THAT(std::abs(global_sum(y)), Catch::Matchers::WithinRel(2., 1.e-12));
}

// -----------------------------------------------------------------------------
// Phases on different finite element spaces (general path of the integrator).
// Phase 0 on Q1, phase 1 on Q2, mobility M = m0 phi0 (1 - phi1) of both phases.
// With fields exact in both spaces (phi0 = 0.3 + 0.1 y, phi1 = 0.4 + 0.2 x), the
// residual is the one computed with both phases on Q1:
//   block 0 (Q1 in both cases): equal dof by dof;
//   block 1: equal moments (R_1, 1) and (R_1, x), 1 and x being in Q1 and Q2.
// The Jacobian (rectangular off-diagonal blocks) matches the residual.
// -----------------------------------------------------------------------------
namespace {
double phi0_affine(const mfem::Vector& p) { return 0.3 + 0.1 * p[1]; }
double phi1_affine(const mfem::Vector& p) { return 0.4 + 0.2 * p[0]; }

/// @brief Interpolant of f on a space, as true dofs
mfem::Vector true_dofs(mfem::ParFiniteElementSpace* fes, double (*f)(const mfem::Vector&)) {
  mfem::ParGridFunction gf(fes);
  mfem::FunctionCoefficient c(f);
  gf.ProjectCoefficient(c);
  mfem::Vector v;
  gf.GetTrueDofs(v);
  return v;
}

/// @brief Two-phase Allen-Cahn form on the spaces (s0, s1); ctx must outlive the form
std::unique_ptr<mfem::ParBlockNonlinearForm> two_phase_form(mfem::ParFiniteElementSpace* s0,
                                                            mfem::ParFiniteElementSpace* s1,
                                                            Context& ctx) {
  for (auto* s : {s0, s1}) {
    mfem::ParGridFunction old(s);
    mfem::FunctionCoefficient c(phi_old);
    old.ProjectCoefficient(c);
    ctx.u_old.push_back(old);
  }
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(2);
  spaces[0] = s0;
  spaces[1] = s1;
  auto form = std::make_unique<mfem::ParBlockNonlinearForm>(spaces);
  form->AddDomainIntegrator(
      maker<AllenCahnNLFormIntegrator>({two_phase_mobility_coefficients<MobilityCross>(),
                                        two_phase_mobility_coefficients<MobilityCross>()})(ctx));
  return form;
}

/// @brief True dofs of (phi0, phi1) on the spaces (s0, s1), blocks one after the other
mfem::Vector two_phase_fields(mfem::ParFiniteElementSpace* s0, mfem::ParFiniteElementSpace* s1) {
  const mfem::Vector x0 = true_dofs(s0, phi0_affine), x1 = true_dofs(s1, phi1_affine);
  mfem::Vector x(x0.Size() + x1.Size());
  x.SetVector(x0, 0);
  x.SetVector(x1, x0.Size());
  return x;
}
}  // namespace

TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "AllenCahnNLFormIntegrator: phases on different finite element spaces",
                 "[Integrators][derivatives]") {
  mfem::H1_FECollection fec2(2, 2);
  mfem::ParFiniteElementSpace fes2(pmesh.get(), &fec2);
  Context ctx_mixed;
  auto mixed = two_phase_form(fes.get(), &fes2, ctx_mixed);
  const int n0 = fes->GetTrueVSize();

  SECTION("residual: same as with both phases on Q1") {
    Context ctx_ref;
    auto ref = two_phase_form(fes.get(), fes.get(), ctx_ref);
    const mfem::Vector x_ref = two_phase_fields(fes.get(), fes.get());
    const mfem::Vector x_mixed = two_phase_fields(fes.get(), &fes2);
    mfem::Vector F_ref(x_ref.Size()), F_mixed(x_mixed.Size());
    ref->Mult(x_ref, F_ref);
    mixed->Mult(x_mixed, F_mixed);

    for (int r = 0; r < n0; ++r) {
      CAPTURE(r);
      CHECK_THAT(F_mixed(r), Catch::Matchers::WithinAbs(F_ref(r), 1.e-12));
    }
    mfem::Vector R1_ref(F_ref.GetData() + n0, F_ref.Size() - n0);
    mfem::Vector R1_mixed(F_mixed.GetData() + n0, F_mixed.Size() - n0);
    CHECK_THAT(global_sum(R1_mixed), Catch::Matchers::WithinAbs(global_sum(R1_ref), 1.e-12));
    const mfem::Vector X_ref = true_dofs(fes.get(), x_of), X_mixed = true_dofs(&fes2, x_of);
    CHECK_THAT(mfem::InnerProduct(MPI_COMM_WORLD, X_mixed, R1_mixed),
               Catch::Matchers::WithinAbs(mfem::InnerProduct(MPI_COMM_WORLD, X_ref, R1_ref),
                                          1.e-12));
  }
  SECTION("Jacobian") {
    mfem::Vector x(n0 + fes2.GetTrueVSize());
    x.Randomize(5);
    x *= 0.6;
    x += 0.2;
    check_gradient(*mixed, x);
  }
}
