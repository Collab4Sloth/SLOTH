/**
 * @file test_QuadraturePoints.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests of the coefficients evaluated at all the quadrature points of the local mesh
 *        (SlothNLFormIntegrator::prepare_on_mesh, prepare_quadrature_points)
 * @version 0.1
 * @date 2026-11-15
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
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Sloth/sloth.hpp"
#include "Sloth/tests.hpp"
#include "sloth_unit.hpp"

// Coefficients generated at configure time from Integrators/MobilityCoefficients.json
#include "CoupledEnergy.hpp"
#include "FirstAuxiliary.hpp"
#include "MobilityCross.hpp"
#include "MobilityPhi.hpp"
#include "MobilityX.hpp"
#include "SemiImplicitEnergy.hpp"

namespace {
using VARS2D = Sloth2D::VARS;
using Base = SlothNLFormIntegrator<VARS2D>;

const Parameters kNoParameters{std::vector<Parameter>{}};

double phi0_old(const mfem::Vector& p) {
  return 0.5 + 0.2 * std::sin(3. * p[0]) * std::cos(2. * p[1]);
}
double phi1_old(const mfem::Vector& p) {
  return 0.4 + 0.15 * std::cos(2. * p[0]) * std::sin(3. * p[1]);
}
double phi0_now(const mfem::Vector& p) { return 0.3 + 0.2 * p[0] * p[1]; }
double phi1_now(const mfem::Vector& p) { return 0.6 - 0.25 * p[0] + 0.1 * p[1]; }

/// @brief Largest absolute difference between two vectors of values
double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  REQUIRE(a.size() == b.size());
  double d = 0.;
  for (std::size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - b[i]));
  return d;
}

/// @brief Largest absolute difference between two vectors of true dofs (all ranks)
double max_difference(const mfem::Vector& a, const mfem::Vector& b) {
  REQUIRE(a.Size() == b.Size());
  double local = 0., global = 0.;
  for (int i = 0; i < a.Size(); ++i) local = std::max(local, std::abs(a(i) - b(i)));
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

/// @brief Integrator declaring no coefficient on the local mesh (Fourier), giving access to the
/// evaluation on the local mesh of the coefficients declared by the test
class MeshProbe : public FourierNLFormIntegrator<VARS2D> {
 public:
  explicit MeshProbe(const std::vector<mfem::ParGridFunction>& u_old)
      : FourierNLFormIntegrator<VARS2D>(
            Geometry::Cartesian, 0.1, u_old, {}, kNoParameters, {},
            std::vector<Coefficients>(u_old.size(),
                                      Coefficients(Coefficient(Glossary::Conductivity, 1.)))) {}
  using Base::CoefficientAtPoints;
  using Base::ValuesAt;
  static constexpr unsigned int kValue = Base::CoefValue;
  static constexpr unsigned int kFirst = Base::CoefFirstDerivative;
  static constexpr unsigned int kSecond = Base::CoefSecondDerivative;

  void declare(Coefficient& coef, const int var_blk, const unsigned int outputs,
               const std::vector<mfem::ParGridFunction>& aux, const ValuesAt values_at) {
    this->prepare_on_mesh(coef, var_blk, outputs, outputs, aux, values_at);
  }

  /**
   * @brief Values of coef on element e: read in the prepared values (prepared = true), or
   * evaluated element by element as during the assembly (prepared = false)
   * @return false if prepared values are requested but not available
   */
  bool at_element(Coefficient& coef, const int var_blk, const unsigned int outputs,
                  const std::vector<mfem::ParGridFunction>& aux, const ValuesAt values_at,
                  mfem::ParFiniteElementSpace& fes, const std::vector<mfem::ParGridFunction>& u,
                  const int e, const bool prepared, CoefficientAtPoints& out) {
    const int nb = static_cast<int>(u.size());
    std::vector<mfem::Vector> dofs(u.size());
    mfem::Array<const mfem::Vector*> elfun(nb);
    mfem::Array<const mfem::FiniteElement*> el(nb);
    for (int b = 0; b < nb; ++b) {
      u[static_cast<std::size_t>(b)].GetElementDofValues(e, dofs[static_cast<std::size_t>(b)]);
      elfun[b] = &dofs[static_cast<std::size_t>(b)];
      el[b] = fes.GetFE(e);
    }
    const mfem::FiniteElement& fe = *fes.GetFE(e);
    mfem::ElementTransformation& Tr = *fes.GetElementTransformation(e);
    const mfem::IntegrationRule& ir = this->integration_rule(fe, Tr);
    if (prepared) {
      return this->read_prepared_coefficient(coef, var_blk, outputs, aux, values_at, e, ir, out);
    }
    // As during the assembly: element tables, then the coefficient
    this->prepare_element_values(el, Tr, elfun, aux, ir);
    this->eval_coefficient_at_points(coef, var_blk, outputs, fe, Tr, ir, elfun, aux, out,
                                     values_at);
    return true;
  }
};

/// @brief One coefficient to prepare: name, coefficient, block, auxiliary field, prepared
struct MeshCase {
  std::string name;
  Coefficient coef;
  int var_blk;
  bool uses_aux;
  bool prepared;  ///< false: not prepared (evaluated element by element during the assembly)
};

Coefficients phase_field_coefficients(Coefficient mobility) {
  const double lambda = 0.15, omega = 120.;
  return Coefficients(Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(omega)),
                      Coefficient(Glossary::Capillary, lambda), mobility,
                      Coefficient(Glossary::GradEnergy, Scheme::Implicit, GradientEnergy(lambda)));
}
}  // namespace

// Values prepared on the local mesh = values evaluated element by element (identical), for the
// residual and the Jacobian outputs, the values of the unknowns and at the end of the step
TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "prepare_quadrature_points: same values as the evaluation element by element",
                 "[Integrators][coefficients]") {
  std::vector<mfem::ParGridFunction> u_old{project(phi0_old), project(phi1_old)};
  std::vector<mfem::ParGridFunction> u{project(phi0_now), project(phi1_now)};
  const std::vector<mfem::ParGridFunction> aux{project(phi1_old)};
  const std::vector<mfem::ParGridFunction> no_aux;

  // Local dofs of the two blocks
  mfem::Array<int> offsets(3);
  offsets[0] = 0;
  offsets[1] = fes->GetVSize();
  offsets[2] = 2 * fes->GetVSize();
  mfem::BlockVector x(offsets);
  x.GetBlock(0) = u[0];
  x.GetBlock(1) = u[1];
  mfem::Array<mfem::FiniteElementSpace*> spaces(2);
  spaces[0] = spaces[1] = fes.get();

  const std::vector<MeshCase> cases{
      {"W, implicit, block 0", Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(3.)), 0, false,
       true},
      {"W, implicit, block 1", Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(3.)), 1, false,
       true},
      {"H, explicit", Coefficient(Glossary::FreeEnergy, Scheme::Explicit, H()), 0, false, true},
      {"MobilityCross, implicit",
       Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityCross()), 0, false, true},
      {"MobilityCross, explicit",
       Coefficient(Glossary::Mobility, Scheme::Explicit, MobilityCross()), 0, false, true},
      {"SemiImplicitEnergy, semi-implicit",
       Coefficient(Glossary::FreeEnergy, Scheme::SemiImplicit, SemiImplicitEnergy()), 0, false,
       true},
      {"MobilityX, coordinates", Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityX()), 1,
       false, true},
      {"FirstAuxiliary, auxiliary field",
       Coefficient(Glossary::Mobility, Scheme::Implicit, FirstAuxiliary()), 0, true, true},
      {"constant", Coefficient(Glossary::Mobility, 3.5), 0, false, false},
  };

  using Assembly = Base::Assembly;
  const unsigned int all = MeshProbe::kValue | MeshProbe::kFirst | MeshProbe::kSecond;
  for (const bool at_end_of_step : {false, true}) {
    const auto values_at =
        at_end_of_step ? MeshProbe::ValuesAt::EndOfStep : MeshProbe::ValuesAt::Unknowns;
    for (const MeshCase& c : cases) {
      INFO(c.name << (at_end_of_step ? " (end of step)" : ""));
      MeshProbe probe(u_old);
      probe.init();
      REQUIRE(!probe.uses_quadrature_points());
      Coefficient coef = c.coef;
      const std::vector<mfem::ParGridFunction>& a = c.uses_aux ? aux : no_aux;
      probe.declare(coef, c.var_blk, all, a, values_at);
      CHECK(probe.uses_quadrature_points() == c.prepared);

      const int ne = fes->GetNE();
      std::vector<MeshProbe::CoefficientAtPoints> on_mesh(static_cast<std::size_t>(ne));
      probe.prepare_quadrature_points(spaces, x, Assembly::Jacobian);
      for (int e = 0; e < ne; ++e) {
        CAPTURE(e);
        auto& out = on_mesh[static_cast<std::size_t>(e)];
        REQUIRE(probe.at_element(coef, c.var_blk, all, a, values_at, *fes, u, e, true, out) ==
                c.prepared);
        // Other block: not read
        MeshProbe::CoefficientAtPoints unused;
        CHECK(!probe.at_element(coef, 1 - c.var_blk, all, a, values_at, *fes, u, e, true, unused));
      }

      // After the assembly: nothing is read, the coefficient is evaluated element by element
      probe.clear_quadrature_points();
      MeshProbe::CoefficientAtPoints cleared;
      CHECK(!probe.at_element(coef, c.var_blk, all, a, values_at, *fes, u, 0, true, cleared));
      if (!c.prepared) continue;
      for (int e = 0; e < ne; ++e) {
        CAPTURE(e);
        const auto& prepared = on_mesh[static_cast<std::size_t>(e)];
        MeshProbe::CoefficientAtPoints by_element;
        REQUIRE(
            probe.at_element(coef, c.var_blk, all, a, values_at, *fes, u, e, false, by_element));
        CHECK(prepared.nb == by_element.nb);
        CHECK(prepared.is_explicit == by_element.is_explicit);
        CHECK(max_difference(prepared.F, by_element.F) == 0.);
        CHECK(max_difference(prepared.dF, by_element.dF) == 0.);
        CHECK(max_difference(prepared.d2F, by_element.d2F) == 0.);
      }
    }
  }
}

// Residual outputs only: the outputs of the Jacobian are not read
TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "prepare_quadrature_points: outputs of the residual and of the Jacobian",
                 "[Integrators][coefficients]") {
  std::vector<mfem::ParGridFunction> u_old{project(phi0_old)};
  std::vector<mfem::ParGridFunction> u{project(phi0_now)};
  const std::vector<mfem::ParGridFunction> no_aux;
  mfem::Array<int> offsets(2);
  offsets[0] = 0;
  offsets[1] = fes->GetVSize();
  mfem::BlockVector x(offsets);
  x.GetBlock(0) = u[0];
  mfem::Array<mfem::FiniteElementSpace*> spaces(1);
  spaces[0] = fes.get();

  MeshProbe probe(u_old);
  probe.init();
  Coefficient coef(Glossary::FreeEnergy, Scheme::Implicit, W(3.));
  const auto unknowns = MeshProbe::ValuesAt::Unknowns;
  const unsigned int first = MeshProbe::kFirst;
  const unsigned int second = MeshProbe::kSecond;
  probe.declare(coef, 0, first, no_aux, unknowns);  // same outputs for both assemblies

  MeshProbe::CoefficientAtPoints out;
  probe.prepare_quadrature_points(spaces, x, Base::Assembly::Residual);
  CHECK(probe.at_element(coef, 0, first, no_aux, unknowns, *fes, u, 0, true, out));
  CHECK(!probe.at_element(coef, 0, first | second, no_aux, unknowns, *fes, u, 0, true, out));
}

// -----------------------------------------------------------------------------
// Integrators assembled by SlothParBlockNonlinearForm (coefficients prepared on the local mesh):
// same residual and Jacobian as with mfem::ParBlockNonlinearForm (evaluation element by element)
// -----------------------------------------------------------------------------
namespace {
/// @brief One integrator to assemble: name, number of blocks, how to build it
struct AssemblyCase {
  std::string name;
  int nb_blocks;
  std::function<SlothNLFormIntegrator<VARS2D>*(const std::vector<mfem::ParGridFunction>&)> make;
};

/// @brief Builds an integrator with its coefficients (one Coefficients per block)
template <template <class> class INTEGRATOR>
std::function<SlothNLFormIntegrator<VARS2D>*(const std::vector<mfem::ParGridFunction>&)> maker(
    std::vector<Coefficients> coefficients) {
  return [coefficients](const std::vector<mfem::ParGridFunction>& u_old) {
    SlothNLFormIntegrator<VARS2D>* integrator = new INTEGRATOR<VARS2D>(
        Geometry::Cartesian, 0.1, u_old, {}, kNoParameters, {}, coefficients);
    integrator->init();
    return integrator;
  };
}

/**
 * @brief Assembles the integrators built by make with mfem::ParBlockNonlinearForm (evaluation
 * element by element: the hook is not called) and with SlothParBlockNonlinearForm (coefficients
 * prepared on the local mesh), and checks that the residuals and the products J v are identical.
 */
void check_same_assembly(mfem::Array<mfem::ParFiniteElementSpace*>& spaces,
                         const std::function<SlothNLFormIntegrator<VARS2D>*()>& make) {
  mfem::ParBlockNonlinearForm reference(spaces);
  reference.AddDomainIntegrator(make());
  SlothParBlockNonlinearForm<VARS2D> form(spaces);
  form.AddDomainIntegrator(make());

  // Values in [0.2, 0.8]
  int size = 0;
  for (int b = 0; b < spaces.Size(); ++b) size += spaces[b]->GetTrueVSize();
  mfem::Vector x(size), v(size);
  x.Randomize(5);
  x *= 0.6;
  x += 0.2;
  v.Randomize(7);

  mfem::Vector y_reference(size), y(size);
  reference.Mult(x, y_reference);
  form.Mult(x, y);
  CHECK(max_difference(y_reference, y) == 0.);

  mfem::Vector jv_reference(size), jv(size);
  reference.GetGradient(x).Mult(v, jv_reference);
  form.GetGradient(x).Mult(v, jv);
  CHECK(max_difference(jv_reference, jv) == 0.);
}

std::vector<AssemblyCase> assembly_cases() {
  const Coefficient mobility_phi(Glossary::Mobility, Scheme::Implicit, MobilityPhi());
  const Coefficient mobility_cross(Glossary::Mobility, Scheme::Implicit, MobilityCross());
  const Coefficients coupled(Coefficient(Glossary::FreeEnergy, Scheme::Implicit, CoupledEnergy()),
                             Coefficient(Glossary::Capillary, 0.15), mobility_cross,
                             Coefficient(Glossary::GradEnergy, Scheme::Implicit,
                                         GradientEnergy(0.15)));
  const auto h_mobility = [](const double p) {
    return phase_field_coefficients(Coefficient(Glossary::Mobility, Scheme::Implicit, H(p)));
  };
  const auto diffusivity = [](const Coefficient& d) { return Coefficients(d); };
  return {
      {"AllenCahn, 1 phase, M constant", 1,
       maker<AllenCahnNLFormIntegrator>(
           {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.))})},
      {"AllenCahn, 1 phase, M(phi)", 1,
       maker<AllenCahnNLFormIntegrator>({phase_field_coefficients(mobility_phi)})},
      {"AllenCahn, 2 phases, M(phi0, phi1), coupled energy", 2,
       maker<AllenCahnNLFormIntegrator>({coupled, coupled})},
      {"AllenCahn, 2 phases, W(phi_b) and M = H(phi_b)", 2,
       maker<AllenCahnNLFormIntegrator>({h_mobility(2.), h_mobility(2.)})},
      {"Fick, D constant", 1,
       maker<FickNLFormIntegrator>({diffusivity(Coefficient(Glossary::Diffusivity, 2.))})},
      {"Fick, D(u)", 1,
       maker<FickNLFormIntegrator>(
           {diffusivity(Coefficient(Glossary::Diffusivity, Scheme::Implicit, MobilityPhi()))})},
      {"Fick, 2 blocks, D = H(u_b) per block", 2,
       maker<FickNLFormIntegrator>(
           {diffusivity(Coefficient(Glossary::Diffusivity, Scheme::Implicit, H(2.))),
            diffusivity(Coefficient(Glossary::Diffusivity, Scheme::Implicit, H(2.)))})},
      {"Fick, 2 blocks, D(u_0, u_1) explicit", 2,
       maker<FickNLFormIntegrator>(
           {diffusivity(Coefficient(Glossary::Diffusivity, Scheme::Explicit, MobilityCross())),
            diffusivity(Coefficient(Glossary::Diffusivity, Scheme::Implicit, MobilityCross()))})},
      {"Fourier, lambda(T)", 1,
       maker<FourierNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::Conductivity, Scheme::Implicit, H(3.)))})},
      // Time-derivative integrators: the unknown is k = du/dt, implicit coefficients evaluated
      // at the end of the step u = u_old + dt k
      {"TimeDerivative, A = B = 1", 1, maker<TimeNLFormIntegrator>({Coefficients()})},
      {"TimeDerivative, 2 blocks, implicit A(u_b) and B(u_b)", 2,
       maker<TimeNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::ExplicitTime_A, Scheme::Implicit, H(2.)),
                         Coefficient(Glossary::ExplicitTime_B, Scheme::Implicit, W(5.))),
            Coefficients(Coefficient(Glossary::ExplicitTime_A, Scheme::Implicit, H(3.)),
                         Coefficient(Glossary::ExplicitTime_B, 2.))})},
      {"TimeDerivative, 2 blocks, A(u_0, u_1) implicit", 2,
       maker<TimeNLFormIntegrator>(
           {Coefficients(
                Coefficient(Glossary::ExplicitTime_A, Scheme::Implicit, MobilityCross())),
            Coefficients(
                Coefficient(Glossary::ExplicitTime_A, Scheme::Explicit, MobilityCross()))})},
      {"HeatTimeDerivative, explicit Cp(T)", 1,
       maker<HeatTimeNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::Concentration, 2.),
                         Coefficient(Glossary::Cp, Scheme::Explicit, W(5.)))})},
      {"HeatTimeDerivative, implicit rho(T) and Cp(T)", 1,
       maker<HeatTimeNLFormIntegrator>(
           {Coefficients(Coefficient(Glossary::Concentration, Scheme::Implicit, H(2.)),
                         Coefficient(Glossary::Cp, Scheme::Implicit, W(5.)))})},
      {"SplitTimeDerivative, A = 1", 2,
       maker<TimeCHNLFormIntegrator>({Coefficients(), Coefficients()})},
      {"CahnHilliard, W(phi), M constant", 2,
       maker<CahnHilliardNLFormIntegrator>(
           {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.)),
            phase_field_coefficients(Coefficient(Glossary::Mobility, 1.))})},
      {"CahnHilliard, mobility M(phi)", 2,
       maker<CahnHilliardNLFormIntegrator>(
           {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.)),
            phase_field_coefficients(mobility_phi)})},
      {"CahnHilliard, degenerate mobility H(phi)", 2,
       maker<CahnHilliardNLFormIntegrator>(
           {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.)), h_mobility(2.)})},
      {"SplitAllenCahn, W(phi), M constant", 2,
       maker<BlockAllenCahnNLFormIntegrator>(
           {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.)),
            phase_field_coefficients(Coefficient(Glossary::Mobility, 1.))})},
      {"SplitAllenCahn, mobility M(phi)", 2,
       maker<BlockAllenCahnNLFormIntegrator>(
           {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.)),
            phase_field_coefficients(mobility_phi)})},
      {"SplitTimeDerivative, implicit A(phi)", 2,
       maker<TimeCHNLFormIntegrator>(
           {Coefficients(),
            Coefficients(Coefficient(Glossary::ExplicitTime_A, Scheme::Implicit, H(2.)))})},
  };
}
}  // namespace

TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "Integrators: same residual and Jacobian with the coefficients prepared on the mesh",
                 "[Integrators][coefficients]") {
  for (const AssemblyCase& c : assembly_cases()) {
    INFO(c.name);
    const int nb = c.nb_blocks;
    std::vector<mfem::ParGridFunction> u_old{project(phi0_old)};
    if (nb == 2) u_old.push_back(project(phi1_old));
    mfem::Array<mfem::ParFiniteElementSpace*> spaces(nb);
    for (int b = 0; b < nb; ++b) spaces[b] = fes.get();
    check_same_assembly(spaces, [&c, &u_old]() { return c.make(u_old); });
  }
}

// -----------------------------------------------------------------------------
// Integrators with auxiliary variables: SLOTH Variables recognised by their additional
// information, built on the mesh of a SpatialDiscretization
// -----------------------------------------------------------------------------
namespace {
struct AuxiliarySquare {
  using SPA2D = Sloth2D::SPA;
  using BCS2D = Sloth2D::BCS;
  using VAR2D = Sloth2D::VAR;
  static constexpr int DIM2 = 2;

  std::unique_ptr<SPA2D> spatial =
      std::make_unique<SPA2D>("InlineSquareWithQuadrangles", 1, 0, std::make_tuple(4, 4, 1., 1.));
  std::unique_ptr<BCS2D> bcs = std::make_unique<BCS2D>(
      spatial.get(), Boundary("lower", 0, "Neumann"), Boundary("right", 1, "Neumann"),
      Boundary("upper", 2, "Neumann"), Boundary("left", 3, "Neumann"));

  mfem::ParFiniteElementSpace* fes() const { return spatial->get_finite_element_space(); }

  mfem::ParGridFunction project(double (*f)(const mfem::Vector&)) const {
    mfem::ParGridFunction gf(fes());
    mfem::FunctionCoefficient c(f);
    gf.ProjectCoefficient(c);
    return gf;
  }

  /// @brief Auxiliary variable with its additional informations (e.g. {"U", "mu"})
  VARS2D aux_variable(const std::string& name, const std::vector<std::string>& infos,
                      GlossaryQuantity qty,
                      std::function<double(const mfem::Vector&, double)> f) const {
    auto v = VAR2D(spatial.get(), *bcs, name, qty, 2, AnalyticalFunctions<DIM2>(f));
    for (const auto& info : infos) v.set_additional_information(info);
    return VARS2D(v);
  }
};
}  // namespace

// Mass diffusion flux M grad(mu) with a stabilisation coefficient D (prepared on the mesh)
TEST_CASE_METHOD(AuxiliarySquare,
                 "MassDiffusionFlux: same residual and Jacobian with the coefficients prepared on "
                 "the mesh",
                 "[Integrators][coefficients]") {
  VARS2D mu = aux_variable("mu_U", {"U", "mu"}, Glossary::ChemicalPotential,
                           [](const mfem::Vector& p, double) { return p[0] * p[1]; });
  VARS2D mob = aux_variable("M_U", {"U", "inter_mob"}, Glossary::InterDiffusionMobility,
                            [](const mfem::Vector& p, double) { return 2. + p[0]; });
  const std::vector<VARS2D*> auxvars{&mu, &mob};
  const std::vector<mfem::ParGridFunction> aux_old{mu.getVariables()[0].get_gf(),
                                                   mob.getVariables()[0].get_gf()};
  const std::vector<mfem::ParGridFunction> u_old{project(phi0_old)};
  mfem::Array<mfem::ParFiniteElementSpace*> spaces(1);
  spaces[0] = fes();

  const std::vector<std::pair<std::string, Coefficient>> stabilisations{
      {"D constant", Coefficient(Glossary::Diffusivity, 1.)},
      {"D = H(u)", Coefficient(Glossary::Diffusivity, Scheme::Implicit, H(3.))},
      {"D(u) generated", Coefficient(Glossary::Diffusivity, Scheme::Implicit, MobilityPhi())},
  };
  for (const auto& stabilisation : stabilisations) {
    INFO(stabilisation.first);
    const Coefficient& diffusivity = stabilisation.second;
    check_same_assembly(spaces, [&]() {
      SlothNLFormIntegrator<VARS2D>* integrator = new MassDiffusionFluxNLFormIntegrator<VARS2D>(
          Geometry::Cartesian, 0.1, u_old, aux_old, kNoParameters, auxvars,
          {Coefficients(diffusivity)});
      integrator->init();
      return integrator;
    });
  }
}
