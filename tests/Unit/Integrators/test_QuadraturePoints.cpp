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
#include <string>
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

/// @brief Allen-Cahn integrator giving access to the evaluation on the local mesh
class MeshProbe : public AllenCahnNLFormIntegrator<VARS2D> {
 public:
  using AllenCahnNLFormIntegrator<VARS2D>::AllenCahnNLFormIntegrator;
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
    // As AllenCahnNLFormIntegrator: element tables, then the coefficient
    this->prepare_element_values(el, Tr, elfun, aux, ir);
    this->eval_coefficient_at_points(coef, var_blk, outputs, fe, Tr, ir, elfun, aux, out,
                                     values_at);
    return true;
  }
};

/// @brief Allen-Cahn integrator declaring its coefficients on the local mesh
class AllenCahnOnMesh : public AllenCahnNLFormIntegrator<VARS2D> {
 public:
  using AllenCahnNLFormIntegrator<VARS2D>::AllenCahnNLFormIntegrator;

  void init() override {
    AllenCahnNLFormIntegrator<VARS2D>::init();
    const unsigned int value = Base::CoefValue;
    const unsigned int first = Base::CoefFirstDerivative;
    const unsigned int second = Base::CoefSecondDerivative;
    for (int blk = 0; blk < static_cast<int>(this->nb_blk_); ++blk) {
      this->prepare_on_mesh(this->mobility[blk], blk, value, value | first, this->vaux_gf_);
      this->prepare_on_mesh(this->lambda[blk], blk, value, value, this->vaux_gf_);
      this->prepare_on_mesh(this->double_well_energy[blk], blk, first, first | second,
                            this->vaux_gf_);
    }
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
      {"W, implicit, block 0", Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(3.)), 0,
       false, true},
      {"W, implicit, block 1", Coefficient(Glossary::FreeEnergy, Scheme::Implicit, W(3.)), 1,
       false, true},
      {"H, explicit", Coefficient(Glossary::FreeEnergy, Scheme::Explicit, H()), 0, false, true},
      {"MobilityCross, implicit",
       Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityCross()), 0, false, true},
      {"MobilityCross, explicit",
       Coefficient(Glossary::Mobility, Scheme::Explicit, MobilityCross()), 0, false, true},
      {"SemiImplicitEnergy, semi-implicit",
       Coefficient(Glossary::FreeEnergy, Scheme::SemiImplicit, SemiImplicitEnergy()), 0, false,
       true},
      {"MobilityX, coordinates", Coefficient(Glossary::Mobility, Scheme::Implicit, MobilityX()),
       1, false, true},
      {"FirstAuxiliary, auxiliary field",
       Coefficient(Glossary::Mobility, Scheme::Implicit, FirstAuxiliary()), 0, true, true},
      {"constant", Coefficient(Glossary::Mobility, 3.5), 0, false, false},
  };

  using Assembly = Base::Assembly;
  const unsigned int all =
      MeshProbe::kValue | MeshProbe::kFirst | MeshProbe::kSecond;
  for (const bool at_end_of_step : {false, true}) {
    const auto values_at = at_end_of_step ? MeshProbe::ValuesAt::EndOfStep
                                          : MeshProbe::ValuesAt::Unknowns;
    for (const MeshCase& c : cases) {
      INFO(c.name << (at_end_of_step ? " (end of step)" : ""));
      MeshProbe probe(Geometry::Cartesian, 0.1, u_old, {}, kNoParameters, {},
                      {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.)),
                       phase_field_coefficients(Coefficient(Glossary::Mobility, 1.))});
      probe.init();
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
        CHECK(!probe.at_element(coef, 1 - c.var_blk, all, a, values_at, *fes, u, e, true,
                                unused));
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
        REQUIRE(probe.at_element(coef, c.var_blk, all, a, values_at, *fes, u, e, false,
                                 by_element));
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

  MeshProbe probe(Geometry::Cartesian, 0.1, u_old, {}, kNoParameters, {},
                  {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.))});
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

// Assembly with the coefficients of Allen-Cahn prepared on the local mesh: same residual and
// Jacobian as mfem::ParBlockNonlinearForm with the evaluation element by element
TEST_CASE_METHOD(sloth_unit::UnitSquare,
                 "AllenCahn: same residual and Jacobian with the coefficients prepared on the mesh",
                 "[Integrators][coefficients]") {
  struct AllenCahnCase {
    std::string name;
    std::vector<Coefficients> coefficients;
  };
  const Coefficient mobility_phi(Glossary::Mobility, Scheme::Implicit, MobilityPhi());
  const Coefficient mobility_cross(Glossary::Mobility, Scheme::Implicit, MobilityCross());
  const Coefficients coupled(Coefficient(Glossary::FreeEnergy, Scheme::Implicit, CoupledEnergy()),
                             Coefficient(Glossary::Capillary, 0.15), mobility_cross,
                             Coefficient(Glossary::GradEnergy, Scheme::Implicit,
                                         GradientEnergy(0.15)));
  const std::vector<AllenCahnCase> cases{
      {"1 phase, M constant", {phase_field_coefficients(Coefficient(Glossary::Mobility, 1.))}},
      {"1 phase, M(phi)", {phase_field_coefficients(mobility_phi)}},
      {"2 phases, M(phi0, phi1), coupled energy", {coupled, coupled}},
      {"2 phases, W(phi_b) and M = H(phi_b)",
       {phase_field_coefficients(Coefficient(Glossary::Mobility, Scheme::Implicit, H(2.))),
        phase_field_coefficients(Coefficient(Glossary::Mobility, Scheme::Implicit, H(2.)))}},
  };

  for (const AllenCahnCase& c : cases) {
    INFO(c.name);
    const int nb = static_cast<int>(c.coefficients.size());
    std::vector<mfem::ParGridFunction> u_old{project(phi0_old)};
    if (nb == 2) u_old.push_back(project(phi1_old));
    mfem::Array<mfem::ParFiniteElementSpace*> spaces(nb);
    for (int b = 0; b < nb; ++b) spaces[b] = fes.get();

    mfem::ParBlockNonlinearForm reference(spaces);
    auto* by_element = new AllenCahnNLFormIntegrator<VARS2D>(
        Geometry::Cartesian, 0.1, u_old, {}, kNoParameters, {}, c.coefficients);
    by_element->init();
    reference.AddDomainIntegrator(by_element);

    SlothParBlockNonlinearForm<VARS2D> form(spaces);
    auto* on_mesh = new AllenCahnOnMesh(Geometry::Cartesian, 0.1, u_old, {}, kNoParameters, {},
                                        c.coefficients);
    on_mesh->init();
    REQUIRE(on_mesh->uses_quadrature_points());
    form.AddDomainIntegrator(on_mesh);

    // Values in [0.2, 0.8]
    mfem::Vector x(nb * fes->GetTrueVSize()), v(x.Size());
    x.Randomize(5);
    x *= 0.6;
    x += 0.2;
    v.Randomize(7);

    mfem::Vector y_reference(x.Size()), y(x.Size());
    reference.Mult(x, y_reference);
    form.Mult(x, y);
    CHECK(max_difference(y_reference, y) == 0.);

    mfem::Vector jv_reference(x.Size()), jv(x.Size());
    reference.GetGradient(x).Mult(v, jv_reference);
    form.GetGradient(x).Mult(v, jv);
    CHECK(max_difference(jv_reference, jv) == 0.);
  }
}
