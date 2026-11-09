/**
 * @file UserIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Extension point: nonlinear form integrators defined outside the SLOTH kernel
 * @version 0.1
 * @date 2026-11-09
 *
 * @copyright CEA (C) 2026
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
#pragma once

#include <functional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "Coefficients/Coefficients.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/ProblemsOptions.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Arguments provided by an operator to build a user integrator: those of the constructor
 *        of the SLOTH integrators.
 *
 * @details The references are only valid during the construction of the integrator.
 */
template <class VARS>
struct UserIntegratorArgs {
  Geometry geometry;
  double time_step;
  const std::vector<mfem::ParGridFunction>& u_old;
  const std::vector<mfem::ParGridFunction>& aux_old;
  const Parameters& params;
  std::vector<VARS*> auxvars;
  const std::vector<Coefficients>& coefficients;
};

/**
 * @brief Integrator defined outside the SLOTH kernel: a name and a function building it.
 *
 * @details The name identifies the integrator in the messages of SLOTH. It must not be the name
 *          of a built-in integrator.
 */
template <class VARS>
class UserIntegrator {
 public:
  using Builder = std::function<SlothNLFormIntegrator<VARS>*(const UserIntegratorArgs<VARS>&)>;

  UserIntegrator(std::string name, Builder builder)
      : name_(std::move(name)), builder_(std::move(builder)) {}

  const std::string& name() const { return name_; }

  SlothNLFormIntegrator<VARS>* build(const UserIntegratorArgs<VARS>& args) const {
    return builder_(args);
  }

 private:
  std::string name_;
  Builder builder_;
};

/**
 * @brief User integrator whose class follows the constructor convention of the SLOTH
 *        integrators (geometry, time step, previous unknowns, previous auxiliary fields,
 *        parameters, auxiliary variables, coefficients).
 *
 * @code
 *   TransientOPE oper(spatials, {"AllenCahn", user_integrator<MyNLFormIntegrator, VARS>("Mine")},
 *                     TimeScheme::EulerImplicit, "TimeDerivative");
 * @endcode
 */
template <template <class> class Integrator, class VARS>
UserIntegrator<VARS> user_integrator(const std::string& name) {
  return UserIntegrator<VARS>(
      name, [](const UserIntegratorArgs<VARS>& args) -> SlothNLFormIntegrator<VARS>* {
        return new Integrator<VARS>(args.geometry, args.time_step, args.u_old, args.aux_old,
                                    args.params, args.auxvars, args.coefficients);
      });
}

/**
 * @brief Integrator of an operator: the name of a built-in integrator, or a user integrator.
 *
 * @details A list of names ({"AllenCahn", "Fick"}) is a valid list of entries.
 */
template <class VARS>
using IntegratorEntry = std::variant<std::string, UserIntegrator<VARS>>;

/**
 * @brief Name of an integrator entry (built-in or user integrator).
 */
template <class VARS>
const std::string& integrator_name(const IntegratorEntry<VARS>& entry) {
  if (const auto* user = std::get_if<UserIntegrator<VARS>>(&entry)) {
    return user->name();
  }
  return std::get<std::string>(entry);
}
