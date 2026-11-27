/**
 * @file EvalAtPoints.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Evaluation of a generated coefficient at a set of points, on the host or on the device
 * @version 0.1
 * @date 2026-11-07
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
#pragma once

#include <string>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Evaluates a generated coefficient at npts points.
 *
 * EVAL is a structure <class_name>Eval written by GenerateCoefficient.py from the JSON file of
 * the coefficients. The loop over the points is an mfem::forall: it runs on the host or on the
 * device according to the mfem::Device configuration, with the same code.
 *
 * Layout of the data: values of one point are contiguous (point by point), e.g.
 *   values(k, q) = values[q * m + k],  m = values.Size() / npts
 *   hessian(i, j, q) = hessian[q * n * n + i * n + j],  n = EVAL::nb_variables
 * This is the layout of an mfem::QuadratureFunction whose vdim is the number of components.
 *
 * An input may have more components per point than the evaluation structure reads (e.g. all the
 * blocks or all the auxiliary fields of an integrator): the first ones are used, as in the point
 * by point evaluation. The outputs have the size of the structure: n derivatives per point.
 *
 * Each output is computed only if its pointer is not null.
 *
 * @tparam EVAL Evaluation structure generated from the JSON file
 * @param eval Evaluation structure (constants and prefactor)
 * @param npts Number of points
 * @param time Time of the evaluation
 * @param values Variables at the points (at least nb_variables per point)
 * @param explicit_values Explicit variables at the points (at least nb_explicit_variables per
 *        point)
 * @param auxiliary_values Auxiliary variables at the points (at least nb_auxiliary_variables per
 *        point)
 * @param coordinates Physical coordinates of the points (size npts * dim, with dim at least
 *        nb_coordinates: e.g. the dimension of the mesh)
 * @param F Values of the coefficient (size npts), or nullptr
 * @param gradient Gradients of the coefficient (size npts * n), or nullptr
 * @param hessian Hessians of the coefficient (size npts * n * n), or nullptr
 */
template <class EVAL>
void eval_at_points(const EVAL& eval, const int npts, const double time, const mfem::Vector& values,
                    const mfem::Vector& explicit_values, const mfem::Vector& auxiliary_values,
                    const mfem::Vector& coordinates, mfem::Vector* F,
                    mfem::Vector* gradient = nullptr, mfem::Vector* hessian = nullptr) {
  constexpr int n = EVAL::nb_variables;
  constexpr int n_exp = EVAL::nb_explicit_variables;
  constexpr int n_aux = EVAL::nb_auxiliary_variables;
  constexpr int n_coord = EVAL::nb_coordinates;

  // Number of components per point of an input, at least the number read by the structure
  auto components = [npts](const mfem::Vector& v, const int required, const char* what) {
    const int m = npts > 0 ? v.Size() / npts : 0;
    MFEM_VERIFY(v.Size() == npts * m && m >= required,
                std::string("eval_at_points: wrong size of ") + what);
    return m;
  };
  const int m_values = components(values, n, "values");
  const int m_exp = components(explicit_values, n_exp, "explicit values");
  const int m_aux = components(auxiliary_values, n_aux, "auxiliary values");
  const int dim = components(coordinates, n_coord,
                             "coordinates (spatial coordinates not set or inconsistent with the "
                             "coefficient)");
  MFEM_VERIFY(F == nullptr || F->Size() == npts, "eval_at_points: wrong size of F");
  MFEM_VERIFY(gradient == nullptr || gradient->Size() == npts * n,
              "eval_at_points: wrong size of gradient");
  MFEM_VERIFY(hessian == nullptr || hessian->Size() == npts * n * n,
              "eval_at_points: wrong size of hessian");

  // Data are read and written where the loop runs (host or device). An empty vector (no
  // explicit variable, no point...) is not accessed: nullptr.
  auto read = [](const mfem::Vector& v) -> const double* {
    if (v.Size() == 0) return nullptr;
    v.UseDevice(true);
    return v.Read();
  };
  auto write = [](mfem::Vector* v) -> double* {
    if (v == nullptr || v->Size() == 0) return nullptr;
    v->UseDevice(true);
    return v->Write();
  };
  const double* d_values = read(values);
  const double* d_explicit = read(explicit_values);
  const double* d_auxiliary = read(auxiliary_values);
  const double* d_coordinates = read(coordinates);
  double* d_F = write(F);
  double* d_gradient = write(gradient);
  double* d_hessian = write(hessian);

  // Copied by value in the loop (host or device)
  const EVAL e = eval;
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
    e(d_values + q * m_values, d_explicit + q * m_exp, d_auxiliary + q * m_aux,
      d_coordinates + q * dim, time, d_F ? d_F + q : nullptr,
      d_gradient ? d_gradient + q * n : nullptr, d_hessian ? d_hessian + q * n * n : nullptr);
  });
}
