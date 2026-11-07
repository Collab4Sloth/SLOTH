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

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Evaluates a generated coefficient at npts points.
 *
 * EVAL is a structure <class_name>Eval written by GenerateCoefficient.py from the JSON file of
 * the coefficients. The loop over the points is an mfem::forall: it runs on the host or on the
 * device according to the mfem::Device configuration, with the same code.
 *
 * Layout of the data: values of one point are contiguous (point by point), e.g.
 *   values(k, q) = values[q * EVAL::nb_variables + k]
 *   hessian(i, j, q) = hessian[q * n * n + i * n + j],  n = EVAL::nb_variables
 * This is the layout of an mfem::QuadratureFunction whose vdim is the number of components.
 *
 * Each output is computed only if its pointer is not null.
 *
 * @tparam EVAL Evaluation structure generated from the JSON file
 * @param eval Evaluation structure (constants and prefactor)
 * @param npts Number of points
 * @param time Time of the evaluation
 * @param values Variables at the points (size npts * nb_variables)
 * @param explicit_values Explicit variables at the points (size npts * nb_explicit_variables)
 * @param auxiliary_values Auxiliary variables at the points (size npts * nb_auxiliary_variables)
 * @param coordinates Physical coordinates of the points (size npts * nb_coordinates)
 * @param F Values of the coefficient (size npts), or nullptr
 * @param gradient Gradients of the coefficient (size npts * n), or nullptr
 * @param hessian Hessians of the coefficient (size npts * n * n), or nullptr
 */
template <class EVAL>
void eval_at_points(const EVAL& eval, const int npts, const double time,
                    const mfem::Vector& values, const mfem::Vector& explicit_values,
                    const mfem::Vector& auxiliary_values, const mfem::Vector& coordinates,
                    mfem::Vector* F, mfem::Vector* gradient = nullptr,
                    mfem::Vector* hessian = nullptr) {
  constexpr int n = EVAL::nb_variables;
  constexpr int n_exp = EVAL::nb_explicit_variables;
  constexpr int n_aux = EVAL::nb_auxiliary_variables;
  constexpr int n_coord = EVAL::nb_coordinates;

  MFEM_VERIFY(values.Size() == npts * n, "eval_at_points: wrong size of values");
  MFEM_VERIFY(explicit_values.Size() == npts * n_exp,
              "eval_at_points: wrong size of explicit values");
  MFEM_VERIFY(auxiliary_values.Size() == npts * n_aux,
              "eval_at_points: wrong size of auxiliary values");
  MFEM_VERIFY(coordinates.Size() == npts * n_coord, "eval_at_points: wrong size of coordinates");
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
    e(d_values + q * n, d_explicit + q * n_exp, d_auxiliary + q * n_aux,
      d_coordinates + q * n_coord, time, d_F ? d_F + q : nullptr,
      d_gradient ? d_gradient + q * n : nullptr, d_hessian ? d_hessian + q * n * n : nullptr);
  });
}
