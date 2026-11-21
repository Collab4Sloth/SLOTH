/**
 * @file VariationalEval.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Evaluation at a set of points of a variational coefficient f(u, grad u)
 * @version 0.1
 * @date 2026-11-21
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
 * @brief Outputs of the evaluation of a variational coefficient at one point.
 *
 * A variational coefficient ("variational": true in the JSON file) is an energy density
 * f(u, p), u being its n variables and p = grad u the n*d components of their gradients,
 * p_{i,k} = du_i/dx_k stored at the index i*d + k (d: dimension).
 *
 * Each output is computed only if its pointer is not null:
 *   F               f                        1
 *   df_dvar         df/du_i                  n
 *   df_dgrad        df/dp_{i,k}              n*d
 *   d2f_dvar2       d2f/du_i du_j            n*n          (row i, column j)
 *   d2f_dvar_dgrad  d2f/du_i dp_{j,k}        n*(n*d)      (row i, column j*d + k)
 *   d2f_dgrad2      d2f/dp_{i,k} dp_{j,l}    (n*d)*(n*d)  (row i*d + k, column j*d + l)
 */
struct VariationalOutputs {
  double* F = nullptr;
  double* df_dvar = nullptr;
  double* df_dgrad = nullptr;
  double* d2f_dvar2 = nullptr;
  double* d2f_dvar_dgrad = nullptr;
  double* d2f_dgrad2 = nullptr;
};

/**
 * @brief Outputs of the evaluation of a variational coefficient at a set of points (see
 * VariationalOutputs): values of one point are contiguous, e.g. with m = n*d
 *   df_dgrad(r, q) = df_dgrad[q * m + r],  d2f_dgrad2(r, s, q) = d2f_dgrad2[(q * m + r) * m + s]
 * Each output is computed only if its pointer is not null.
 */
struct VariationalBatch {
  mfem::Vector* F = nullptr;
  mfem::Vector* df_dvar = nullptr;
  mfem::Vector* df_dgrad = nullptr;
  mfem::Vector* d2f_dvar2 = nullptr;
  mfem::Vector* d2f_dvar_dgrad = nullptr;
  mfem::Vector* d2f_dgrad2 = nullptr;
};

/**
 * @brief Evaluates a generated variational coefficient at npts points.
 *
 * EVAL is a structure <class_name>VariationalEval written by GenerateCoefficient.py from the JSON
 * file of the coefficients. The loop over the points is an mfem::forall: it runs on the host or
 * on the device according to the mfem::Device configuration, with the same code.
 *
 * Layout of the inputs, as in eval_at_points (Coefficients/EvalAtPoints.hpp): values of one point
 * are contiguous, an input may have more components per point than the structure reads (the
 * first ones are used).
 *   values               variables u_i                        (at least n per point)
 *   explicit_values      explicit variables                    (at least nb_explicit_variables)
 *   auxiliary_values     auxiliary variables                   (at least nb_auxiliary_variables)
 *   gradients            grad u_i, component k at i*d + k      (at least n*d per point)
 *   explicit_gradients   gradients of the explicit variables   (at least nb_explicit_variables*d
 *                        per point if EVAL::uses_explicit_gradients, may be empty otherwise)
 *   auxiliary_gradients  gradients of the auxiliary variables  (at least nb_auxiliary_variables*d
 *                        per point if EVAL::uses_auxiliary_gradients, may be empty otherwise)
 *   coordinates          physical coordinates of the points   (at least nb_coordinates per point)
 *
 * @tparam EVAL Variational evaluation structure generated from the JSON file
 * @param eval Evaluation structure (constants and prefactor)
 * @param npts Number of points
 * @param time Time of the evaluation
 * @param outputs Outputs (see VariationalBatch), each one of size npts times its size per point
 */
template <class EVAL>
void eval_variational_at_points(const EVAL& eval, const int npts, const double time,
                                const mfem::Vector& values, const mfem::Vector& explicit_values,
                                const mfem::Vector& auxiliary_values, const mfem::Vector& gradients,
                                const mfem::Vector& explicit_gradients,
                                const mfem::Vector& auxiliary_gradients,
                                const mfem::Vector& coordinates, const VariationalBatch& outputs) {
  constexpr int n = EVAL::nb_variables;
  constexpr int d = EVAL::dimension;
  constexpr int nd = n * d;
  constexpr int n_exp = EVAL::nb_explicit_variables;
  constexpr int n_aux = EVAL::nb_auxiliary_variables;
  constexpr int n_coord = EVAL::nb_coordinates;

  // Number of components per point of an input, at least the number read by the structure
  auto components = [npts](const mfem::Vector& v, const int required, const char* what) {
    const int m = npts > 0 ? v.Size() / npts : 0;
    MFEM_VERIFY(v.Size() == npts * m && m >= required,
                std::string("eval_variational_at_points: wrong size of ") + what);
    return m;
  };
  const int m_values = components(values, n, "values");
  const int m_exp = components(explicit_values, n_exp, "explicit values");
  const int m_aux = components(auxiliary_values, n_aux, "auxiliary values");
  const int m_grad = components(gradients, nd, "gradients");
  const int m_exp_grad =
      components(explicit_gradients, EVAL::uses_explicit_gradients ? n_exp * d : 0,
                 "gradients of the explicit variables");
  const int m_aux_grad =
      components(auxiliary_gradients, EVAL::uses_auxiliary_gradients ? n_aux * d : 0,
                 "gradients of the auxiliary variables");
  const int dim = components(coordinates, n_coord,
                             "coordinates (spatial coordinates not set or inconsistent with the "
                             "coefficient)");

  auto check_output = [npts](const mfem::Vector* v, const int size, const char* what) {
    MFEM_VERIFY(v == nullptr || v->Size() == npts * size,
                std::string("eval_variational_at_points: wrong size of ") + what);
  };
  check_output(outputs.F, 1, "F");
  check_output(outputs.df_dvar, n, "df_dvar");
  check_output(outputs.df_dgrad, nd, "df_dgrad");
  check_output(outputs.d2f_dvar2, n * n, "d2f_dvar2");
  check_output(outputs.d2f_dvar_dgrad, n * nd, "d2f_dvar_dgrad");
  check_output(outputs.d2f_dgrad2, nd * nd, "d2f_dgrad2");

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
  const double* d_gradients = read(gradients);
  const double* d_explicit_gradients = read(explicit_gradients);
  const double* d_auxiliary_gradients = read(auxiliary_gradients);
  const double* d_coordinates = read(coordinates);
  double* d_F = write(outputs.F);
  double* d_df_dvar = write(outputs.df_dvar);
  double* d_df_dgrad = write(outputs.df_dgrad);
  double* d_d2f_dvar2 = write(outputs.d2f_dvar2);
  double* d_d2f_dvar_dgrad = write(outputs.d2f_dvar_dgrad);
  double* d_d2f_dgrad2 = write(outputs.d2f_dgrad2);

  // Copied by value in the loop (host or device)
  const EVAL e = eval;
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
    VariationalOutputs out;
    out.F = d_F ? d_F + q : nullptr;
    out.df_dvar = d_df_dvar ? d_df_dvar + q * n : nullptr;
    out.df_dgrad = d_df_dgrad ? d_df_dgrad + q * nd : nullptr;
    out.d2f_dvar2 = d_d2f_dvar2 ? d_d2f_dvar2 + q * n * n : nullptr;
    out.d2f_dvar_dgrad = d_d2f_dvar_dgrad ? d_d2f_dvar_dgrad + q * n * nd : nullptr;
    out.d2f_dgrad2 = d_d2f_dgrad2 ? d_d2f_dgrad2 + q * nd * nd : nullptr;
    e(d_values + q * m_values, d_explicit + q * m_exp, d_auxiliary + q * m_aux,
      d_gradients + q * m_grad, d_explicit_gradients + q * m_exp_grad,
      d_auxiliary_gradients + q * m_aux_grad, d_coordinates + q * dim, time, out);
  });
}
