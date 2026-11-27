/**
 * @file QuadraturePointsKernels.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Loops (host or device) of the evaluation at the quadrature points of the local mesh
 * @version 0.1
 * @date 2026-11-17
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
 * @brief mfem::forall loops of SlothNLFormIntegrator (prepare_quadrature_points).
 *
 * Free functions: with CUDA, a host / device lambda cannot be defined in a private or protected
 * member function. The pointers are read or written where the loops run (Read / Write of the
 * mfem::Vector by the caller). Layout point by point, npts points.
 */
namespace QuadraturePointsKernels {

/// {u, u_old} of block b at the points: u_uold[q * 2 nb + b], u_uold[q * 2 nb + nb + b].
/// End of step (integrators of the time derivative, u = k): u_old + dt k
inline void fill_u_uold(const int npts, const int nb, const int b, const bool end_of_step,
                        const double dt, const double* u, const double* u_old, double* u_uold) {
  const int width = 2 * nb;
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
    u_uold[q * width + b] = end_of_step ? u_old[q] + dt * u[q] : u[q];
    u_uold[q * width + nb + b] = u_old[q];
  });
}

/// Component k of ncomp at the points: to[q * ncomp + k] = from[q]
inline void fill_component(const int npts, const int ncomp, const int k, const double* from,
                           double* to) {
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) { to[q * ncomp + k] = from[q]; });
}

/// Coordinates point by point (dim per point) from the geometric factors (NQ x SDIM x NE)
inline void coordinates_by_points(const int npts, const int nq, const int dim, const double* x,
                                  double* x_points) {
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int p) {
    const int e = p / nq;
    const int q = p % nq;
    for (int d = 0; d < dim; ++d) {
      x_points[p * dim + d] = x[q + nq * (d + dim * e)];
    }
  });
}

/// Values given to a coefficient (n per point) selected from {u, u_old}: u (implicit) or u_old
/// (explicit), and u_old for the explicit values of a semi-implicit coefficient; block var_blk
/// only for a coefficient of one variable
inline void select_values(const int npts, const int nb, const int n, const bool single,
                          const int var_blk, const bool previous, const bool semi,
                          const double* u_uold, double* v, double* v_explicit) {
  const int width = 2 * nb;
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
    for (int k = 0; k < n; ++k) {
      const int b = single ? var_blk : k;
      v[q * n + k] = previous ? u_uold[q * width + nb + b] : u_uold[q * width + b];
      if (semi) v_explicit[q * n + k] = u_uold[q * width + nb + b];
    }
  });
}

/// First derivatives with respect to the nc variables of the coefficient stored with respect to
/// the nb blocks: dF[q * nb + b] (dF set to zero by the caller)
inline void scatter_first_derivatives(const int npts, const int nb, const int nc, const bool single,
                                      const int var_blk, const double* gradient, double* dF) {
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
    for (int k = 0; k < nc; ++k) {
      const int b = single ? var_blk : k;
      dF[q * nb + b] = gradient[q * nc + k];
    }
  });
}

/// Second derivatives stored with respect to the blocks: d2F[(q * nb + i) * nb + j] (d2F set to
/// zero by the caller)
inline void scatter_second_derivatives(const int npts, const int nb, const int nc,
                                       const bool single, const int var_blk, const double* hessian,
                                       double* d2F) {
  mfem::forall(npts, [=] MFEM_HOST_DEVICE(int q) {
    for (int k = 0; k < nc; ++k) {
      const int bk = single ? var_blk : k;
      for (int l = 0; l < nc; ++l) {
        const int bl = single ? var_blk : l;
        d2F[(q * nb + bk) * nb + bl] = hessian[(q * nc + k) * nc + l];
      }
    }
  });
}

}  // namespace QuadraturePointsKernels