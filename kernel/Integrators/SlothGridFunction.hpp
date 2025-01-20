/**
 * @file SlothGridFunction.hpp
 * @author  cp273896  (clement.plumecocq@cea.fr)
 * @brief Class dedicated to the calculation of the gradient of a GridFunction
 * @version 0.1
 * @date 2025-01-14
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <memory>
#include <tuple>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class SlothGridFunction : public mfem::ParGridFunction {
 private:
 protected:
 public:
  SlothGridFunction();
  explicit SlothGridFunction(const mfem::ParGridFunction &orig);
  void GetGradient(mfem::ElementTransformation &Tr, mfem::DenseMatrix &gradPsi, mfem::Vector &grad);
  ~SlothGridFunction();
};

/**
 * @brief Construct a new Sloth Grid Function:: Sloth Grid Function object
 *
 */
SlothGridFunction::SlothGridFunction() : mfem::ParGridFunction() {}

/**
 * @brief Construct a new Sloth Grid Function:: Sloth Grid Function object
 *
 * @param gf
 */
SlothGridFunction::SlothGridFunction(const mfem::ParGridFunction &gf) : mfem::ParGridFunction(gf) {}

/**
 * @brief Compute the gradient of a given GridFunction
 *
 * @param Tr
 * @param gradPsi
 * @param grad
 */
void SlothGridFunction::GetGradient(mfem::ElementTransformation &Tr, mfem::DenseMatrix &gradPsi,
                                    mfem::Vector &grad) {
  mfem::Vector elfun;
  const int &num = Tr.ElementNo;
  GetElementDofValues(num, elfun);
  gradPsi.MultTranspose(elfun, grad);
}

/**
 * @brief Destroy the Sloth Grid Function:: Sloth Grid Function object
 *
 */
SlothGridFunction::~SlothGridFunction() {}
