/**
 * @file SlothGridFunction.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-12-5
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <chrono>
#include <memory>
#include <tuple>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class SlothGridFunction : public mfem::ParGridFunction {
 private:
 protected:
 public:
  using mfem::ParGridFunction::ParGridFunction;
  SlothGridFunction() : mfem::ParGridFunction() {}
  SlothGridFunction(const mfem::ParGridFunction &orig) : mfem::ParGridFunction(orig) {}
  void GetGradient(mfem::ElementTransformation &Tr, mfem::DenseMatrix &gradPsi, mfem::Vector &grad);
  // ~SlothGridFunction();
};


/**
 * @brief
 *
 * 
 */
void SlothGridFunction::GetGradient(mfem::ElementTransformation &Tr, mfem::DenseMatrix &gradPsi,
                                    mfem::Vector &grad) {
  mfem::Vector elfun;
  int num = Tr.ElementNo;
  GetElementDofValues(num, elfun);
  gradPsi.MultTranspose(elfun, grad);
}

// SlothGridFunction::~SlothGridFunction() { mfem::~ParGridFunction()}
