/*
 * AllenCahnSpecializedBilinearFormIntegrator.cpp
 *
 *  Created on: 8 déc. 2021
 *      Author: ci230846
 */

#include <Solvers/AllenCahnSpecializedBilinearFormIntegrator.h>

AllenCahnSpecializedBilinearFormIntegrator::
    AllenCahnSpecializedBilinearFormIntegrator() {
  // TODO Auto-generated constructor stub
}

void AllenCahnSpecializedBilinearFormIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el, mfem::ElementTransformation &Tr,
    mfem::DenseMatrix &elmat) {
  // Linearized DoubleWell contribution
  // W' = W'(phi_n)+W''(phi_n)*(phi_n+1 - phi_n)
  // with :
  // W = phi^2 * (1-phi)^2
  // W' = 2 phi * (1-phi) * (1-2*phi)
  // w'' = 2 * (1 - 6*phi + 6 * phi^2)
}

void AllenCahnSpecializedBilinearFormIntegrator::
    ~AllenCahnSpecializedBilinearFormIntegrator() {
  // TODO Auto-generated destructor stub
}
