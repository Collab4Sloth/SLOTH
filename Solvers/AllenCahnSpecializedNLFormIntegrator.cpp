/*
 * AllenCahnSpecializedNLFormIntegrator.cpp
 *
 *  Created on: 8 déc. 2021
 *      Author: ci230846
 */

#include <Solvers/AllenCahnSpecializedNLFormIntegrator.h>

void AllenCahnSpecializedNLFormIntegrator::
    AllenCahnSpecializedNLFormIntegrator() {}

void AllenCahnSpecializedNLFormIntegrator::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
    const mfem::Vector& elfun, mfem::Vector& elvect) {
  // Linearized DoubleWell contribution
  // W' = W'(phi_n)+W''(phi_n)*(phi_n+1 - phi_n)
  // with :
  // W = phi^2 * (1-phi)^2
  // W' = 2 phi * (1-phi) * (1-2*phi)
  // w'' = 2 * (1 - 6*phi + 6 * phi^2)
  auto e = 0;
}

void AllenCahnSpecializedNLFormIntegrator::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
    const mfem::Vector& elfun, mfem::DenseMatrix& elmat) {
  // Linearized DoubleWell contribution
  // W' = W'(phi_n)+W''(phi_n)*(phi_n+1 - phi_n)
  // with :
  // W = phi^2 * (1-phi)^2
  // W' = 2 phi * (1-phi) * (1-2*phi)
  // w'' = 2 * (1 - 6*phi + 6 * phi^2)
  auto e = 0;
}

AllenCahnSpecializedNLFormIntegrator::~AllenCahnSpecializedNLFormIntegrator() {
  // TODO Auto-generated destructor stub
}
