/*
 * AllenCahnSpecializedBilinearFormIntegrator.h
 *
 *  Created on: 8 déc. 2021
 *      Author: ci230846
 */
#include "mfem.hpp"

#ifndef SOLVERS_ALLENCAHNSPECIALIZEDBILINEARFORMINTEGRATOR_H_
#define SOLVERS_ALLENCAHNSPECIALIZEDBILINEARFORMINTEGRATOR_H_

class AllenCahnSpecializedBilinearFormIntegrator
    : public mfem::BilinearFormIntegrator {
 private:
  mfem::DenseMatrix dshape, gshape, Jinv, V_ir, Q_ir;
  mfem::Coefficient *Q;

 public:
  AllenCahnSpecializedBilinearFormIntegrator(mfem::Coefficient &q) : Q(&q){};
  virtual void AssembleElementMatrix(const mfem::FiniteElement &el,
                                     mfem::ElementTransformation &Tr,
                                     mfem::DenseMatrix &elmat);
  virtual ~AllenCahnSpecializedBilinearFormIntegrator();
};

#endif /* SOLVERS_ALLENCAHNSPECIALIZEDBILINEARFORMINTEGRATOR_H_ \
          */
