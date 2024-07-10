/**
 * @file AllenCahnNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <tuple>

#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseFieldMobilities.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]

using FuncType = std::function<double(const double&, const double&)>;
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
class AllenCahnNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  mfem::GridFunction u_old_;
  mfem::DenseMatrix dshape, dshapedxt, invdfdx;
  mfem::Vector shape, vec, pointflux;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;

  MobilityFunctions<MOBI> mobility_function_;
  FuncType laplacian();
  FuncType double_well_derivative(const int& order_derivative);

 protected:
  double omega_, lambda_, mob_;

  virtual FuncType energy_derivatives(const int& order_derivative);

 public:
  AllenCahnNLFormIntegrator(const mfem::GridFunction& u_old, const double& omega,
                            const double& lambda, const double& mob);
  ~AllenCahnNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun, mfem::Vector& elvect);

  virtual void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun, mfem::DenseMatrix& elmat);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Energy derivatives contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return FuncType
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FuncType AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::energy_derivatives(
    const int& order_derivative) {
  return [this, order_derivative](const double& u, const double& un) {
    return this->double_well_derivative(order_derivative)(u, un);
  };
}

/**
 * @brief double_well_derivative contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return std::function<double(const double&, const double&)>
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FuncType AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::double_well_derivative(
    const int& order_derivative) {
  return FuncType([this, order_derivative](const double& u, const double& un) {
    std::function<double(const double&)> W_derivative;
    if (order_derivative == 1) {
      W_derivative = this->energy_first_derivative_potential_.getPotentialFunction(un);
    } else if (order_derivative == 2) {
      W_derivative = this->energy_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const auto& Mphi = this->mobility_function_.getMobilityFunction(un);
    const auto& w_prime = Mphi(this->mob_) * this->omega_ * W_derivative(u);
    return w_prime;
  });
}

/**
 * @brief Laplacian coefficient
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @return std::function<double(const double&, const double&)>
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FuncType AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::laplacian() {
  return FuncType([this](const double& u, const double& un) {
    const auto& Mphi = this->mobility_function_.getMobilityFunction(un);
    const auto& laplacian = Mphi(this->mob_) * this->lambda_;
    return laplacian;
  });
}

/**
 * @brief Construct a new AllenCahnNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AllenCahnNLFormIntegrator(
    const mfem::GridFunction& u_old, const double& omega, const double& lambda, const double& mob)
    : u_old_(u_old), omega_(omega), lambda_(lambda), mob_(mob) {}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  Catch_Time_Section("AllenCahnNLFormIntegrator:AssembleElementVector");

  int nd = el.GetDof();
  int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();
  dshape.SetSize(nd, dim);
  shape.SetSize(nd);
  invdfdx.SetSize(dim, spaceDim);
  vec.SetSize(dim);
  pointflux.SetSize(spaceDim);

  elvect.SetSize(nd);
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());
  elvect = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcDShape(ip, dshape);  // dphi
    el.CalcShape(ip, shape);    // phi
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * shape;
    const auto& un = this->u_old_.GetValue(Tr, ip);

    // Given phi, compute (w'(phi), v), v is shape function
    const double& ww = ip.weight * Tr.Weight() * this->energy_derivatives(1)(u, un);
    add(elvect, ww, shape, elvect);

    // Laplacian : given u, compute (grad(u), grad(v)), v is shape function.
    CalcAdjugate(Tr.Jacobian(), invdfdx);  // invdfdx = adj(J)
    dshape.MultTranspose(elfun, vec);
    invdfdx.MultTranspose(vec, pointflux);
    double w;
    w = this->laplacian()(u, un) * ip.weight / Tr.Weight();
    pointflux *= w;
    invdfdx.Mult(pointflux, vec);

    //
    dshape.AddMult(vec, elvect);
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::DenseMatrix& elmat) {
  Catch_Time_Section("AllenCahnNLFormIntegrator::AssembleElementGrad");

  int nd = el.GetDof();
  int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();
  bool square = (dim == spaceDim);
  double w;

  shape.SetSize(nd);
  dshape.SetSize(nd, dim);
  dshapedxt.SetSize(nd, spaceDim);
  elmat.SetSize(nd);

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  elmat = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcDShape(ip, dshape);  // dphi
    const auto& u = elfun * shape;
    const auto& un = this->u_old_.GetValue(Tr, ip);

    Tr.SetIntPoint(&ip);
    w = Tr.Weight();  // det(J)
    // std::cout << " SQUARE  ? " << square << std::endl;
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise

    // Tr.AdjugateJacobian() det(J)J-1

    w *= this->laplacian()(u, un);

    // dshapedxt =  det(J)J-1 dshape
    Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);
    // elmat += w * dshapedxt * dshapedxt^T
    AddMult_a_AAt(w, dshapedxt, elmat);

    // Compute w'(u)*(du,v), v is shape function ( // w''(u))
    double fun_val = this->energy_derivatives(2)(u, un) * ip.weight * Tr.Weight();
    // elmat += fun_val * shape * shape^T
    AddMult_a_VVt(fun_val, shape, elmat);  // w'(u)*(du, v)
  }
}

/**
 * @brief Destroy the AllenCahnNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::~AllenCahnNLFormIntegrator() {}
