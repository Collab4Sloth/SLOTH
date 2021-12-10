#include "mfem.hpp"

using namespace mfem;

// Class used to defined AllenCahn specific BilinearFormIntegrator

bool maxtau = false;  // take a maximum on each element
double vA = 1.0;
double ALPHA = 0.1;  // the parameter in stabilization B terms
bool dtfloor = false;
double dtmin = 0.025;
double dtfactor = 10.;  // a scale factor in tau of dt/dtfactor
// Integrator for (tau * Q.grad func , V.grad v)
// Here we always assume V is the advection speed
class StabConvectionIntegrator : public BilinearFormIntegrator {
 private:
  DenseMatrix dshape, gshape, Jinv, V_ir, Q_ir;
  Coefficient *nuCoef;
  VectorCoefficient *V, *Q;
  double dt;
  int itau;
  double alpha;

 public:
  StabConvectionIntegrator(double dt_, VectorCoefficient &q, int itau_ = 1,
                           double alpha_ = 0.)
      : V(&q), Q(NULL), dt(dt_), itau(itau_), alpha(alpha_){};

  StabConvectionIntegrator(double dt_, VectorCoefficient &q,
                           VectorCoefficient &v, int itau_ = 1,
                           double alpha_ = 0.)
      : V(&v), Q(&q), dt(dt_), itau(itau_), alpha(alpha_){};

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Tr,
                                     DenseMatrix &elmat);

  virtual ~StabConvectionIntegrator(){};
};

void StabConvectionIntegrator::AssembleElementMatrix(const FiniteElement &el,
                                                     ElementTransformation &Tr,
                                                     DenseMatrix &elmat) {
  double norm, tau;
  double Unorm, invtau;
  int dim = 2;
  int nd = el.GetDof();
  Vector advGrad(nd), advGrad2(nd), vec1(dim), vec2(dim);

  dshape.SetSize(nd, dim);
  gshape.SetSize(nd, dim);
  Jinv.SetSize(dim);

  elmat.SetSize(nd);
  elmat = 0.;

  double eleLength = sqrt(Geometry::Volume[el.GetGeomType()] * Tr.Weight());
  // integration order is el.order + grad.order-1 (-1 due to another derivative
  // taken in V)
  int intorder = 1 * (el.GetOrder() + Tr.OrderGrad(&el) - 1);
  const IntegrationRule &ir = IntRules.Get(el.GetGeomType(), intorder);

  V->Eval(V_ir, Tr, ir);
  if (Q != NULL) Q->Eval(Q_ir, Tr, ir);

  for (int i = 0; i < ir.GetNPoints(); i++) {
    const IntegrationPoint &ip = ir.IntPoint(i);
    el.CalcDShape(ip, dshape);

    Tr.SetIntPoint(&ip);
    norm = ip.weight * Tr.Weight();
    CalcInverse(Tr.Jacobian(), Jinv);
    Mult(dshape, Jinv, gshape);

    V_ir.GetColumnReference(i, vec1);

    // compute tau
    double nu = 0.;
    Unorm = vec1.Norml2();
    if (itau == 1) {
      invtau = sqrt(pow(2.0 * Unorm / eleLength,
                        2) /*+ pow(4.0*nu/(eleLength*eleLength),2)*/);
    } else if (itau == 2) {
      // invtau = sqrt( pow(2./.1,2) +  pow(2.0 * Unorm / eleLength, 2) );
      invtau = sqrt(pow(2. / dt, 2) +
                    pow(2.0 * Unorm / eleLength,
                        2) /*+ pow(4.0*nu/(eleLength*eleLength),2)*/);
    } else if (itau == 3) {
      if (alpha > 0.)
        invtau = sqrt(pow(alpha / dt, 2) +
                      pow(2.0 * Unorm / eleLength,
                          2) /*+ pow(4.0*nu/(eleLength*eleLength),2)*/);
      else
        invtau = sqrt(pow(dtfactor / dt, 2) +
                      pow(2.0 * Unorm / eleLength,
                          2) /*+ pow(4.0*nu/(eleLength*eleLength),2)*/);
    } else if (itau == 4) {
      invtau = 2. / dt;
    } else
      invtau = 2.0 / dt +
               2.0 * Unorm / eleLength /*+ 4.0 * nu / (eleLength * eleLength)*/;
    tau = 1.0 / invtau;

    norm *= tau;

    gshape.Mult(vec1, advGrad);
    if (Q == NULL) {
      AddMult_a_VVt(-norm, advGrad, elmat);
    } else {
      Q_ir.GetColumnReference(i, vec2);
      gshape.Mult(vec2, advGrad2);
      AddMult_a_VWt(norm, advGrad, advGrad2, elmat);
    }
  }
}

// Integrator for (tau * rhs , V\cdot grad u)
class StabDomainLFIntegrator : public LinearFormIntegrator {
 private:
  DenseMatrix dshape, gshape, Jinv, V_ir;
  VectorCoefficient *V;
  Coefficient &Q;
  double dt;
  int itau;
  double alpha;

 public:
  StabDomainLFIntegrator(double dt_, double visc, VectorCoefficient &q,
                         Coefficient &QF, int itau_ = 1, double alpha_ = 0.)
      : V(&q), Q(QF), dt(dt_), itau(itau_), alpha(alpha_){};

  virtual void AssembleRHSElementVect(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      Vector &elvect);

  virtual ~StabDomainLFIntegrator(){};
};

void StabDomainLFIntegrator::AssembleRHSElementVect(const FiniteElement &el,
                                                    ElementTransformation &Tr,
                                                    Vector &elvect) {
  double norm, tau;
  double Unorm, invtau;
  int dim = 2;
  int nd = el.GetDof();
  Vector advGrad(nd), vec1(dim);

  dshape.SetSize(nd, dim);
  gshape.SetSize(nd, dim);
  Jinv.SetSize(dim);

  // here we assume 2d quad
  double eleLength = sqrt(Geometry::Volume[el.GetGeomType()] * Tr.Weight());

  elvect.SetSize(nd);
  elvect = 0.;

  int intorder = 2 * el.GetOrder() + Tr.OrderGrad(&el) - 1;
  const IntegrationRule &ir = IntRules.Get(el.GetGeomType(), intorder);

  V->Eval(V_ir, Tr, ir);

  // compare maximum tau
  // cout <<"tauMax="<<tauMax<<" length="<<eleLength<<" "<<"dt="<<dt<<"
  // U="<<Unorm<<" ";

  for (int i = 0; i < ir.GetNPoints(); i++) {
    const IntegrationPoint &ip = ir.IntPoint(i);

    el.CalcDShape(ip, dshape);

    Tr.SetIntPoint(&ip);
    norm = ip.weight * Tr.Weight();
    CalcInverse(Tr.Jacobian(), Jinv);
    Mult(dshape, Jinv, gshape);

    V_ir.GetColumnReference(i, vec1);
    // compute tau
    double nu = 0.;
    Unorm = vec1.Norml2();
    if (itau == 1) {
      invtau = sqrt(pow(2.0 * Unorm / eleLength, 2) +
                    pow(4.0 * nu / (eleLength * eleLength), 2));
      // invtau = 2.0 * Unorm / eleLength + 4.0 * nu / (eleLength * eleLength);
    } else if (itau == 2)
      invtau = sqrt(pow(2. / dt, 2) + pow(2.0 * Unorm / eleLength, 2) +
                    pow(4.0 * nu / (eleLength * eleLength), 2));
    else if (itau == 3) {
      if (alpha > 0.)
        invtau = sqrt(pow(alpha / dt, 2) + pow(2.0 * Unorm / eleLength, 2) +
                      pow(4.0 * nu / (eleLength * eleLength), 2));
      else
        invtau = sqrt(pow(dtfactor / dt, 2) + pow(2.0 * Unorm / eleLength, 2) +
                      pow(4.0 * nu / (eleLength * eleLength), 2));
    } else if (itau == 4)
      invtau = 2. / dt;
    else
      invtau = 2.0 / dt + 2.0 * Unorm / eleLength +
               4.0 * nu / (eleLength * eleLength);
    tau = 1.0 / invtau;

    norm *= (tau * Q.Eval(Tr, ip));

    gshape.Mult(vec1, advGrad);

    add(elvect, norm, advGrad, elvect);
  }
}
