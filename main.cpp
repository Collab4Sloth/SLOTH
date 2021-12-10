
#include <fstream>
#include <iostream>
#include "Solvers/AllenCahnSpecializedNLFormIntegrator.h"
#include "mfem.hpp"

using namespace std;
using namespace mfem;

////
/// Includes à mettre ailleurs, à généraliser, à optimiser...
///
double InitialOrderParameter(const Vector &x);
double ExactSolution(const Vector &x);
double radius;
/// TODO
/// // vérifier l'existence commune et mettre dans un fichier différent sinon
class NLOperator : public Operator {
 private:
  NonlinearForm *N;
  mutable SparseMatrix *Jacobian;
  // f in F(u) = -Laplace u + u^2

 public:
  NLOperator(NonlinearForm *N_, int size)
      : Operator(size), N(N_), Jacobian(NULL) {}

  // - Residual
  virtual void Mult(const Vector &x, Vector &y) const { N->Mult(x, y); }

  // Jacobian
  virtual Operator &GetGradient(const Vector &x) const {
    Jacobian = dynamic_cast<SparseMatrix *>(&N->GetGradient(x));
    return *Jacobian;
  }
};

int main(int argc, char *argv[]) {
  // 1. Parse command-line options.
  const char *mesh_file = "mesh_for_test.msh";
  int order = 1;
  bool static_cond = false;
  bool pa = false;
  const char *device_config = "cpu";
  bool visualization = true;
  bool algebraic_ceed = false;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree) or -1 for"
                 " isoparametric space.");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                 "--no-partial-assembly", "Enable Partial Assembly.");
  args.AddOption(&device_config, "-d", "--device",
                 "Device configuration string, see Device::Configure().");
#ifdef MFEM_USE_CEED
  args.AddOption(&algebraic_ceed, "-a", "--algebraic", "-no-a",
                 "--no-algebraic", "Use algebraic Ceed solver");
#endif
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(cout);
    return 1;
  }
  args.PrintOptions(cout);

  // 1. Enable hardware devices such as GPUs, and programming models such as
  //    CUDA, OCCA, RAJA and OpenMP based on command line options.
  Device device(device_config);
  device.Print();

  // 2. Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
  //    the same code.
  Mesh mesh(mesh_file, 1, 1);
  int dim = mesh.Dimension();
  std::cout << "Dimension " << dim << std::endl;

  // 4. Refine the mesh to increase the resolution. In this example we do
  //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
  //    largest number that gives a final mesh with no more than 50,000
  //    elements.
  {
    int ref_levels = 6;
    for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
    }
  }

  // 5. Define a finite element space on the mesh. Here we use continuous
  //    Lagrange finite elements of the specified order. If order < 1, we
  //    instead use an isoparametric/isogeometric space.
  H1_FECollection fec(order, dim);
  FiniteElementSpace fespace(&mesh, &fec);
  auto size = fespace.GetVSize();
  mfem::Vector minCoord, maxCoord;
  mesh.GetBoundingBox(minCoord, maxCoord);
  std::cout << "mesh size x = " << minCoord[0] << " " << maxCoord[0]
            << std::endl;
  std::cout << "mesh size y = " << minCoord[1] << " " << maxCoord[1]
            << std::endl;
  std::cout << "Number of finite element unknowns: " << fespace.GetTrueVSize()
            << std::endl;
  radius = maxCoord[0] - minCoord[0];

  // Coefficients
  auto _mobility(1.e4);
  ConstantCoefficient mobility(_mobility);
  auto _sigma(0.1);
  auto _epsilon(2.e-1);
  auto _omega = 24. * _sigma / _epsilon;
  ConstantCoefficient omega(_omega);
  auto _lambda = 1.5 * _sigma * _epsilon;
  ConstantCoefficient lambda(_lambda);
  auto _meltingForce(0.);
  ConstantCoefficient calphadContribution(_meltingForce);

  // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
  //    In this example, the boundary conditions are defined by marking
  //    the boundary attributes from the mesh as essential (Dirichlet), Natural
  //    (Neumann) and converting them to a list of true dofs.

  // Define two arrays containing boundar attributes numbers stored
  // "marker arrays" are used to define the portions of boundary associated
  //    with each type of boundary condition. These arrays have an entry
  //    corresponding to each boundary attribute.  Placing a '1' in entry i
  //    marks attribute i+1 as being active, '0' is inactive.
  // Here 1,2 NeumanHomogène, 3 Dirichlet
  Array<int> Dirichlet_bdr(mesh.bdr_attributes.Max());
  //  std::cout << " mesh.bdr_attributes.Max() " << mesh.bdr_attributes.Max()
  //            << std::endl;
  Dirichlet_bdr = 0;
  //  Dirichlet_bdr[0] = 1;
  //  Dirichlet_bdr[1] = 1;
  //  Dirichlet_bdr[2] = 1;
  // Dirichlet_bdr[3] = 1;
  auto _phiD = 1.;
  // Essential BdC = Dirichlet
  // Set the Dirichlet values in the solution vector
  //  Vector p_init(mesh.bdr_attributes.Max());
  //  p_init(0) = 1.00;
  //  p_init(1) = _phiD;
  //  p_init(2) = _phiD;
  // p_init(3) = _phiD;
  //  PWConstCoefficient dbr_value(p_init);
  //   4. Extract the list of all the boundary DOFs. These will be marked as
  //      Dirichlet in order to enforce  boundary conditions.
  Array<int> essential_tdof_list(0);
  fespace.GetEssentialTrueDofs(Dirichlet_bdr, essential_tdof_list);

  /////
  /// dphi/dt = -M [ Calphad + W'(phi) - Div(lambda Grad phi) ]
  /// -----------------------
  /// Variational formulation :
  /// -----------------------
  /// < dphi/dt, v >  + <M W'(phi), v > + < M lambda Grad phi, Grad v > =
  /// <-M Calphad, v> + < M  (n . lambda Grad phi), v >
  ///=============
  /// BilinearForm
  /// ------------
  /// < M lambda Grad phi, Grad v > := Diffusion Integrator
  /// <M W'(phi), v > :=
  ///=============
  /// LinearForm
  /// <-M Calphad, v> := DomainLFIntegrator(-M Calphad)
  ///============
  ///

  // GridFunctions definition
  GridFunction phi_GF(&fespace);

  // Set initial conditions
  FunctionCoefficient phi_0(InitialOrderParameter);
  phi_GF.ProjectCoefficient(phi_0);

  Vector phi;
  phi_GF.GetTrueDofs(phi);  // CCI pourquoi cette étape ?

  //////
  ParaViewDataCollection *pd = NULL;

  // 15. Save data in the ParaView format
  ParaViewDataCollection paraview_dc("MainPST", &mesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order);
  paraview_dc.SetCycle(0);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetTime(0.0);  // set the time
  paraview_dc.RegisterField("phi", &phi_GF);
  paraview_dc.Save();
  /////////////////////////
  // Start Newton
  /////////////////////////
  NonlinearForm N(&fespace);
  auto NN = new AllenCahnSpecializedNLFormIntegrator(calphadContribution,
                                                     mobility, lambda, omega);
  N.AddDomainIntegrator(NN);
  N.SetEssentialTrueDofs(essential_tdof_list);

  NLOperator N_oper(&N, size);

  Solver *J_solver;
  Solver *J_prec = new DSmoother(1);
  MINRESSolver *J_minres = new MINRESSolver;
  J_minres->SetPreconditioner(*J_prec);
  J_solver = J_minres;

  NewtonSolver newton_solver;
  newton_solver.SetSolver(*J_solver);
  newton_solver.SetOperator(N_oper);
  newton_solver.SetPrintLevel(2);

  GridFunction phi_h(&fespace);
  Vector zero;
  newton_solver.Mult(zero, phi_h);

  /////////////////////////
  // End Newton
  /////////////////////////
  FunctionCoefficient phi_exact_coeff(ExactSolution);
  cout << "L2 error norm: " << phi_h.ComputeL2Error(phi_exact_coeff) << endl;
  //

  paraview_dc.RegisterField("phi", &phi_GF);
  paraview_dc.SetCycle(1);
  paraview_dc.SetTime(1.0);
  paraview_dc.Save();

  ////////////////////////////////////////
  // Time loop
  //////////////////
  ///
  /// Assuming oper be a TimeDependent operator
  TimeDependentOperator oper;
  /// output of oper -> phi
  /// TODO add calculation of phi
  phi_GF.SetFromTrueDofs(phi);

  auto ode_solver_type = 1;
  ODESolver *ode_solver;
  switch (ode_solver_type) {
    // Implicit method
    case 1:
      ode_solver = new BackwardEulerSolver;
      break;
    // Explicit methods
    case 2:
      ode_solver = new ForwardEulerSolver;
      break;
    case 3:
      ode_solver = new RK4Solver;
      break;
    default:
      cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
      delete mesh;
      return 3;
  }

  double t = 0.0;
  // TODO mettre ces variables en arguments, ou dans un jeu de données...
  auto dt = 0.1;
  auto t_final = 1.0;
  // Initialization
  oper.SetTime(t);
  ode_solver->Init(oper);
  // 8. Perform time-integration (looping over the time iterations, ti, with a
  //    time-step dt).
  bool last_step = false;
  for (int ti = 1; !last_step; ti++) {
    auto dt_real = std::min(dt, t_final - t);

    ode_solver->Step(phi, t, dt_real);

    last_step = (t >= t_final - 1e-8 * dt);

    // Save at step ti
    paraview_dc.SetCycle(ti);
    paraview_dc.SetTime(t);
    paraview_dc.Save();
  }
  // 10. Free the used memory.
  delete ode_solver;
  delete mesh;
  delete paraview_dc;

  return 0;
}
///////////////////////////////////////////////////////////
/// TODO mettre dans un autre fichier avec des options
/// pour choisir entre heaviside, uniform....
///

double InitialOrderParameter(const Vector &x) {
  if (x.Norml2() < 0.5 * radius) {
    return 0.0;
  } else {
    return 1.0;
  }
}

// TODO changer pour la tangente hyperbolique
double ExactSolution(const Vector &x) {
  if (x.Norml2() < 0.5 * radius) {
    return 0.0;
  } else {
    return 1.0;
  }
}
