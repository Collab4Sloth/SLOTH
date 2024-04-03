# 2D Allen-Cahn simulation in a periodic square 


## Statement of the problem

$$
\frac{\partial \phi}{\partial t}=\nabla \cdot{} \lambda \nabla \phi\text{ in }\Omega
$$

$$
\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{top}
$$

$$
\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{bottom}
$$

$$
\phi_{\Gamma_{left}}=\phi_{\Gamma_{right}}
$$

(to be finished)
## Input file description


```CPP
  const auto DIM = 2;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI>;
  using TIME = TimeDiscretization<PST, OPE, VAR>;
```

(to be finished)

## Running 

```SHELL
AllenCahn2DPeriodictest2
```
(to be finished)

## Post-processing

(to be written)

## Files & Dependencies


(to be written)

## References


(to be written)

## Intellectual Property

(to be written)