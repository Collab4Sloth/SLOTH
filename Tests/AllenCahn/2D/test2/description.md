# 2D Allen-Cahn simulation in a fragment  


## Statement of the problem

$$
\frac{\partial \phi}{\partial t}=\nabla \cdot{} \lambda \nabla \phi\text{ in }\Omega
$$

$$
\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma
$$

(to be finished)
## Input file description


```CPP
  const auto DIM = 2;
  using NLFI =
      AllenCahnMeltingNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                       ThermodynamicsPotentials::W, ThermodynamicsPotentials::H,
                                       Mobility::Constant, PhaseChange::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperatorMelting<FECollection, DIM, NLFI>;
  using TIME = TimeDiscretization<PST, OPE, VAR>;

```

(to be finished)

## Running 

```SHELL
AllenCahn2Dtest2
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