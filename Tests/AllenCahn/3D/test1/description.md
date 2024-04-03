# 3D Allen-Cahn simulation in a square 


## Statement of the problem

$$\frac{\partial \phi}{\partial t}=\nabla \cdot{} \lambda \nabla \phi\text{ in }\Omega$$
$$\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{top}$$
$$\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{bottom}$$
$$\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{front}$$
$$\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{rear}$$
$$\phi=0 \text{ on }\Gamma_{left}$$
$$\phi=1 \text{ on }\Gamma_{right}$$

(to be finished)
## Input file description


```CPP
  const auto DIM = 3;
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
AllenCahn3Dtest1
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