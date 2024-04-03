# 1D Allen-Cahn simulation  


## Statement of the problem

$$\frac{\partial \phi}{\partial t}=\nabla \cdot{} \lambda \nabla \phi\text{ in }\Omega$$
$$\bf{n} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma$$

## Input file description


```CPP
  const auto DIM = 1;
  using NLFI = AllenCahnNLFormIntegrator<ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
  using FECollection = mfem::H1_FECollection;
  using PSTCollection = mfem::ParaViewDataCollection;
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  using VAR = Variables<FECollection, DIM>;
  using OPE = PhaseFieldOperator<FECollection, DIM, NLFI>;
  using TIME = TimeDiscretization<PST, OPE, VAR>;


```


## Running 

```SHELL
AllenCahn1Dtest1
```

## Post-processing

(to be written)

## Files & Dependencies


(to be written)

## References


(to be written)

## Intellectual Property

(to be written)