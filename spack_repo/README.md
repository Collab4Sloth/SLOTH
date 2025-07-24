# Spack repo for SLOTH

## how to use it ?

```
spack env create Sloth
spack env activate Sloth
spack repo add spack_repo
spack install sloth
```

See info:

```
spack info sloth
```

Variant (only petsc):

```
spack install sloth~petsc 
```

Use another compiler

```
spack install sloth%gcc@11.4.0
```
