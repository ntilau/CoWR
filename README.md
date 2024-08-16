## Computing Residuals a.k.a. CORE
Generalized version of FES to extend to any Weighted Residuals Methods framework
The goal is to have an as modular as possible implementation minimizing the need
for external libraries.
- Preprocessing and Postprocessing with VTK file formats
  + Mesh manipulations with TRIANGLE & TETGEN
- Wave Equation solvers in both FE or MOM formulations
  + BLAS, ARPACK & MUMPS solvers to be standardized on
- Low frequency stabilization to be addressed

# Revision history
To-do: 
- Shared libraries conversion
- Pre preprocessing information in modular fashion
- Modular post processing

# Releases
Download [x86_64-linux-gnu](https://github.com/ntilau/core/raw/master/bin/x86_64-linux-gnu/fes) build or type in your terminal:
```shell
wget https://github.com/ntilau/core/raw/master/bin/x86_64-linux-gnu/fes
```
