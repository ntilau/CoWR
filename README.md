## Compute Weighted Residuals a.k.a. CoWR
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
Download [x86_64-linux-gnu](https://github.com/ntilau/CoWR/raw/master/bin/x86_64-linux-gnu/core) build or type in your terminal:
```shell
wget https://github.com/ntilau/CoWR/raw/master/bin/x86_64-linux-gnu/core
```

# References
Zhizhang Chen and Michel Ney, [The Method of Weighted Residuals: A General Approach to Deriving Time- and Frequency-Domain Numerical Methods](https://www.cs.buap.mx/~sandoval/MetodosNumericos/MetNumCEROS/MarcoAntonioRamires.pdf), IEEE Antennas and Propagation Magazine, Vol. 51, No.1, February 2009
