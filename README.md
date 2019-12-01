# spmchf
A B-spline based MCHF code using matrix methods for solving eigenvalue problems for bound states.

Charlotte Froese Fischer
September, 2010

This code is "a work-in-porgress".  A this point, the effort has been to
obtain convergence and maintain accuracy, the emphasis being accuracy.
This is particularly evident in the calculation for Ac where the full
NR method is overkill except that there are lagrange multipliers.

Further effort is needed for efficiency.

### Installation
   - Go to the `src/` directory, open Makefile and, if needed, edit compiler: `FC`, compilation options: `FC_FLAGS` and if needed the location of your Lapack and Blas installation: `LAPACK_DIR` (typically not necessary using modern GNU/Linux/Mac systems).
Default is `gfortran` with the debug flag, `g`.
```bash
FC = gfortran
FC_FLAGS = -g
#LAPACK_DIR = -L <path_to_lapack and blas libs>
```
   - Run `make` to compile the codes, which should produce the executable `spmchf` in the `/bin` folder. Add this folder to your (bash) `PATH` variable is you want to.
   - To clean up the compilation files, just run `make clean` as usual.

### Tests
At this point, a number of different tests can be selected.
   - `test/HF` tests the basic HF cases.
   - `test/He`, `test/He_cas`, and `test/He_cas2` test correlation studies for He 1S.
