# spmchf
A B-spline based MCHF code using matrix methods for solving eigenvalue problems for bound states.

Charlotte Froese Fischer
September, 2010

This code is "a work-in-porgress".  A this point, the effort has been to
obtain convergence and maintain accuracy, the emphasis being accuracy.
This is particularly evident in the calculation for Ac where the full
NR method is overkill except that there are lagrange multipliers.

Further effort is needed for efficiency.

To compile and test
```bash
cd src
make
cp spmchf ../test
cd ../test
```
At this point, a number of different tests can be selected.
   - `test/HF` tests the basic HF cases.
   - `test/He`, `test/He_cas`, and `test/He_cas2` test correlation studies for He 1S.

