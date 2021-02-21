# Quadratic programming solver

## Description
Solve quadratic programming problem with linear constraints.
![Equation](./qp_equation.svg)

Same calling convention as Matlab counterpart.
Equality constraints and boundary constraints are already bake in to the solver.

## Details
The C version of this quadprog algorithm was written back in October, 2019, seem to be fun if it go open source.

The algorithm use the nonlinear optimizer [SolvOpt](https://imsc.uni-graz.at/kuntsevich/solvopt/index.html), it implement Shor's r-algorithm.

Solving a QP problem with nonlinear optimizer would probably overkill, it will just works anyway.

## Accuracy
I'm not sure, but according to my real world application which uses current quadprog library as digital filter design, the final computed result of my real world application does not show accuracy is a major concern.