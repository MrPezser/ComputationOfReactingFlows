//
// Created by tskoepli on 3/20/2024.
//

#ifndef PROJECT2_LINEARSOLVE_H
#define PROJECT2_LINEARSOLVE_H

int LUPDecompose(double **A, int N, double Tol, int *P);
void LUPSolve(double **A, int *P, double *b, int N, double *x);


#endif //PROJECT2_LINEARSOLVE_H
