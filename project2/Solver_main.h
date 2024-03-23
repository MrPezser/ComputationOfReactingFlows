//
// Created by Tsail on 3/10/2024.
//

#ifndef PROJECT2_SOLVER_MAIN_H
#define PROJECT2_SOLVER_MAIN_H

#include "Chemistry.h"

#define MXITER 1e5
#define RESTOL 5e-5
void solve_nonreacting(int nelem, double dx,double CFL, Chem &air, double* u0, double* u, double* xcc,double* Acc,double* Afa,double* dAdx);

#endif //PROJECT2_SOLVER_MAIN_H
