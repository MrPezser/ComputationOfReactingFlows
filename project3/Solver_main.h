//
// Created by Tsail on 3/10/2024.
//

#ifndef PROJECT2_SOLVER_MAIN_H
#define PROJECT2_SOLVER_MAIN_H

#include "Chemistry.h"


int solve(int& ireact, int nelem, double dx, double CFL, double pb, Chem &air, double* u0, double* u, double* xcc,
                      const double* Acc,const double* Afa,const double* dAdx);

#endif //PROJECT2_SOLVER_MAIN_H
