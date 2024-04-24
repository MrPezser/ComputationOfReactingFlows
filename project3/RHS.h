//
// Created by Tsail on 3/10/2024.
//

#ifndef PROJECT2_RHS_H
#define PROJECT2_RHS_H

#include "Indexing.h"
#include "Chemistry.h"
#include "StateVariables.h"

void CalcRes(int ireact, int nelem, double dx, double CFL, double pb, Chem &air, State* ElemVar, double* u0, double* u,
             const double* Acc,const double* Afa,const double* dAdx, double* res);

#endif //PROJECT2_RHS_H
