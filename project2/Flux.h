//
// Created by Tsail on 3/10/2024.
//

#ifndef PROJECT2_FLUX_H
#define PROJECT2_FLUX_H

#include "Indexing.h"
#include "Chemistry.h"

//void LDFSS(double A, double* uL, double* uR, Chem &air, double* flux);
double wavespeed(const double* u, Chem& air);
void CalcRes(int nelem, double dx,double CFL, Chem &air, double* u0, double* u,double* Acc,double* Afa,double* dAdx, double* res);

#endif //PROJECT2_FLUX_H
