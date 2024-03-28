//
// Created by tskoepli on 3/27/2024.
//

#ifndef PROJECT2_CHEMSOURCE_H
#define PROJECT2_CHEMSOURCE_H
#include "Chemistry.h"
#include "StateVariables.h"

void CalcOmega(const double* unk, Chem& air, State& var, double* omega);

#endif //PROJECT2_CHEMSOURCE_H
