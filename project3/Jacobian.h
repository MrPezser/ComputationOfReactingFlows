//
// Created by tskoepli on 3/18/2024.
//

#ifndef PROJECT2_JACOBIAN_H
#define PROJECT2_JACOBIAN_H

#include "Chemistry.h"
#include "Indexing.h"
#include "StateVariables.h"

void BuildJacobian(int isource, double dt, const double* unk, Chem &air, State& var, double** D);

#endif //PROJECT2_JACOBIAN_H
