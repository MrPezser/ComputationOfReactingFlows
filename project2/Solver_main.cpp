//
// Created by Tsail on 3/10/2024.
//

#include "Solver_main.h"
#include "Indexing.h"
#include "Flux.h"


void solve_nonreacting(int nelem, double dx,double CFL, Chem &air, double* u,double* Acc,double* Afa,double* dAdx) {
    //Solve the nonreacting / chemically frozen problem
    double res[(nelem)*(NSP+3)];

    for (int iter=0; iter<MXITER; iter++){
        //Loop through elements and conduct local timestepping

        //Calculate residual
        CalcRes(nelem, dx,CFL, air, u,Acc,Afa,dAdx, res);

    }


}