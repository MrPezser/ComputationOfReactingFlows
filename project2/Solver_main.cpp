//
// Created by Tsail on 3/10/2024.
//

#include "Solver_main.h"
#include "Indexing.h"
#include "Flux.h"
#include "VariableTransform.h"


void solve_nonreacting(int nelem, double dx,double CFL, Chem &air, double* u,double* Acc,double* Afa,double* dAdx) {
    //Solve the nonreacting / chemically frozen problem
    double res[(nelem)*(NSP+3)];

    for (int iter=0; iter<MXITER; iter++){
        //Loop through elements and conduct local timestepping

        //Calculate residual
        CalcRes(nelem, dx,CFL, air, u,Acc,Afa,dAdx, res);

        // Solve linear system on each element
        auto D = (double**)malloc((NSP+3) * sizeof(double*));
        for (int isp = 0; isp < NSP+3; isp++)
            D[isp] = (double*)malloc( (NSP+3) * sizeof(double));s

        for (int ielem=0; ielem<nelem; ielem++) {
            double* unk = &(u[uIJK(ielem,0,0)]);
            BuildDudv(unk, air, D);
        }
    }


}