//
// Created by Tsail on 3/10/2024.
//

#include "Solver_main.h"
#include "Indexing.h"
#include "Flux.h"


void solve_nonreacting(int nelem, int nsp, double dx,double CFL,double* u,double* Acc,double* Afa,double* dAdx) {
    //Solve the nonreacting / chemically frozen problem
    for (int iter=0; iter<MXITER; iter++){
        //Loop through elements and conduct local timestepping

        //Find the common flux at each face
        double flux_comm[(nelem+1)*(nsp+3)];
        //Interior faces
        for (int iface=1; iface<nelem; iface++){
            double *uL, *uR, *flx;
            uL = &(u[uIJK(iface-1,0,0)]);
            uR = &(u[uIJK(iface,0,0)]);
            flx = &flux_comm[aIJ(iface,0)];

            LDFSS(uL, uR, flx);
        }

        //Calculate residual RHS and diagonal element sof implicit operator
        for (int ielem=0; ielem<nelem; ielem++) {

        }
    }

}