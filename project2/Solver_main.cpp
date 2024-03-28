//
// Created by Tsail on 3/10/2024.
//

#include <errno.h>

#include "Solver_main.h"
#include "Indexing.h"
#include "Flux.h"
#include "Jacobian.h"
#include "LinearSolve.h"
#include "StateVariables.h"

int IterUpdate(int iter, int nelem, const double* res){
    printf("Iter:%5d\t",iter);

    double ressum[NSP+3]{0.0};

    for (int ie=0; ie<nelem; ie++){
        for (int iv=0; iv<NSP+3; iv++) {
            ressum[iv] += (res[fIJ(ie,iv)])*(res[fIJ(ie,iv)]);
        }
    }
    double maxres = sqrt(ressum[0]);
    for (int iv=0; iv<NSP+3; iv++) { // NOLINT(modernize-loop-convert)
        ressum[iv] = sqrt(ressum[iv]);
        maxres = fmax(maxres, ressum[iv]);
    }
    printf("Res:%8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e ",
           ressum[0],ressum[1],ressum[2],ressum[3],ressum[4],ressum[5],ressum[6],ressum[7]);

    printf("\n");

    if (maxres < RESTOL){
        return 1;
    } else {
        return 0;
    }
}

int solve_nonreacting(int nelem, double dx, double CFL, double pb, Chem &air, double* u0, double* u, double* xcc,
                      const double* Acc,const double* Afa,const double* dAdx) {
    //Solve the nonreacting / chemically frozen problem

    //Create some arrays
    int iconv = 0;
    double res[(nelem)*(NSP+3)*NDEGR];
    double dv[(nelem)*(NSP+3)*NDEGR];
    auto D = (double**)malloc((NSP+3) * sizeof(double*)); //Same memory to be used for each local matrix (chg this if making parallel)
    for (int isp = 0; isp < NSP+3; isp++)
        D[isp] = (double*)malloc( (NSP+3) * sizeof(double));

    State ElemVar[nelem+1];

    //Set up structures for calculating/containing non-state variables on each element
    for (int ielem=0; ielem<nelem; ielem++){
        int id = uIJK(ielem,0,0);
        ElemVar[ielem].Initialize(&(u[id]));
        ElemVar[ielem].UpdateState(air);
    }
    ElemVar[nelem].Initialize(u0);
    ElemVar[nelem].UpdateState(air);


    for (int iter=0; iter<MXITER; iter++){
        int flg{};
        //Loop through elements and conduct local timestepping

        //========== Calculate residual
        CalcRes(nelem, dx, CFL, pb, air, ElemVar, u0, u, Acc, Afa, dAdx, res);

        //========== Solve linear system on each element
        for (int ielem=0; ielem<nelem; ielem++) {
            double* unk = &(u[uIJK(ielem,0,0)]);
            double tol = 1e-12;
            int P[NSP+3]{}; //permutation vector for pivoting
            int N = NSP+3;

            //Evaluate the jacobian / Implicit matrix
            BuildJacobian(unk, air, ElemVar[ielem], D);

            for (int i=0; i<NSP+3; i++) {
                for (int j=0; j<NSP+3; j++) {
                    D[i][j] *= Acc[ielem];
                }
            }
            //get the rhs block needed
            double* b = &(res[uIJK(ielem,0,0)]);
            double* x = &(dv[uIJK(ielem,0,0)]);

            flg = LUPDecompose(D, N, tol, P);
            if (flg == 0){printf("LU Decomp Fail\n");break;}
            LUPSolve(D, P, b, N, x);
        }
        if (flg == 0){return 0;}




        //Carry out Euler Timestep
        for (int ielem=0; ielem<nelem; ielem++){
            double* ui = &(u[uIJK(ielem,0,0)]);
            double a = wavespeed(ui, air);
            double dt = dx*CFL/(a + fabs(ui[NSP]));

            for (int kvar=0; kvar<NSP+3; kvar++){
                int id = uIJK(ielem,0,kvar);
                u[id] += dt*dv[id];
                dv[id] = 0.0; //reset
            }

            ElemVar[ielem].UpdateState(air);
        }



        if (iter%1000 ==0) {
            iconv = IterUpdate(iter, nelem, res);
            //save soln file
            FILE* fout = fopen("waveout.tec", "w");
            if (fout == nullptr)
            {
                printf("~~~~~~~~~ Failed to save output file, error:%d\n", errno);
                //printf("Oh dear, something went wrong with read()! %s\n", strerror(errno));
            } else {
                fprintf(fout, "x\trhoN2\trhoO2\trhoNO\trhoN\trhoO\tu\tT\tTv\n");

                for (int i=0; i<nelem; i++) {
                    fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",xcc[i],
                            u[uIJK(i,0,0)], u[uIJK(i,0,1)], u[uIJK(i,0,2)], u[uIJK(i,0,3)], u[uIJK(i,0,4)], u[uIJK(i,0,5)], u[uIJK(i,0,6)], u[uIJK(i,0,7)]);
                }
            }
            fclose(fout);


        }
        if (iconv == 1) {break;}
    }

    free(D);
    return 1;
}