#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-loop-convert"
//
// Created by Tsail on 3/10/2024.
//

#include <errno.h>

#include "Solver_main.h"
#include "Indexing.h"
#include "RHS.h"
#include "Jacobian.h"
#include "LinearSolve.h"
#include "StateVariables.h"

void ResNorm(int nelem, const double* res, double* resout) {
    // find an "unweighted" L2 norm, since dx is constant it will get cancelled out when the res is normalized by init
    double ressum[NVAR]{0.0};
    for (int ie = 0; ie < nelem; ie++) {
        for (int iv = 0; iv < NVAR; iv++) {
            ressum[iv] += (res[fIJ(ie, iv)]) * (res[fIJ(ie, iv)]);
        }
    }

    for (int ivar=0; ivar<NVAR; ivar++){
        resout[ivar] = sqrt(ressum[ivar]);
    }
}


int IterUpdate(int& isource, int iter, int nelem,double& CFL, const double* res, double* res0){
    printf("Iter:%8d",iter);

    double ressum[NVAR]{0.0}, resnorm{0.0};
    ResNorm(nelem, res, ressum);

    /* potential ignoring of multiphase source term
    if (isource==1 && res0[2]<1e-16) {
        res0[2] = ressum[2];
    }*/

    for (int ivar=0; ivar<NVAR; ivar++) {
        resnorm +=  ressum[ivar]/ res0[ivar];
    }

    printf("\tRes/Res0:%8.1e", resnorm);

    if (IDETAIL)
    {
        for (int ie=0; ie<nelem; ie++){
            for (int iv=0; iv<NVAR; iv++) {
                ressum[iv] += (res[fIJ(ie,iv)])*(res[fIJ(ie,iv)]);
            }
        }
        double maxres = sqrt(ressum[0]);
        for (int iv=0; iv<NVAR; iv++) { // NOLINT(modernize-loop-convert)
            ressum[iv] = sqrt(ressum[iv]);
            maxres = fmax(maxres, ressum[iv]);
        }
        printf("\t\tRes:%8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e ",
               ressum[0],ressum[1],ressum[2],ressum[3],ressum[4],ressum[5],ressum[6]);
    }
    printf("\n");

    if (resnorm < RESTOL){
        return 1;
    } else{
        /*
        if (resnorm < RXTOL && isource==0) {
            isource = 1;
            CFL = CFL*CFLTCNE;
            //return 1;///ONLY DO FULLY FROZEN FLOW
            printf("ENABLING MULTIPHASE SOURCE TERMS\n");
        }
         */
        return 0;
    }
}

int solve(int& isource, int nelem, double dx, double CFL, double pb, Chem &air, double* u0, double* u, double* xcc,
                      const double* Acc,const double* Afa,const double* dAdx) {
    //Solve the nonreacting / chemically frozen problem

    //Create some arrays
    int iconv = 0;
    double res[(nelem)*(NVAR)*NDEGR], dv[(nelem)*(NVAR)*NDEGR];
    double res0[NVAR] ;
    double dt;
    auto D = (double**)malloc((NVAR) * sizeof(double*)); //Same memory to be used for each local matrix (chg this if making parallel)
    for (int isp = 0; isp < NVAR; isp++)
        D[isp] = (double*)malloc( (NVAR) * sizeof(double));
    //FILE* fres = fopen("reshist.tec", "w");
    State ElemVar[nelem+1];

    //Set up structures for calculating/containing non-state variables on each element
    for (int ielem=0; ielem<nelem; ielem++){
        int id = uIJK(ielem,0,0);
        ElemVar[ielem].Initialize(&(u[id]));
        ElemVar[ielem].UpdateState(air, isource);
    }
    ElemVar[nelem].Initialize(u0);
    ElemVar[nelem].UpdateState(air, isource);

    //save residual history
    FILE* fres = fopen("res.tec", "w");
    if (fres == nullptr) {
        printf("~~~~~~~~~ Failed to save residual file, error:%d\n", errno);}
    else {
        fprintf(fres, "Residual history\n");
        fprintf(fres, "%d,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n",0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);
    }

    for (int iter=0; iter<MXITER; iter++){
        int flg{};
        //Loop through elements and conduct local timestepping
        //========== Calculate Flow Properties in each element
        for (int iel=0; iel<nelem; iel++) {
            ElemVar[iel].UpdateState(air, isource);
        }

        //========== Calculate residual
        CalcRes(isource, nelem, dx, CFL, pb, air, ElemVar, u0, u, Acc, Afa, dAdx, res);

        if (iter==0){
            ResNorm(nelem, res, res0);

            //normalize the first two residuals by the total change in vapor species fraction
            //double rhores = res0[0] + res0[1];
            //for (int isp=0; isp<2; isp++){
            //    res0[isp] = rhores;
            //}
        }

        //========== Solve linear system on each element
        for (int ielem=0; ielem<nelem; ielem++) {
            double* unk = &(u[uIJK(ielem,0,0)]);
            double tol = 1e-20;
            int P[NVAR+1]{}; //permutation vector for pivoting
            int N = NVAR;
            dt = dx*CFL/(ElemVar[ielem].a + fabs(unk[4]));

            //Evaluate the jacobian / Implicit matrix
            BuildJacobian(isource, dt, unk, air, ElemVar[ielem], D);

            for (int i=0; i<NVAR; i++) {
                for (int j=0; j<NVAR; j++) {
                    D[i][j] *= Acc[ielem];
                }
            }
            //get the rhs block needed
            double* b = &(res[uIJK(ielem,0,0)]);
            double* x = &(dv[uIJK(ielem,0,0)]);

            (void)LUPDecompose(D, N, tol, P);
            LUPSolve(D, P, b, N, x);
        }




        //Carry out Timestep
        for (int ielem=0; ielem<nelem; ielem++){
            double* ui = &(u[uIJK(ielem,0,0)]);

            for (int kvar=0; kvar<NVAR; kvar++){
                int id = uIJK(ielem,0,kvar);
                u[id] += dv[id];
                dv[id] = 0.0; //reset dv array
            }

            //temperature limiting
            ui[5] = fmax(ui[5], 201);
            ui[5] = fmin(ui[5], 6001);

        }

        //Save to residual log file
        double ressum[NVAR];
        ResNorm(nelem, res, ressum);
        if ( res0[2] > 1e-16) {
            fprintf(fres, "%d,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", iter+1, ressum[0] / res0[0], ressum[1] / res0[1],
                    ressum[2] / res0[2], ressum[3] / res0[3], ressum[4] / res0[4], ressum[5] / res0[5], ressum[6] / res0[6]);
        } else {
            fprintf(fres, "%d,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", iter+1, ressum[0] / res0[0], ressum[1] / res0[1],
                    1.0, ressum[3] / res0[3], ressum[4] / res0[4], ressum[5] / res0[5], ressum[6] / res0[6]);
        }
        //Printout Solution and residual
        if (iter%100 == 0) {
            iconv = IterUpdate(isource, iter, nelem, CFL, res, res0);

            //save soln file
            FILE* fout = fopen("waveout.tec", "w");
            if (fout == nullptr) {
                printf("~~~~~~~~~ Failed to save output file, error:%d\n", errno);}
            else {
                fprintf(fout, "TITLE = \"%s\"\n", "title");
                fprintf(fout, "VARIABLES = \"X\",\"yO2v\",\"yN2v\",\"Yv\",\"P\",\"u\","
                              "\"T\",\"n_tilde\",\"rho_mix\",\"Mach\"\n");
                fprintf(fout, "ZONE I=%d, DATAPACKING=POINT\n", nelem);

                for (int i=0; i<nelem; i++) {
                    double rm = ElemVar[i].rho_mix;
                    double* unk = &(u[uIJK(i,0,0)]);
                    fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",xcc[i],
                            unk[0], unk[1], unk[2], unk[3], unk[4], unk[5], unk[6],rm, unk[4]/ElemVar[i].a);
                }
            }
            fclose(fout);
        }

        if (iconv == 1) {break;}
    }

    free(D);
    fclose(fres);
    return 1;
}
#pragma clang diagnostic pop