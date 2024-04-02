#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-loop-convert"
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

void ResNorm(int nelem, const double* res, double* resout) {
    // find an "unweighted" L2 norm, since dx is constant it will get cancelled out when the res is normalized by init
    double ressum[NSP + 3]{0.0};
    for (int ie = 0; ie < nelem; ie++) {
        for (int iv = 0; iv < NSP + 3; iv++) {
            ressum[iv] += (res[fIJ(ie, iv)]) * (res[fIJ(ie, iv)]);
        }
    }

    for (int ivar=0; ivar<NSP+3; ivar++){
        resout[ivar] = sqrt(ressum[ivar]);
    }
}


int IterUpdate(int& ireact, int iter, int nelem, const double* res, double* res0){
    printf("Iter:%8d",iter);

    double ressum[NSP+3]{0.0}, resnorm{0.0};
    ResNorm(nelem, res, ressum);

    if (ireact==1 && res0[NSP+2]<1e-16) {
        res0[NSP+2] = ressum[NSP+2];
    }

    for (int ivar=0; ivar<NSP+2+ireact; ivar++) {
        resnorm +=  ressum[ivar]/ res0[ivar];
    }

    printf("\tRes/Res0:%8.1e", resnorm);

    if (IDETAIL)
    {
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
        printf("\t\tRes:%8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e ",
               ressum[0],ressum[1],ressum[2],ressum[3],ressum[4],ressum[5],ressum[6],ressum[7]);
    }
    printf("\n");

    if (resnorm < RESTOL){
        return 1;
    } else{
        if (resnorm < RXTOL && ireact==0) {
            ireact = 1;
            printf("ENABLING THERMOCHEMICAL SOURCE TERMS\n");
        }
        return 0;
    }
}

int solve(int& ireact, int nelem, double dx, double CFL, double pb, Chem &air, double* u0, double* u, double* xcc,
                      const double* Acc,const double* Afa,const double* dAdx) {
    //Solve the nonreacting / chemically frozen problem

    //Create some arrays
    int iconv = 0;
    double res[(nelem)*(NSP+3)*NDEGR], dv[(nelem)*(NSP+3)*NDEGR];
    double res0[NSP+3] ;
    double dt;
    auto D = (double**)malloc((NSP+3) * sizeof(double*)); //Same memory to be used for each local matrix (chg this if making parallel)
    for (int isp = 0; isp < NSP+3; isp++)
        D[isp] = (double*)malloc( (NSP+3) * sizeof(double));
    //FILE* fres = fopen("reshist.tec", "w");
    State ElemVar[nelem+1];

    //Set up structures for calculating/containing non-state variables on each element
    for (int ielem=0; ielem<nelem; ielem++){
        int id = uIJK(ielem,0,0);
        ElemVar[ielem].Initialize(&(u[id]));
        ElemVar[ielem].UpdateState(air, ireact);
    }
    ElemVar[nelem].Initialize(u0);
    ElemVar[nelem].UpdateState(air, ireact);

    for (int iter=0; iter<MXITER; iter++){
        int flg{};
        //Loop through elements and conduct local timestepping
        //========== Calculate Flow Properties in each element
        for (int iel=0; iel<nelem; iel++) {
            ElemVar[iel].UpdateState(air, ireact);
        }

        //========== Calculate residual
        CalcRes(ireact, nelem, dx, CFL, pb, air, ElemVar, u0, u, Acc, Afa, dAdx, res);

        if (iter==0){
            ResNorm(nelem, res, res0);
            double rhores = res0[0] + res0[1] + res0[2] + res0[3] + res0[4];
            for (int isp=0; isp<NSP; isp++){
                res0[isp] = rhores;
            }
        }

        //========== Solve linear system on each element
        for (int ielem=0; ielem<nelem; ielem++) {
            double* unk = &(u[uIJK(ielem,0,0)]);
            double tol = 1e-32;
            int P[NSP+3]{}; //permutation vector for pivoting
            int N = NSP+2+ireact;
            dt = dx*CFL/(ElemVar[ielem].a + fabs(unk[NSP]));

            //Evaluate the jacobian / Implicit matrix
            BuildJacobian(ireact, dt, unk, air, ElemVar[ielem], D);

            for (int i=0; i<NSP+3; i++) {
                for (int j=0; j<NSP+3; j++) {
                    D[i][j] *= Acc[ielem];
                }
            }
            //get the rhs block needed
            double* b = &(res[uIJK(ielem,0,0)]);
            double* x = &(dv[uIJK(ielem,0,0)]);

            flg = LUPDecompose(D, N, tol, P);
            if (flg >= 0){printf("LU Decomp Fail on row %d\n", flg+1);break;}
            LUPSolve(D, P, b, N, x);
        }
        if (flg >= 0){return 0;}




        //Carry out Euler Timestep
        for (int ielem=0; ielem<nelem; ielem++){
            double* ui = &(u[uIJK(ielem,0,0)]);

            for (int kvar=0; kvar<NSP+3; kvar++){
                int id = uIJK(ielem,0,kvar);
                u[id] += dv[id];
                dv[id] = 0.0; //reset
            }

            //temperature limiting
            ui[NSP+2] = fmax(ui[NSP+2], 201);
            ui[NSP+1] = fmax(ui[NSP+1], 201);

            ui[NSP+2] = fmin(ui[NSP+2], 9500);
            ui[NSP+1] = fmin(ui[NSP+1], 9500);

            //frozen thermo
            if (ireact==0){
                ui[NSP+2] = u0[NSP+1];
            }
        }



        if (iter%1000 ==0) {
            iconv = IterUpdate(ireact, iter, nelem, res, res0);


            //save soln file
            FILE* fout = fopen("waveout.tec", "w");
            if (fout == nullptr) {
                printf("~~~~~~~~~ Failed to save output file, error:%d\n", errno);}
            else {
                fprintf(fout, "TITLE = \"%s\"\n", "title");
                fprintf(fout, "VARIABLES = \"X\",\"rhoN2\",\"rhoO2\",\"rhoNO\",\"rhoN\",\"rhoO\","
                              "\"u\",\"T\",\"Tv\",\"P\",\"M\",\"YN2\",\"YO2\",\"YNO\",\"YN\",\"YO\"\n");
                fprintf(fout, "ZONE I=%d, DATAPACKING=POINT\n", nelem);

                for (int i=0; i<nelem; i++) {
                    double rm = ElemVar[i].rho_mix;
                    double *r = &(u[uIJK(i,0,0)]);
                    fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",xcc[i],
                            r[0], r[1], r[2], r[3], r[4], u[uIJK(i,0,5)], u[uIJK(i,0,6)], u[uIJK(i,0,7)],
                            ElemVar[i].p, u[uIJK(i,0,5)]/ElemVar[i].a, r[0]/rm, r[1]/rm, r[2]/rm, r[3]/rm, r[4]/rm);
                }
            }
            fclose(fout);


        }
        if (iconv == 1) {break;}
    }

    free(D);

    /*
    //save soln file
    FILE* fout = fopen("waveout.tec", "w");
    if (fout == nullptr) {
        printf("~~~~~~~~~ Failed to save output file, error:%d\n", errno);}
    else {
        fprintf(fout, "TITLE = \"%s\"\n", "title");
        fprintf(fout, "VARIABLES = \"X\",\"rhoN2\",\"rhoO2\",\"rhoNO\",\"rhoN\",\"rhoO\","
                      "\"u\",\"T\",\"Tv\",\"P\",\"M\",\"YN2\",\"YO2\",\"YNO\",\"YN\",\"YO\"\n");
        fprintf(fout, "ZONE I=%d, DATAPACKING=POINT\n", nelem);

        for (int i=0; i<nelem; i++) {
            double rm = ElemVar[i].rho_mix;
            double *r = &(u[uIJK(i,0,0)]);
            fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",xcc[i],
                    r[0], r[1], r[2], r[3], r[4], u[uIJK(i,0,5)], u[uIJK(i,0,6)], u[uIJK(i,0,7)],
                    ElemVar[i].p, u[uIJK(i,0,5)]/ElemVar[i].a, r[0]/rm, r[1]/rm, r[2]/rm, r[3]/rm, r[4]/rm);
        }
    }
    fclose(fout);
     */
    return 1;
}
#pragma clang diagnostic pop