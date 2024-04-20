#include <iostream>
#include <cmath>

#include "Indexing.h"
#include "Chemistry.h"
#include "Solver_main.h"

double surf_height(const double x){
    return fmin(2.0, 1 + 4*(x-0.5)*(x-0.5));
}

double area_slope(const double x){
    if (x > 1.0){
        return 0.0;
    } else {
        return (2*M_PI*surf_height(x))*8*(x-0.5);
    }
}

void restart(int nelem, double* u){
    FILE* frst = fopen("restart.dat","r"); //restart_therm_equlib.dat | restart_good.dat
    if (frst==nullptr) {
        printf("No restart file found.\n");
        return;
    }

    fscanf(frst,"%*[^\n]\n"); //skip header
    for (int ielem=0; ielem<nelem; ielem++) {
        int id = uIJK(ielem, 0, 0);
        fscanf(frst, "%*lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf%*[^\n]\n", &(u[id]), &(u[id+1]), &(u[id+2]), &(u[id+3]),
               &(u[id+4]), &(u[id+5]), &(u[id+6]), &(u[id+7]));
    }

    fclose(frst);

}

int main() {
    //Inputs
    int nelem, isource;
    double CFL;
    double pb_ratio, pb;

    nelem = 150;
    CFL = 0.4;    // will be reduced by multiplication with CFLTCNE(=0.25) when source terms are enabled
    pb_ratio = 50.0;
    isource = 0;

    ///MAKE SURE THAT THERE IS CONSISTENCE AMONG THE INDEXING FOR DIFFEREENT FACE-VALUED VARIABLES/ARRAYS
    Chem air = Chem();
    double M0, p0, T0, yO20, yN20, yRP0, Yv0, dp0;
    //Initial state as given
    M0 = 2.5;
    p0 = 102325.0;
    T0 = 300;
    yN20 = 0.7643;
    yRP0 = 0.0;
    Yv0 = 0.94;
    dp0 = 10.0e-6;

    double dx, xmin, xmax, xcc[nelem], xfa[nelem+1],
            Acc[nelem],Afa[nelem+1], dAdx[nelem];

    //Define and calculate geometric values
    xmin = 0.0;
    xmax = 1.5;
    dx = (xmax - xmin) / nelem;
    for (int i=0; i<nelem; i++) {
        xcc[aIJ(i,0)] = i*dx + 0.5*dx; //cell centers
        xfa[i] = i*dx; //cell faces

        double s;
        //Area at cell ceneters
        s = surf_height(xcc[aIJ(i,0)]);
        Acc[aIJ(i,0)] = M_PI*s*s;

        //Area at faces of cel
        s = surf_height(xfa[i]);
        Afa[i] = M_PI*s*s;

        //Slope of area at cell centers
        dAdx[aIJ(i,0)] = area_slope(xcc[aIJ(i,0)]);
    }
    xfa[nelem] = nelem*dx;
    Afa[nelem] = M_PI*surf_height(xfa[nelem])*surf_height(xfa[nelem]);

    //Find initial flow variables
    yO20 = 1.0 - yN20 - yRP0;
    double rhol = air.RHOLREF + (p0 - air.PREF)/(air.AREF*air.AREF);
    double rhov = p0 / (T0*(yO20*air.Rs[0] + yN20*air.Rs[1] + yRP0*air.Rs[2]));
    double rho_mix = 1.0 / ((Yv0/rhov) + ((1.0-Yv0)/rhol));

    double rhoR = p0/T0;
    double rhoCv = rhov * (yO20*(air.Calc_cp_Curve(0,T0)-air.Rs[0]) +
                                   yN20*(air.Calc_cp_Curve(1,T0)-air.Rs[1]) +
                                   yRP0*(air.Calc_cp_Curve(2,T0)-air.Rs[2]));
    double a0 = sqrt((p0/rho_mix)*(1.0 + rhoR/rhoCv));

    double u0[8], u0back[8];
    u0[0] = yO20;       //fraction of vapor which is O2
    u0[1] = yN20;       // same for N2
    u0[2] = Yv0;        //Fraction of mas which is vapor
    u0[3] = p0;         //pressure
    u0[4] = M0 * a0;    //velocity
    u0[5] = T0;         //Temperature
    u0[6] = 6.0*(1.0-Yv0) / (M_PI*rhol*dp0*dp0*dp0);         //Number density


    pb = p0*pb_ratio;

    //Approximate post shock conditions
    double p_shock = 7.125;
    double T_shock = 2.1;
    double u_shock = 5.0 * sqrt(T_shock);
    u0back[0] = u0[0];
    u0back[1] = u0[1];
    u0back[2] = u0[2];
    u0back[3] = u0[3]*p_shock;
    u0back[4] = u0[4]*u_shock;
    u0back[5] = u0[5]*T_shock;
    u0back[6] = u0[6];

    //Intialize flow
    auto u = (double*)malloc(nelem*NDEGR*(NSP+3)*sizeof(double));
    //restart(nelem,u);
    for(int ielem=0; ielem<nelem; ielem++) {
        for (int jdegr=0; jdegr<NDEGR; jdegr++){
            for (int kvar=0; kvar<(NSP+3); kvar++){

                if (xfa[ielem]<0.8) {   //Freestream conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0[kvar];

                } else {                //Post-Shock conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0back[kvar];
                }
            }
        }
    }
    //restart(nelem,u);

    int success = solve(isource, nelem, dx, CFL, pb, air, u0, u, xcc, Acc, Afa, dAdx);



    if (!success){
        printf("WRONG! Try again but don't code any bugs this time idiot...\n");
    }

    return 0;
}
