#include <iostream>
#include <cmath>

#include "Indexing.h"
#include "Chemistry.h"
#include "Solver_main.h"
#include "StateVariables.h"

double surf_height(const double x){
    return fmin(2.0, 1.0 + 4*(x-0.5)*(x-0.5));
}

double area_slope(const double x){
    if (x > 1.0){
        return 0.0;
    } else {
        return (2*M_PI*surf_height(x))*8*(x-0.5);
    }
}

void restart(int nelem, double* u){
    FILE* frst = fopen("restart_decoupled.dat","r"); //restart_therm_equlib.dat | restart_good.dat
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

    nelem = 300;
    CFL = 0.05;//0.4;
    pb_ratio = 10.0;
    isource = 1;

    ///MAKE SURE THAT THERE IS CONSISTENCE AMONG THE INDEXING FOR DIFFEREENT FACE-VALUED VARIABLES/ARRAYS
    Chem air = Chem();
    double M0, p0, T0, yO20, yN20, yRP0, Yv0, dp0, deltau0;
    //Initial state as given
    M0 = 2.5;
    p0 = 102325.0;
    T0 = 300.0;
    yN20 = 0.7643;
    yRP0 = 0.0;
    Yv0 = 0.94;
    dp0 = 20.0e-6;
    deltau0 = 0.0;

    ///Mach number higher than prompt   |||   Temperature lower than prompt

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
        Acc[aIJ(i,0)] = s   ;//  M_PI*s*s;

        //Area at faces of cel
        s = surf_height(xfa[i]);
        Afa[i] = s      ;// M_PI*s*s;

        //Slope of area at cell centers
        dAdx[aIJ(i,0)] = area_slope(xcc[aIJ(i,0)]);
    }
    xfa[nelem] = nelem*dx;
    Afa[nelem] = surf_height(xfa[nelem])    ;//  M_PI*surf_height(xfa[nelem])*surf_height(xfa[nelem]);

    //Find initial flow variables
    yO20 = 1.0 - yN20 - yRP0;
    double rhol = air.RHOLREF + (p0 - air.PREF)/(air.AREF*air.AREF);

    double u0[8], u0back[8];
    u0[0] = yO20;       //fraction of vapor which is O2
    u0[1] = yN20;       // same for N2
    u0[2] = Yv0;        //Fraction of mas which is vapor
    u0[3] = p0;         //pressure
    u0[4] = M0 * 300.0; //velocity placeholder
    u0[5] = T0;         //Temperature
    u0[6] = 6.0*(1.0-Yv0) / (M_PI*rhol*dp0*dp0*dp0);         //Number density
    u0[7] = deltau0;    // initial drift velocity

    State var;
    var.Initialize(u0);
    var.UpdateState(air,isource);
    u0[4] = M0 * var.a;    //actual velocity
    var.UpdateState(air,isource);

    pb = p0*pb_ratio;

    //Approximate post shock conditions
    double p_shock = 7.125;
    double T_shock = 2.1;
    double u_shock = 5.0 * sqrt(T_shock);
    u0back[0] = u0[0];
    u0back[1] = u0[1];
    u0back[2] = u0[2];
    u0back[3] = u0[3]*p_shock;
    u0back[4] = u0[4]/u_shock;
    u0back[5] = u0[5]*T_shock;
    u0back[6] = u0[6];
    u0back[7] = u0[7];

    //Intialize flow
    auto u = (double*)malloc(nelem*NDEGR*(NVAR)*sizeof(double));
    //restart(nelem,u);
    for(int ielem=0; ielem<nelem; ielem++) {
        for (int jdegr=0; jdegr<NDEGR; jdegr++){
            for (int kvar=0; kvar<(NVAR); kvar++){

                if (xfa[ielem]<0.8) {   //Freestream conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0[kvar];

                } else {                //Post-Shock conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0back[kvar];
                }

                if (kvar==7){
                    //add intial drif velocity for when testing non-pressure terms
                    u[uIJK(ielem, jdegr, kvar)] = 1e-2;
                }
            }
        }
    }
    restart(nelem,u);

    int success = solve(isource, nelem, dx, CFL, pb, air, u0, u, xcc, Acc, Afa, dAdx);



    if (!success){
        printf("WRONG! Try again but do it right this time idiot...\n");
    }

    return 0;
}
