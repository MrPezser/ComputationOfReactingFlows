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
    FILE* frst = fopen("restart_therm_equlib.dat","r"); //restart_therm_equlib.dat | restart_good.dat
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
    int nelem, ireact;
    double CFL;
    double pb_ratio, pb,p0{};

    nelem = 200;
    CFL = 0.4;    // will be reduced by multiplication with CFLTCNE when source terms are enabled
    pb_ratio = 200.0;
    ireact = 0;

    ///MAKE SURE THAT THERE IS CONSISTENCE AMONG THE INDEXING FOR DIFFEREENT FACE-VALUED VARIABLES/ARRAYS
    Chem air = Chem();
    double vel0, rho0, T0, YN20, YNO0, YO0, YN0, YO20;
    //Initial state as given
    vel0 = 3500.0;
    rho0 = 1.0 * 3.074e-4;
    T0 = 350;
    YN20 = 0.7643;
    YNO0 = 0.0;
    YO0 = 0.0;
    YN0 = 0.0;
    YO20 = 1.0 - YN20 - YNO0 - YO0 - YN0;

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
    double u0[8], u0back[8];
    u0[0] = YN20*rho0;
    u0[1] = YO20*rho0;
    u0[2] = YNO0*rho0;
    u0[3] = YN0*rho0;
    u0[4] = YO0*rho0;
    u0[5] = vel0; //u_freestream
    u0[6] = T0;  //trns-rot temperature
    u0[7] = T0;  //vib temperautre

    //
    for (int isp=0;isp<NSP;isp++) {
        p0 += u0[isp] * (air.Ruv/air.Mw[isp]) * T0;
    }
    pb = p0*pb_ratio;

    //Approximate post shock conditions
    double rho_shock = 5.8;
    u0back[0] = u0[0]*rho_shock;
    u0back[1] = u0[1]*rho_shock;
    u0back[2] = u0[2]*rho_shock;
    u0back[3] = u0[3]*rho_shock;
    u0back[4] = u0[4]*rho_shock;
    u0back[5] = u0[5]/rho_shock;
    u0back[6] = u0[6] * 17.0;
    u0back[7] = u0[7] * 17.0;

    //Intialize flow
    auto u = (double*)malloc(nelem*NDEGR*(NSP+3)*sizeof(double));
    //restart(nelem,u);
    for(int ielem=0; ielem<nelem; ielem++) {
        for (int jdegr=0; jdegr<NDEGR; jdegr++){
            //u[uIJK(ielem,0,7)] = u[uIJK(ielem,0,6)];
            for (int kvar=0; kvar<(NSP+3); kvar++){

                if (xfa[ielem]<0.50) {   //Freestream conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0[kvar];

                } else {                //Post-Shock conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0back[kvar];
                }
            }
        }
    }
    //restart(nelem,u);

    int success = solve(ireact, nelem, dx, CFL, pb, air, u0, u, xcc, Acc, Afa, dAdx);



    if (!success){
        printf("WRONG! Try again but don't code any bugs this time idiot...\n");
    }

    return 0;
}
