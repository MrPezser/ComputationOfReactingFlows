#include <iostream>
#include <cmath>

#include "Indexing.h"

double surf_height(const double x){
    return fmin(2.0, 1 + 4*(x-0.5)*(x-0.5));
}

double area_slope(const double x){
    if (x > 1.0){
        return 0.0;
    } else {
        return 8*(x-0.5);
    }
}

int main() {
    //Inputs
    int nelem;
    nelem = 200;

    ///MAKE SURE THAT THERE IS CONSISTENCE AMONG THE INDEXING FOR DIFFEREENT FACE-VALUED VARIABLES/ARRAYS

    double vel0, rho0, T0, pb, YN20, YNO0, YO0, YN0, YO20;
    //Initial state as given
    vel0 = 3500.0;
    rho0 = 0.0003074;
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
    u0[0] = YO0;
    u0[1] = YN0;
    u0[2] = YO20;
    u0[3] = YN20;
    u0[4] = YNO0;
    u0[5] = 3500.0; //u_freestream
    u0[6] = 350.0;  //trns-rot temperature
    u0[6] = 350.0;  //vib temperautre

    //Calculate post shock conditions (approx Mach 9.33333 normal shock)
    u0back[0] = u0[0];
    u0back[1] = u0[1];
    u0back[2] = u0[2];
    u0back[3] = u0[3];
    u0back[4] = u0[4];
    u0back[5] = 616.808;
    u0back[6] = 6258.0;
    u0back[7] = 6258.0;

    //Intialize flow
    auto u = (double*)malloc(nelem*NDEGR*(NSP+3)*sizeof(double));
    for(int ielem=0; ielem<nelem; ielem++) {
        for (int jdegr=0; jdegr<NDEGR; jdegr++){
            for (int kvar=0; kvar<(NSP+2); kvar++){

                if (xfa[ielem]<1.0) {   //Freestream conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0[kvar];

                } else {                //Post-Shock conditions
                    u[uIJK(ielem, jdegr, kvar)] = u0back[kvar];
                }
            }
        }
    }

    solve_nonreacting();




    return 0;
}
