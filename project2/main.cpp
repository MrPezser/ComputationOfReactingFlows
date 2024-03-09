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
    double vel0, rho0, T0, pb, YN20, YNO0, YO0, YN0;
    //Initial state as given
    vel0 = 3500.0;
    rho0 = 0.0003074;
    T0 = 350;
    YN20 = 0.7643;
    YNO0 = 0.0;
    YO0 = 0.0;
    YN0 = 0.0;

    //Turn into initial flow variables


    //Intialize flow


    nelem = 200;
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




    return 0;
}
