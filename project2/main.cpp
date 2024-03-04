#include <iostream>
#include <cmath>

int main() {
    int nx;
    double dx, xmin, xmax, xc[nx], Ax[nx];

    nx = 200;

    xmin = 0.0;
    xmax = 1.5;
    dx = (xmax - xmin) / nx;
    for (int i=0; i<nx; i++) {
        xc[i] = i*dx + 0.5*dx; //cell centers

        double s = fmin(2.0, 1 + 4*(xc[i]-0.5)*(xc[i]-0.5));
        Ax[i] = M_PI*s*s;
    }

    return 0;
}
