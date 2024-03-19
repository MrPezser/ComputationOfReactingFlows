//
// Created by Tsail on 3/10/2024.
//

#include <valarray>
#include "Flux.h"


#define sign(x)  ((std::signbit(x) ?  -1 : 1))

double wavespeed(const double* u, Chem& air){
    double rhoR, rhoCv, rho, p, a2, cp;

    rhoR = 0.0;
    rhoCv = 0.0;
    rho = 0.0;
    for (int isp=0; isp<NSP; isp++){
        cp = air.Calc_cp_curve(isp, u[NSP]);

        rhoR  += u[isp]*air.Ruv/air.Mw[isp];
        rhoCv += u[isp]*(cp - air.Ruv/air.Mw[isp]);
        rho   += u[isp];
    }

    p = rhoR*u[NSP+1];
    a2 = (p/rho)*(1.0 + rhoR/rhoCv);

    return sqrt(a2);
}

void LDFSS(const double A, const double* uL, const double* uR, Chem &air, double* flux) {
    /*
c --------------------------------------------------------------------
c ----- inviscid flux contribution (LDFSS)
c
c     rho - density
c     p - pressure
c     u - velocity
c     ho - stagnation enthalpy
c     ys - mass fractions
c     a - sound speed
c     area - interface area
c     dx   - mesh spacing for cell
c     res  - residual vector
c     ev - interface flux vector
c --------------------------------------------------------------------
    */
    ///MOVE CELL SPECIFIC CALCULATIONS (PRESSURE, DENSITY, ENTHALPY,....) OUTSIDE OF FLUX LOOP

    //Calculate wavespeed
    double aL, aR, ahalf;
    aL = wavespeed(uL, air);
    aR = wavespeed(uR, air);
    ahalf = 0.5*(aL+aR);

    //Calculate Pressure
    double pL, pR;
    pL = 0.0;
    pR = 0.0;
    for (int s=0; s<NSP; s++){
        pL += uL[s]*air.Ruv/air.Mw[s];
        pR += uR[s]*air.Ruv/air.Mw[s];
    }
    pL = pL * uL[NSP+1];
    pR = pR * uR[NSP+1];

    //Calculate Density
    double rhoL, rhoR;
    rhoL = 0.0; rhoR = 0.0;
    for (int s=0; s<NSP; s++) {
        rhoL += uL[s];
        rhoR += uR[s];
    }

    //Calculate Total Enthaply
    double hoL, hoR;
    hoL = air.Calc_h_Mix(uL) + pL + 0.5*uL[NSP]*uL[NSP]*rhoL;
    hoR = air.Calc_h_Mix(uR) + pR + 0.5*uR[NSP]*uR[NSP]*rhoR;

    // Flux Calculation
    double xml = uL[NSP]/ahalf;
    double xmr = uR[NSP]/ahalf;

    double all = 0.5*(1.0 + sign(xml));
    double alr = 0.5*(1.0 - sign(xmr));

    double btl = -fmax(0.0,1.0-double(int(fabs(xml))));
    double btr = -fmax(0.0,1.0-double(int(fabs(xmr))));

    double xmml =  0.25*(xml+1.0)*(xml+1.0);
    double xmmr = -0.25*(xmr-1.0)*(xmr-1.0);

    double xmhalf = sqrt(0.5*(xml*xml + xmr*xmr));
    double xmc = 0.25*btl*btr*(xmhalf - 1.0)*(xmhalf - 1.0);

    double delp = pL - pR;
    double psum = pL + pR;

    double xmcp = xmc*fmax(0.0,(1.0 - (delp/psum + 2.0*fabs(delp)/pL)));
    double xmcm = xmc*fmax(0.0,(1.0 + (delp/psum - 2.0*fabs(delp)/pR)));
    double cvlp = all*(1.+btl)*xml - btl*xmml;
    double cvlm = alr*(1.+btr)*xmr - btr*xmmr;
    double cep = cvlp - xmcp;
    double cem = cvlm + xmcm;

    double fml = A*rhoL*ahalf*cep;
    double fmr = A*rhoR*ahalf*cem;

    double ppl = 0.25*(xml+1.)*(xml+1.)*(2.-xml);
    double ppr = 0.25*(xmr-1.)*(xmr-1.)*(2.+xmr);

    double pnet = (all*(1.+btl) - btl*ppl)*pL
            + (alr*(1.+btr) - btr*ppr)*pR;

    for (int isp=0; isp<NSP; isp++) {
        flux[isp] = fml * uL[isp] +fmr * uR[isp];       //species  density
    }
    flux[NSP] = fml*uL[NSP]  + fmr*uR[NSP] + A*pnet;    //momentum
    flux[NSP+1] = fml*hoL + fmr*hoR;                    //total energy
    flux[NSP+2] = 0.0;                                  //vibrational energy
}

void CalcRes(int nelem, double dx,double CFL, Chem &air, double* u,double* Acc,double* Afa,double* dAdx, double* res) {
    //Find the common flux at each face
    double flux_comm[(nelem+1)*(NSP+3)];
    double *uL, *uR, *flux;

    //Interior faces
    for (int iface=1; iface<nelem; iface++) {
        uL = &(u[uIJK(iface - 1, 0, 0)]);
        uR = &(u[uIJK(iface, 0, 0)]);
        flux = &flux_comm[fIJ(iface, 0)];

        LDFSS(Afa[iface], uL, uR, air, flux);
    }
    //Boundary faces - just extrapolate both for now I guess - need to do this better
    uL = &(u[uIJK(0, 0, 0)]);
    uR = &(u[uIJK(0, 0, 0)]);
    flux = &flux_comm[fIJ(0, 0)];
    LDFSS(Afa[0], uL, uR, air, flux);

    uL = &(u[uIJK(nelem - 1, 0, 0)]);
    uR = uL;
    flux = &flux_comm[fIJ(nelem, 0)];
    LDFSS(Afa[nelem], uL, uR, air, flux);


    //Calculate RHS residual from flux scheme and pressure source term
    for (int ielem=0; ielem<nelem; ielem++) {
        for (int kvar=0; kvar<NSP+3; kvar++){
            res[fIJ(ielem, kvar)] = (flux[fIJ(ielem+1,kvar)] - flux[fIJ(ielem,kvar)]) / dx;
        }

        //Pressure source term
        double p;
        p = 0.0;
        for (int s=0; s<NSP; s++){
            p += u[uIJK(ielem,0,s)]*air.Ruv/air.Mw[s];
        }
        res[fIJ(ielem, NSP)] -= p*(Afa[ielem+1] - Afa[ielem])/dx;
    }
}