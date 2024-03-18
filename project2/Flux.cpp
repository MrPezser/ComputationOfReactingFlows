//
// Created by Tsail on 3/10/2024.
//

#include <valarray>
#include "Flux.h"
#include "Indexing.h"

#define sign(x)  ((std::signbit(x) ?  -1 : 1))

double wavespeed(const double* u, const double Ruv, const double* Mw, const double *Cv){
    double rhoR, rhoCv, rho, p, a2;

    rhoR = 0.0;
    rhoCv = 0.0;
    rho = 0.0;
    for (int isp=0; isp<NSP; isp++){
        rhoR  += u[isp]*Ruv/Mw[isp];
        rhoCv += u[isp]*Cv[isp];
        rho   += u[isp];
    }

    p = rhoR*u[NSP+1];
    a2 = (p/rho)*(1.0 + rhoR/rhoCv);

    return sqrt(a2);
}

void LDFSS(const double Ruv, const double A, const double* Mw, const double* Cv, const double* uL, const double* uR, double* flux) {
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
    aL = wavespeed(uL, Ruv, Mw, Cv);
    ahalf = 0.5*(aL+aR);

    //Calculate Pressure
    double pL, pR;
    pL = 0.0;
    pR = 0.0;
    for (int s=0; s<NSP; s++){
        pL += uL[s]*Ruv/Mw[s];
        pR += uR[s]*Ruv/Mw[s];
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
    //CalcEnthalpy(Ruv, Mw, uL, hoL);
    //CalcEnthalpy(Ruv, Mw, uR, hoR);

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
}


/*

c ---- Residual calculation

      do i=2,imx
       res(i,1:nsp+2) = (ev(i,1:nsp+2)-ev(i-1,1:nsp+2))/dx(i)
      enddo

c ---- add pressure source term

      do i=2,imx
       res(i,nsp+1) = res(i,nsp+1) - p(i)*(area(i)-area(i-1))/dx(i)
      enddo
 */