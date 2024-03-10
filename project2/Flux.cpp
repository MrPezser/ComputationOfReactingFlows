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

void LDFSS(const double Ruv, const double* Mw, const double* Cv, const double* uL, const double* uR, double* flux) {
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
    double aL, aR, ahalf;
    // Flux Calculation

    //Calculate wavespeed
    aL = wavespeed(uL, Ruv, Mw, Cv);
    ahalf = 0.5*(aL+aR);
    double xml = uL[NSP]/ahalf;
    double xmr = uR[NSP]/ahalf;

    double all = 0.5*(1.0 + sign(xml));
    double alr = 0.5*(1.0 - sign(xmr));

    double btl = -fmax(0.0,1.0-double(int(fabs(xml))));
    double btr = -fmax(0.0,1.0-double(int(fabs(xmr))));

    double xmml =  0.25*(xml+1.0)*(xml+1.0);
    double xmmr = -0.25*(xmr-1.0)*(xmr-1.0);

}


/*



      xmhalf = sqrt(0.5*(xml*xml + xmr*xmr))
      xmc = 0.25*btl*btr*(xmhalf - 1.0)**2
      delp = p(i) - p(i+1)
      psum = p(i) + p(i+1)
      xmcp = xmc*max(0.0,(1.0 - (delp/psum + 2.0*abs(delp)/p(i))))
      xmcm = xmc*max(0.0,(1.0 + (delp/psum - 2.0*abs(delp)/p(i+1))))
      cvlp = all*(1.+btl)*xml - btl*xmml
      cvlm = alr*(1.+btr)*xmr - btr*xmmr
      cep = cvlp - xmcp
      cem = cvlm + xmcm

      fml = area(i)*rho(i)*ahalf*cep
      fmr = area(i)*rho(i+1)*ahalf*cem

      ppl = 0.25*(xml+1.)**2*(2.-xml)
      ppr = 0.25*(xmr-1.)**2*(2.+xmr)

      pnet(i) = (all*(1.+btl) - btl*ppl)*p(i)
     c        + (alr*(1.+btr) - btr*ppr)*p(i+1)

      ev(i,1:nsp) = fml*ys(i,1:nsp) + fmr*ys(i+1,1:nsp)           !species density
      ev(i,nsp+1) = fml*u(i)  + fmr*u(i+1) + area(i)*pnet(i)      !momentum
      ev(i,nsp+2) = fml*ho(i) + fmr*ho(i+1)                       !total energy

      enddo

c ---- Residual calculation

      do i=2,imx
       res(i,1:nsp+2) = (ev(i,1:nsp+2)-ev(i-1,1:nsp+2))/dx(i)
      enddo

c ---- add pressure source term

      do i=2,imx
       res(i,nsp+1) = res(i,nsp+1) - p(i)*(area(i)-area(i-1))/dx(i)
      enddo
 */