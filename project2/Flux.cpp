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

double pressure(const double* u, Chem& air){
    //routine to calculate pressure given flow
    double p;
    p = 0.0;
    for (int s=0; s<NSP; s++){
        p += u[s]*air.Ruv/air.Mw[s];
    }
    p = p * u[NSP+1];
    return p;
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
    pL = pressure(uL,air);
    pR = pressure(uR,air);

    //Calculate Density
    double rhoL, rhoR;
    rhoL = 0.0; rhoR = 0.0;
    for (int s=0; s<NSP; s++) {
        rhoL += uL[s];
        rhoR += uR[s];
    }

    //Calculate Total Enthaply
    double hoL, hoR;
    hoL = air.Calc_rho_h_Mix(uL)/rhoL + 0.5*uL[NSP]*uL[NSP];//*rhoL;
    hoR = air.Calc_rho_h_Mix(uR)/rhoR + 0.5*uR[NSP]*uR[NSP];//*rhoR;

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

    double xmcp = xmc * fmax(0.0,(1.0 - (delp/psum + 2.0*fabs(delp)/pL)));
    double xmcm = xmc * fmax(0.0,(1.0 + (delp/psum - 2.0*fabs(delp)/pR)));
    double cvlp = all*(1.0+btl)*xml - btl*xmml;
    double cvlm = alr*(1.0+btr)*xmr - btr*xmmr;
    double cep = cvlp - xmcp;
    double cem = cvlm + xmcm;

    double fml = A*rhoL*ahalf*cep;
    double fmr = A*rhoR*ahalf*cem;

    double ppl = 0.25*(xml+1.0)*(xml+1.0)*(2.0-xml);
    double ppr = 0.25*(xmr-1.0)*(xmr-1.0)*(2.0+xmr);

    double pnet = (all*(1.0+btl) - btl*ppl)*pL
                + (alr*(1.0+btr) - btr*ppr)*pR;

    for (int isp=0; isp<NSP; isp++) {
        flux[isp] = fml*uL[isp]/rhoL + fmr*uR[isp]/rhoR;       //species  density
    }
    flux[NSP] = fml*uL[NSP]  + fmr*uR[NSP] + A*pnet;    //momentum
    flux[NSP+1] = fml*hoL + fmr*hoR;                    //total energy
    flux[NSP+2] = 0.0;                                  //vibrational energy
}


void EulerFlux(Chem& air, const double *u, double A, double* flux){
    //Convert to multispecies
    double rho_mix{0.0}, p{};

    for (int k=0; k<NSP+3; k++) {
        flux[k] = 0.0;
    }

    p = pressure(u,air);

    for (int isp=0; isp<NSP; isp++){
        flux[isp] = u[isp]*u[NSP];
        rho_mix += u[isp];
    }


    flux[NSP] = ( rho_mix*u[NSP]*u[NSP] + p);
    flux[NSP+1] = u[NSP]*( air.Calc_rho_h_Mix(u) + 0.5*rho_mix*u[NSP]*u[NSP] );
    flux[NSP+2] = 0.0;

    for (int k=0; k<NSP+3; k++) {
        flux[k] *= A;
    }

}

double F1pm(const int isPlus, const double M, const double rho, const double c){
    //Convert to multispecies
    if (isPlus==1){
        if (M<=-1.0) {
            //F1+
            return 0.0;
        }
        if (M>=1.0) {
            return rho * M * c;
        }
        return 0.25*rho*c*(M+1)*(M+1);
    }

    if (isPlus==0) {
        if (M <= -1.0) {
            //F1-
            return rho * M * c;// *-1
        }
        if (M >= 1.0) {
            return 0.0;
        }
        return -0.25 * rho * c * (M - 1) * (M - 1);
    }

    return NAN;
}

void LeerFlux(const double A, const double* uL, const double* uR, Chem& air, double *flux){
    //Convert to multispecies
    double F1L, F1R, fPlus[NSP+3], fMins[NSP+3];

    for (int k=0; k<NSP+3; k++) {
        flux[k] = 0.0;
    }

    //Calculate wavespeed
    double aL, aR, ahalf;
    aL = wavespeed(uL, air);
    aR = wavespeed(uR, air);
    ahalf = 0.5*(aL+aR);
    //Mach number
    double ML, MR;
    ML = uL[NSP] / ahalf;
    MR = uR[NSP] / ahalf;
    //Density
    double rhoL{0.0}, rhoR{0.0};
    for (int isp=0; isp<NSP; isp++) {
        rhoL += uL[isp];
        rhoR += uR[isp];
    }

    if (ML >= 1.0){
        EulerFlux(air, uL, A, flux);
        return;
    }

    if (MR <= -1.0){
        EulerFlux(air, uR, A, flux);
        return;
    }

    F1L = F1pm(1, ML, rhoL, aL);
    F1R = F1pm(0, MR, rhoR, aR);

    if(_isnan(F1L+F1R)){
        throw std::overflow_error("getting NAN mass flux!");
    }
    for (int isp=0; isp<NSP; isp++) {
        flux[isp] = A * (F1L*uL[isp]/rhoL + F1R*uR[isp]/rhoR);       //species  density
    }

    //should figure out a better way to do this
    double rhoCv, rhoRt, cp, h[NSP];
    rhoCv = 0.0;
    rhoRt = 0.0;
    for (int isp=0; isp<NSP; isp++) {
        cp = air.Calc_cp_curve(isp, uL[NSP]);
        air.Calc_h_Curve(isp, uL[NSP + 1], h);

        rhoRt += uL[isp] * air.Ruv / air.Mw[isp];
        rhoCv += uL[isp] * (cp - air.Ruv / air.Mw[isp]);
    }

    double gam = 1 + rhoRt/rhoCv;

    double A2 = ((gam - 1) * uL[NSP]) + (2.0 * aL); //vL
    fPlus[1] = F1L * A2 / gam;
    fPlus[2] = F1L * A2*A2 * 0.5 / (gam*gam - 1.0);

    A2 = -((gam - 1) * uR[NSP]) - (2.0 * aR); //vR
    fMins[1] = F1R * A2 / gam;
    fMins[2] = F1R * A2*A2 * 0.5 / (gam*gam - 1.0);

    flux[NSP]   = A * (fPlus[1] + fMins[1]);
    flux[NSP+1] = A * (fPlus[2] + fMins[2]);
    flux[NSP+2] = A*0.0;
}

void CalcRes(int nelem, double dx,double CFL, Chem &air, double* u0, double* u,double* Acc,double* Afa,const double* dAdx, double* res) {
    //Find the common flux at each face
    double flux_comm[(nelem+1)*(NSP+3)];
    double *uL, *uR, *flux;

    for (int iu=0; iu<nelem*(NSP+3)*NDEGR; iu++) {
        res[iu] = 0.0;
    }

    //Interior faces
    for (int iface=1; iface<nelem; iface++) {
        uL = &(u[uIJK(iface - 1, 0, 0)]);
        uR = &(u[uIJK(iface, 0, 0)]);
        flux = &flux_comm[fIJ(iface, 0)];

        LeerFlux(Afa[iface], uL, uR, air, flux);
    }
    //Boundary faces - just extrapolate RHS for now - need to do pressure BC
    uL = u0;
    uR = &(u[uIJK(0, 0, 0)]);
    flux = &flux_comm[fIJ(0, 0)];
    LeerFlux(Afa[0], uL, uR, air, flux);

    uL = &(u[uIJK(nelem - 1, 0, 0)]);
    uR = &(u[uIJK(nelem - 1, 0, 0)]);
    flux = &flux_comm[fIJ(nelem, 0)];
    LeerFlux(Afa[nelem], uL, uR, air, flux);


    //Calculate RHS residual from flux scheme and pressure source term
    for (int ielem=0; ielem<nelem; ielem++) {
        for (int kvar=0; kvar<(NSP+3); kvar++){
            res[fIJ(ielem, kvar)] = -(flux_comm[fIJ(ielem+1,kvar)] - flux_comm[fIJ(ielem,kvar)]) / dx;
        }

        //Pressure source term
        double p;
        p = 0.0;
        for (int s=0; s<NSP; s++){
            p += u[uIJK(ielem,0,s)] * (air.Ruv/air.Mw[s]) * u[uIJK(ielem,0,NSP+1)];
        }
        res[fIJ(ielem, NSP)] += p*(Afa[ielem+1] - Afa[ielem])/dx;//dAdx[ielem];
    }
}