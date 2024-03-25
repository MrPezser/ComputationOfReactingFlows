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
        cp = air.Get_cptr(isp, u[NSP+1]);  ///change this just to do total cp???

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
    hoL = 0.5*uL[NSP]*uL[NSP];// *rhoL;
    hoR =  0.5*uR[NSP]*uR[NSP];// *rhoR;

    for (int isp=0;isp<NSP;isp++){
        hoL += uL[isp] * air.Get_htotal(isp,uL[NSP+1],uL[NSP+2])/rhoL;
        hoR += uR[isp] * air.Get_htotal(isp,uR[NSP+1],uR[NSP+2])/rhoR;
    }

    //Calculate vibrational enthalpy
    double hvL{}, hvR{};

    for (int isp=0;isp<NSP;isp++){
        hvL += uL[isp] * air.Get_hv(isp,uL[NSP+1],uL[NSP+2])/rhoL;
        hvR += uR[isp] * air.Get_hv(isp,uR[NSP+1],uR[NSP+2])/rhoR;
    }

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
    flux[NSP+2] = fml*hvL + fmr*hvR;                                  //vibrational energy
}


void EulerFlux(Chem& air, const double *u, double A, double* flux){
    //Convert to multispecies
    double rho_mix{0.0}, p{}, ev[NSP]{}, hv[NSP]{}, rho_ev{0.0};

    for (int k=0; k<NSP+3; k++) {
        flux[k] = 0.0;
    }

    p = pressure(u,air);

    for (int isp=0; isp<NSP; isp++){
        flux[isp] = u[isp]*u[NSP];
        rho_mix += u[isp];
        hv[isp] = air.Get_hv(isp, u[NSP+1], u[NSP+2]);
        rho_ev += u[isp]*hv[isp];
    }


    flux[NSP] = ( rho_mix*u[NSP]*u[NSP] + p);
    flux[NSP+1] = u[NSP]*( air.Calc_rho_htr_Mix(u) + rho_ev + 0.5*rho_mix*u[NSP]*u[NSP] );
    flux[NSP+2] = u[NSP]*rho_ev;

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
    double F1L, F1R, fPlus[3], fMins[3];

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
    double FmyspL[NSP],FmyspR[NSP];
    for (int isp=0; isp<NSP; isp++) {
        FmyspL[isp] = uL[isp]/rhoL;
        FmyspR[isp] = uR[isp]/rhoR;
        flux[isp] = A * (F1L*FmyspL[isp] + F1R*uR[isp]/rhoR);       //species  density
    }

    double rhoCvL{}, rhoRtL{}, rhoCvR{}, rhoRtR{}, EvL{}, EvR{}, cp, hL{},hR{}, pL{}, pR{};
    for (int isp=0; isp<NSP; isp++) {
        cp = air.Get_cptr(isp, uL[NSP+1]);
        rhoRtL  += uL[isp] * air.Ruv / air.Mw[isp] * uL[NSP+1];
        rhoCvL  += uL[isp] * (cp - air.Ruv / air.Mw[isp]);

        EvL += FmyspL[isp] * air.Get_hv(isp,uL[NSP+1], uL[NSP+2]);
        hL += uL[isp] * air.Get_htotal(isp,uL[NSP+1],uL[NSP+2]);
        pL += uL[isp] * (air.Ruv/air.Mw[isp]) * uL[NSP+1];

        cp = air.Get_cptr(isp, uR[NSP+1]);
        rhoRtR  += uR[isp] * air.Ruv / air.Mw[isp] * uR[NSP+1];
        rhoCvR  += uR[isp] * (cp - air.Ruv / air.Mw[isp]);

        EvR += FmyspR[isp] * air.Get_hv(isp,uR[NSP+1], uR[NSP+2]);
        hR += uR[isp] * air.Get_htotal(isp,uR[NSP+1],uR[NSP+2]);
        pR += uR[isp] * (air.Ruv/air.Mw[isp]) * uR[NSP+1];
    }

    //double gamL = 1 + rhoRtL/rhoCvL;
    //double gamR = 1 + rhoRtR/rhoCvR;
    double ppl = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
    double ppr = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);

    double pnet = ppl*pL + ppr*pR;

    //double A2 = ((gamL - 1) * uL[NSP]) + (2.0 * aL); //vL
    fPlus[0] = F1L * (uL[NSP]);     //A2 / gamL;
    fPlus[1] = F1L * (hL/rhoL + 0.5*uL[NSP]*uL[NSP]);            //A2*A2 * 0.5 / (gamL*gamL - 1.0);
    fPlus[2] = F1L * EvL;

    //A2 = -((gamR - 1) * uR[NSP]) - (2.0 * aR); //vR
    fMins[0] = F1R * (uR[NSP]);           //A2 / gamR;
    fMins[1] = F1R * (hR/rhoR + 0.5*uR[NSP]*uR[NSP]);           // F1R * A2*A2 * 0.5 / (gamR*gamR - 1.0);
    fPlus[2] = F1R * EvR;

    flux[NSP]   = A * (fPlus[0] + fMins[0] + pnet);
    flux[NSP+1] = A * (fPlus[1] + fMins[1]);
    flux[NSP+2] = A * (fPlus[2] + fMins[2]);
}

void PressureBC(double pb, Chem& air,const double* u, double* uGhost){
    double T{}, rhoR{}, rho_ratio, denom;
    T = u[NSP+1];
    for (int isp=0;isp<NSP;isp++){
        rhoR += u[isp]*air.Ruv/air.Mw[isp];
    }
    rho_ratio = pb/(rhoR*T);

    for (int isp=0;isp<NSP;isp++){
        uGhost[isp] = (u[isp])*rho_ratio;
    }
    uGhost[NSP] = u[NSP]/rho_ratio;
    uGhost[NSP+1] = u[NSP+1];
    uGhost[NSP+2] = u[NSP+2];

}



void TNERelaxation(int nelem, Chem& air, double* u, double* tne_flux){
    //Calculate the Landau-Teller vib/tr energy relaxation process with Millikan and White timescale
    double tau[NSP]{}, evT, evTv, T, Tv, rho_mix, p, fracsum{};

    for (int ielem=0; ielem<nelem; ielem++) {
        int id = uIJK(ielem, 0, 0);
        tne_flux[id+NSP+2] = 0.0;
        double* ui = &(u[id]);
        T = ui[NSP+1];
        Tv = ui[NSP+2];
        if (T != Tv){
            //printf("Thermal neq\n");
        }

        rho_mix = 0.0;
        p=0.0;
        for (int isp=0; isp<NSP; isp++){
            rho_mix += ui[isp];
            p += ui[isp]* (air.Ruv/air.Mw[isp]) * T;
        }

        for (int isp=0; isp<NSP; isp++){
            double numer{}, denom{};
            fracsum = 0.0;

            if (air.thetav[isp] <= 0.0) {
                tau[isp] = -1.0;
                continue;
            }
            for (int jsp=0; jsp<NSP; jsp++) {
                double frac, Aij, Bij, muij;

                muij = air.Mw[isp]*air.Mw[jsp] / (air.Mw[isp] + air.Mw[jsp]);

                Aij = (1.16e-3) * sqrt(muij) * pow( air.thetav[isp], 4.0/3.0);
                Bij = 0.015 * pow(muij, 0.25);

                frac = ui[jsp]/(air.Mw[jsp]*rho_mix);
                fracsum += frac;
                numer += frac;
                denom += exp(Aij*(pow(T,-1.0/3.0)-Bij) - 18.42)/(p/101325.0);
            }
            tau[isp] = numer / (denom);

            //Park vib relax time
            double sigma, c, n, NA;
            NA = 6.0221408e+23;
            sigma = 3.0e-21*(50000.0/T);
            c = sqrt( 8.0 * air.Ruv * T / (M_PI*air.Mw[isp]) );
            n = rho_mix*NA*fracsum;

            tau[isp] += 1.0 / (sigma * c * n);
        }

        for (int isp = 0; isp < NSP; isp++) {
            evT = air.Get_hv(isp, T, T);
            evTv= air.Get_hv(isp, T, Tv);

            if (tau[isp] <= 0.0) continue;
            tne_flux[id+NSP+2] += ui[isp] * (evT - evTv) / tau[isp];
        }
    }

}

void CalcRes(int nelem, double dx, double CFL, double pb, Chem &air, double* u0, double* u,double* Acc,double* Afa,const double* dAdx, double* res) {
    //Find the common flux at each face
    double flux_comm[(nelem+1)*(NSP+3)], uBack[NSP+3];
    double tne_flux[(nelem+1)*(NSP+3)];
    double *uL, *uR, *flux;

    for (int iu=0; iu<nelem*(NSP+3)*NDEGR; iu++) {
        res[iu] = 0.0;
        flux_comm[iu] = 0.0;
    }



    //Interior faces
    for (int iface=1; iface<nelem; iface++) {
        uL = &(u[uIJK(iface - 1, 0, 0)]);
        uR = &(u[uIJK(iface, 0, 0)]);
        flux = &flux_comm[fIJ(iface, 0)];

        LDFSS(Afa[iface], uL, uR, air, flux);
    }
    //Boundary faces - pressure BC outflow
    //Inflow
    uL = u0;
    uR = &(u[uIJK(0, 0, 0)]);
    flux = &flux_comm[fIJ(0, 0)];
    LDFSS(Afa[0], uL, uR, air, flux);
    //Outflow
    uL = &(u[uIJK(nelem - 1, 0, 0)]);
    PressureBC(pb, air, uL, uBack);
    uR = uBack;
    flux = &flux_comm[fIJ(nelem, 0)];
    LDFSS(Afa[nelem], uL, uR, air, flux);


    //Calculate RHS residual from flux scheme and pressure source term
    TNERelaxation(nelem, air, u, tne_flux);

    for (int ielem=0; ielem<nelem; ielem++) {
        for (int kvar=0; kvar<(NSP+3); kvar++){
            res[fIJ(ielem, kvar)] = (- (flux_comm[fIJ(ielem+1,kvar)] - flux_comm[fIJ(ielem,kvar)]) / dx);
        }

        res[fIJ(ielem, NSP+2)] += Acc[ielem]*tne_flux[fIJ(ielem,NSP+2)];

        //Pressure source term
        double p;
        p = 0.0;
        for (int s=0; s<NSP; s++){
            p += u[uIJK(ielem,0,s)] * (air.Ruv/air.Mw[s]) * u[uIJK(ielem,0,NSP+1)];
        }
        res[fIJ(ielem, NSP)] += p*(Afa[ielem+1] - Afa[ielem])/dx; // dAdx[ielem];
    }
}