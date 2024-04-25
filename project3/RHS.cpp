//
// Created by Tsail on 3/10/2024.
//

#include <valarray>
#include "RHS.h"

#define K  (0.0)
#define sign(x)  ((std::signbit(x) ?  -1 : 1))

void LDFSS(const double A, const double* uL, State& varL, const double* uR, State& varR, Chem &air, double* flux) {
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

    double ahalf = 0.5 * (varL.a + varR.a);

    // Flux Calculation
    double xml = uL[4]/ahalf;
    double xmr = uR[4]/ahalf;

    double all = 0.5*(1.0 + sign(xml));
    double alr = 0.5*(1.0 - sign(xmr));

    double btl = -fmax(0.0,1.0-double(int(fabs(xml))));
    double btr = -fmax(0.0,1.0-double(int(fabs(xmr))));

    double xmml =  0.25*(xml+1.0)*(xml+1.0);
    double xmmr = -0.25*(xmr-1.0)*(xmr-1.0);

    double xmhalf = sqrt(0.5*(xml*xml + xmr*xmr));
    double xmc = 0.25*btl*btr*(xmhalf - 1.0)*(xmhalf - 1.0);

    double delp = uL[3] - uR[3];
    double psum = uL[3] + uR[3];

    double xmcp = xmc * fmax(0.0, 1.0 - (0.5*(delp + fabs(delp))/(varL.rho_mix*ahalf*ahalf)) ); //delp/psum + 2.0*fabs(delp)/uL[3]))
    double xmcm = xmc * fmax(0.0, 1.0 + (0.5*(delp - fabs(delp))/(varR.rho_mix*ahalf*ahalf)) ); //delp/psum - 2.0*fabs(delp)/uR[3]))


    double cvlp = all*(1.0+btl)*xml - btl*xmml;
    double cvlm = alr*(1.0+btr)*xmr - btr*xmmr;
    double cep = cvlp - xmcp;
    double cem = cvlm + xmcm;

    double fml = A*varL.rho_mix*ahalf*cep;
    double fmr = A*varR.rho_mix*ahalf*cem;

    double ppl = 0.25*(xml+1.0)*(xml+1.0)*(2.0-xml);
    double ppr = 0.25*(xmr-1.0)*(xmr-1.0)*(2.0+xmr);

    double pnet = (all*(1.0+btl) - btl*ppl)*uL[3]
                + (alr*(1.0+btr) - btr*ppr)*uR[3];


    double YvL = uL[2];
    double YvR = uR[2];
    const double* yL = uL;
    const double* yR = uR;
    for (int isp=0; isp<2; isp++) {
        flux[isp] = fml*YvL*yL[isp] + fmr*YvR*yR[isp];       //vapor mass fraction
    }
    flux[2] = fml*YvL     + fmr*YvR;
    flux[3] = fml         + fmr;
    flux[4] = fml*uL[4]   + fmr*uR[4] + A*pnet;    //momentum
    flux[5] = fml*varL.h0 + fmr*varR.h0;              //total energy
    flux[6] = fml*uL[6]   + fmr*uR[6];
}

double VaporSource(const double A, const double* unk, State& var, Chem& air){
    double Bvap, Sh, Re, cpmix, muv, dp, DeltaU, Pr, Sc, Sv;
    double T = unk[5];

    Pr = 0.7;
    Sc = 0.7;

    DeltaU = K * unk[4];

    dp = cbrt( 6.0*(1.0 - unk[2])/(M_PI*var.rhol*unk[6]) );
    if (dp < 1e-6) return 0.0;

    muv = (1.716e-5) * pow(T/273.15, 3.0/2.0) * ((273.15 + 110.4)/(T + 110.4));
    cpmix = unk[0]*air.Calc_cp_Curve(0,T) + unk[1]*air.Calc_cp_Curve(1,T) + (1.0-unk[0]-unk[1])*air.Calc_cp_Curve(2,T);
    Re = (var.rhov * DeltaU * dp) / muv;
    Sh = 2.0 + 0.552*sqrt(Re)* cbrt(Sc);
    Bvap = cpmix * fmax(0, T-444.0) / air.HVAP;
    Sv = 2.0 * M_PI * var.rho_mix * unk[6] * dp * Sh * (muv/Pr) * log(1 + Bvap);
    ASSERT(!_isnan(Sv),"NaN source term")

    return A * Sv;

}

void PressureBC(double pb, const double* u, double* uGhost){

    for (int kv=0; kv<NVAR; kv++){
        uGhost[kv] = u[kv];
    }
    uGhost[3] = pb;

}

void CalcRes(int isource, int nelem, double dx, double CFL, double pb, Chem &air, State* ElemVar, double* u0, double* u,
             const double* Acc, const double* Afa, const double* dAdx, double* res) {
    //Find the common flux at each face
    double flux_comm[(nelem+1)*(NVAR)], uBack[NVAR];
    double src{};
    double *uL, *uR, *flux;
    State varL, varR;

    for (int iu=0; iu<nelem*(NVAR)*NDEGR; iu++) {
        res[iu] = 0.0;
        flux_comm[iu] = 0.0;
    }

    //Interior faces
    for (int iface=1; iface<nelem; iface++) {
        uL = &(u[uIJK(iface - 1, 0, 0)]);
        varL = ElemVar[iface-1];
        uR = &(u[uIJK(iface, 0, 0)]);
        varR = ElemVar[iface];

        flux = &flux_comm[fIJ(iface, 0)];

        LDFSS(Afa[iface], uL, varL, uR, varR, air, flux);
    }
    //Boundary faces - pressure BC outflow

    //Inflow
    uL = u0;
    varL.Initialize(u0); //freestream vars
    varL.UpdateState(air, isource);
    uR = &(u[uIJK(0, 0, 0)]);
    varR = ElemVar[0];

    flux = &flux_comm[fIJ(0, 0)];
    LDFSS(Afa[0], uL, varL, uR, varR, air, flux);

    //Outflow
    uL = &(u[uIJK(nelem - 1, 0, 0)]);
    varL = ElemVar[nelem-1];

    PressureBC(pb, uL, uBack);
    uR = uBack;
    varR.Initialize(uR);
    varR.UpdateState(air, isource);

    flux = &flux_comm[fIJ(nelem, 0)];
    LDFSS(Afa[nelem], uL, varL, uR, varR, air, flux);


    //Calculate RHS residual from flux differencing and source terms
    //Thermal nonequilibrium source terms
    //TNERelaxation(nelem, air, ElemVar, u, tne_src);

    for (int ielem=0; ielem<nelem; ielem++) {
        int id = uIJK(ielem,0,0);
        double* ui = &(u[id]);
        //State var = ElemVar[ielem];

        //Flux
        for (int kvar=0; kvar<(NVAR); kvar++){
            int idk = id+kvar;
            res[idk] = -(flux_comm[uIJK(ielem+1,0,kvar)] - flux_comm[idk]) / dx;
        }
        //Pressure source term
        res[id+4] += ui[3]*(Afa[ielem+1] - Afa[ielem])/dx; // dAdx[ielem];

        //vapor source term
        if (isource == 1) {
            res[id + 2] += VaporSource(Acc[ielem], ui, ElemVar[ielem], air);
        }

    }
}