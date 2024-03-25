//
// Created by tskoepli on 3/18/2024.
//

#include "VariableTransform.h"

void BuildDudv(double* unk, Chem &air, double** D) {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    //double D[NSP+3][NSP+3]{0.0};
    double rho_mix=0.0;
    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            D[i][j] = 0.0;
        }
    }
    //mixture density
    for (int isp=0; isp<NSP; isp++){
        rho_mix += unk[isp];
    }

    //d rho_s / d rho_s
    for (int i=0; i<NSP; i++){
        //Top left identity block
        D[i][i] = 1.0;
    }

    // momentum derivatives
    for (int j=0; j<NSP; j++){
        D[NSP][j] = unk[NSP];
    }
    D[NSP][NSP] = rho_mix;
    //D[NSP][NSP+1] = 0.0;
    //D[NSP][NSP+2] = 0.0;

    //total energy derivatives
    double ev[NSP], hv[NSP];
    double rhoCv, rhoCp, etr[NSP]{}, rhoR, cptr, htr[NSP], p{};
    rhoCv = 0.0;
    rhoCp = 0.0;
    rhoR = 0.0;
    for (int isp=0; isp<NSP; isp++){
        cptr = air.Get_cptr(isp, unk[NSP]);
        htr[isp] = air.Get_htr(isp, unk[NSP+1]);
        hv[isp] = air.Get_hv(isp, unk[NSP+1], unk[NSP+2]);

        rhoR  += unk[isp]*air.Ruv/air.Mw[isp];
        rhoCv += unk[isp]*(cptr - air.Ruv/air.Mw[isp]);
        rhoCp += unk[isp]*(cptr);
        etr[isp] = htr[isp] - (air.Ruv/air.Mw[isp])*unk[NSP+1];
        ev[isp] = hv[isp];// - (air.Ruv/air.Mw[isp])*unk[NSP+2];
        p += unk[isp]*air.Ruv/air.Mw[isp]*unk[NSP+1];
    }

    for (int isp=0; isp<NSP; isp++){
        D[NSP+1][isp] = etr[isp] + ev[isp] + 0.5*unk[NSP]*unk[NSP]; //- (rhoCv/rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1]
        D[NSP+1][NSP+2] += unk[isp]*ev[isp];
    }
    D[NSP+1][NSP]   = unk[NSP]*rho_mix;
    D[NSP+1][NSP+1] = rhoCv;


    //vibrational energy derivatives   -   Landau-Teller Source term contributions
    double tau[NSP];
    for (int isp=0; isp<NSP; isp++){
        double numer{}, denom{}, fracsum;
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

            frac = unk[jsp]/(air.Mw[jsp]*rho_mix);
            fracsum += frac;
            numer += frac;
            denom += frac / exp(Aij*(pow(unk[NSP+1],-1.0/3.0)-Bij) - 18.42);
        }
        tau[isp] = numer / ((p/101325.0)*denom);

        //Park vib relax time
        double sigma, c, n, NA;
        NA = 6.0221408e+23;
        sigma = 3.0e-21*(50000.0/unk[NSP+1]);
        c = sqrt( 8.0 * air.Ruv * unk[NSP+1] / (M_PI*air.Mw[isp]) );
        n = rho_mix*NA*fracsum;

        tau[isp] += 1.0 / (sigma * c * n);
    }

    for (int isp=0; isp<NSP; isp++){
        double evT, evTv, cpvT, cpvTv;
        evT = air.Get_hv(isp, unk[NSP+1], unk[NSP+1]);
        evTv= air.Get_hv(isp, unk[NSP+1], unk[NSP+2]);

        D[NSP+2][isp] = ev[isp];
        D[NSP+2][NSP+2] += unk[isp]*air.Get_cpv(isp,unk[NSP+1],unk[NSP+2]);

        cpvT = air.Get_cpv(isp, unk[NSP+1], unk[NSP+1]);
        cpvTv = air.Get_cpv(isp, unk[NSP+1], unk[NSP+2]);

        if (tau[isp] > 0) {
            //D[NSP+2][isp]   += -(evT - evTv) / tau[isp];
            //D[NSP+2][NSP+1] += -unk[isp]*cpvT / tau[isp];
            //D[NSP+2][NSP+2] += unk[isp]*cpvTv / tau[isp];
        }

    }

}