//
// Created by tskoepli on 3/27/2024.
//

#include "ChemSource.h"


double Calc_kf(double T, double Cf, double eta, double thetad){
    double coeff, expon;
    // Calculates the Ahrrenious forward reaction rate

    coeff = Cf * pow(T, eta);
    expon = -thetad / T;

    return coeff * exp( expon );
}

double Calc_keq(double T, double alpha, const double* A){
    double z, z2, z3, z4, coeff, expon;

    z = 10000.0 / T;
    z2 = z*z;
    z3 = z*z2;
    z4 = z*z3;

    coeff = alpha;
    expon = A[0] + A[1]*z + A[2]*z2+ A[3]*z3 + A[4]*z4;
    return coeff * exp( expon );
}


void CalcOmega(const double* unk, Chem& air, State& var, double* omega){
    double rrtb, kf, keq, kb, T, Cf, eta, thetad, alpha, concF, concB, nu, rho_tilde[NSP];
    double* A;
    int rxns, include;

    T = var.Tstar;
    for (int isp=0; isp<NSP; isp++){
        rho_tilde[isp] = unk[isp] / air.Mw[isp];
    }

    //Calculate third body contribution (save)
    for (int irx=0; irx<NRX; irx++){
        rrtb = 0.0;
        if (air.itb[irx] == 1){

            for (int jsp=0; jsp<NSP; jsp++){
                rrtb += unk[jsp] * air.tb[irx][jsp];
            }

        } else {
            rrtb = 1.0;
        }
        var.RRtb[irx] = rrtb;
    }

    //For each reaction calculate the reaction contribution (save)
    for (int irx=0; irx<NRX; irx++){
        //find kf
        Cf = air.Cf[irx];
        eta = air.eta[irx];
        thetad = air.thetad[irx];
        kf = Calc_kf(T,Cf,eta, thetad);

        //find keq/kb
        alpha = air.alpham[irx];
        A = air.Arxn[irx];
        keq = Calc_keq(T, alpha, A);
        kb = kf/keq;

        //calculate RR
        concF = 1.0;
        concB = 1.0;
        for(int jsp=0; jsp<NSP; jsp++){
            nu = air.nu[irx][jsp];
            for(int k=0; k<fabs(nu); k++){

                if (nu < 0){
                    concF *= rho_tilde[jsp];
                } else if (nu > 0){
                    concB *= rho_tilde[jsp];
                }
            }
        }

        var.RRlma[irx] = concF*kf - concB*kb;
    }

    //Now assemble the chemical source term
    for (int isp=0; isp<NSP; isp++){
        omega[isp]=0.0;

        rxns = air.sprxn[isp];
        for (int jrx = 0; jrx < NRX; jrx++) {
            include = (rxns >> jrx) & 1;  //read binary encoded val and det if rxn has contribution.
            //idk why I did this to save a few flops, was probably just avoiding coding the jacobian...

            if (include){
                if (air.itb[jrx]){
                    omega[isp] += var.RRlma[jrx]*rrtb;
                } else {
                    omega[isp] += var.RRlma[jrx];
                }
            }
        }
        omega[isp] *= air.Mw[isp];
    }

}