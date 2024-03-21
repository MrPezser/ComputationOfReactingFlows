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
    D[NSP][NSP+1] = 0.0;
    D[NSP][NSP+2] = 0.0;

    //total energy derivatives
    double rhoCv, rhoCp, e[NSP]{}, rhoR, cp, h[NSP];
    rhoCv = 0.0;
    rhoCp = 0.0;
    rhoR = 0.0;
    for (int isp=0; isp<NSP; isp++){
        cp = air.Calc_cp_curve(isp, unk[NSP]);
        air.Calc_h_Curve(isp, unk[NSP+1], h);

        rhoR  += unk[isp]*air.Ruv/air.Mw[isp];
        rhoCv += unk[isp]*(cp - air.Ruv/air.Mw[isp]);
        rhoCp += unk[isp]*(cp);
        e[isp] = h[isp] - (air.Ruv/air.Mw[isp])*unk[NSP+1];
    }

    for (int isp=0; isp<NSP; isp++){
        D[NSP+1][isp] = e[isp] + 0.5*unk[NSP]*unk[NSP]; //- (rhoCv/rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1]
    }
    D[NSP+1][NSP]   = unk[NSP]*rho_mix;
    D[NSP+1][NSP+1] = rhoCv - rho_mix*0.5*unk[NSP]*unk[NSP] / unk[NSP+1];
    D[NSP+1][NSP+2] = 0.0; //leave 0 for now


    //vibrational energy derivatives
    //leave out for now
    D[NSP+2][NSP+2] = 1.0;

}