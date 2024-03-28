//
// Created by tskoepli on 3/18/2024.
//

#include "Jacobian.h"

void BuildJacobian(double* unk, Chem &air, State& var, double** D) {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    //double D[NSP+3][NSP+3]{0.0};

    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            D[i][j] = 0.0;
        }
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
    D[NSP][NSP] = var.rho_mix;

    //total energy derivatives

    for (int isp=0; isp<NSP; isp++){
        D[NSP+1][isp] = var.etr[isp] + var.evTv[isp] + 0.5*unk[NSP]*unk[NSP];//- (rhoCv/rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1];
        D[NSP+1][NSP+2] += unk[isp]*air.Get_cpv(isp,unk[NSP+1],unk[NSP+2]);
    }
    D[NSP+1][NSP]   = unk[NSP]*var.rho_mix;
    D[NSP+1][NSP+1] = var.rhoCv;


    //vibrational energy derivatives   -   Landau-Teller Source term contributions

    for (int isp=0; isp<NSP; isp++){
        D[NSP+2][isp] = var.evTv[isp];
        D[NSP+2][NSP+2] += unk[isp]*air.Get_cpv(isp,unk[NSP+1],unk[NSP+2]);

        if (var.taui[isp] > 0) {
            D[NSP+2][isp]   += -(var.evT[isp] - var.evTv[isp]) * var.taui[isp];
            D[NSP+2][NSP+1] += -unk[isp]*var.cpvT[isp]         * var.taui[isp];
            D[NSP+2][NSP+2] +=  unk[isp]*var.cpvTv[isp]        * var.taui[isp];
        }

    }

}