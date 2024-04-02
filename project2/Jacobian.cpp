//
// Created by tskoepli on 3/18/2024.
//

#include "Jacobian.h"

void BuildJacobian(int ireact, double dt, const double* unk, Chem &air, State& var, double** D) {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    //double D[NSP+3][NSP+3]{0.0};

    double dti = 1.0/dt;

    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            D[i][j] = 0.0;
        }
    }

    //d rho_s / d rho_s
    for (int i=0; i<NSP; i++){
        //Top left identity block
        D[i][i] = dti * 1.0;
    }

    // momentum derivatives
    for (int j=0; j<NSP; j++){
        D[NSP][j] = unk[NSP];
    }
    D[NSP][NSP] = dti * var.rho_mix;

    //total energy derivatives

    for (int isp=0; isp<NSP; isp++){
        D[NSP+1][isp] = dti * (var.etr[isp] - 0.5*unk[NSP]*unk[NSP]);//- (var.rhoCv/var.rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1]);
    }
    D[NSP+1][NSP]   = dti * unk[NSP]*var.rho_mix;
    D[NSP+1][NSP+1] = dti * var.rhoCv;



    if (ireact==1) {
    for (int isp=0; isp<NSP; isp++){
        //vibrational energy derivatives
        D[NSP+2][isp] = dti * var.evTv[isp];
        D[NSP+2][NSP+2] += dti * unk[isp]*var.cpvTv[isp];

        //other vib contributions
        D[NSP+1][isp] += var.evTv[isp];
        D[NSP+1][NSP+2] += dti * unk[isp]*var.cpvTv[isp];


        // Jacobian for Landau-Teller Relaxation Source Term
        if (var.taui[isp] > 0) { //flag for
            //D[NSP+2][isp]   += -(var.evT[isp] - var.evTv[isp]) * var.taui[isp];
            D[NSP+2][NSP+1] += -unk[isp]*var.cpvT[isp]  * var.taui[isp];
            D[NSP+2][NSP+2] +=  unk[isp]*var.cpvTv[isp] * var.taui[isp];
        }
    }


        //Chemical Source Term Jacobians
        //need array of r derivatives wrt each species density
        // 0,N2     1,O2    2,NO    3,N     4,O
        double dRdrho[NRX][NSP]{0.0}, rho_tilde[NSP];
        int irx{};
        for (int isp = 0; isp < NSP; isp++) {
            rho_tilde[isp] = unk[isp] / air.Mw[isp];
        }

        irx = 0;
        dRdrho[irx][0] = var.RRlma[irx] * air.tb[irx][0] / air.Mw[0];
        dRdrho[irx][1] = (var.kf[irx] / air.Mw[1]) * var.RRtb[irx] + var.RRlma[irx] * (air.tb[irx][1] / air.Mw[1]);
        dRdrho[irx][2] = var.RRlma[irx] * air.tb[irx][2] / air.Mw[2];
        dRdrho[irx][3] = var.RRlma[irx] * air.tb[irx][3] / air.Mw[3];
        dRdrho[irx][4] = -2.0 * var.kb[irx] * rho_tilde[4] / air.Mw[4] + var.RRlma[irx] * air.tb[irx][4] / air.Mw[4];

        irx = 1;
        dRdrho[irx][0] = (var.kf[irx] / air.Mw[0]) * var.RRtb[irx] + var.RRlma[irx] * (air.tb[irx][0] / air.Mw[0]);
        dRdrho[irx][1] = var.RRlma[irx] * air.tb[irx][1] / air.Mw[1];
        dRdrho[irx][2] = var.RRlma[irx] * air.tb[irx][2] / air.Mw[2];
        dRdrho[irx][3] = -2.0 * var.kb[irx] * rho_tilde[3] / air.Mw[3] + var.RRlma[irx] * air.tb[irx][3] / air.Mw[3];
        dRdrho[irx][4] = var.RRlma[irx] * air.tb[irx][4] / air.Mw[4];

        irx = 2;
        dRdrho[irx][0] = var.RRlma[irx] * air.tb[irx][0] / air.Mw[0];
        dRdrho[irx][1] = var.RRlma[irx] * air.tb[irx][1] / air.Mw[1];
        dRdrho[irx][2] = (var.kf[irx] / air.Mw[2]) * var.RRtb[irx] + var.RRlma[irx] * air.tb[irx][2] / air.Mw[2];
        dRdrho[irx][3] =
                -(var.kb[irx] * rho_tilde[4] / air.Mw[3]) * var.RRtb[irx] + var.RRlma[irx] * air.tb[irx][3] / air.Mw[3];
        dRdrho[irx][4] =
                -(var.kb[irx] * rho_tilde[3] / air.Mw[4]) * var.RRtb[irx] + var.RRlma[irx] * air.tb[irx][4] / air.Mw[4];

        irx = 3;
        dRdrho[irx][0] = var.kf[irx] * rho_tilde[4] / air.Mw[0];
        dRdrho[irx][1] = 0.0;
        dRdrho[irx][2] = -var.kb[irx] * rho_tilde[3] / air.Mw[2];
        dRdrho[irx][3] = -var.kb[irx] * rho_tilde[2] / air.Mw[3];
        dRdrho[irx][4] = var.kf[irx] * rho_tilde[0] / air.Mw[4];

        irx = 4;
        dRdrho[irx][0] = 0.0;
        dRdrho[irx][1] = -var.kb[irx] * rho_tilde[3] / air.Mw[1];
        dRdrho[irx][2] = var.kf[irx] * rho_tilde[4] / air.Mw[2];
        dRdrho[irx][3] = -var.kb[irx] * rho_tilde[1] / air.Mw[3];
        dRdrho[irx][4] = var.kf[irx] * rho_tilde[2] / air.Mw[4];
        // 0,N2     1,O2    2,NO    3,N     4,O

        double SChemJac;

        for (int isp = 0; isp < NSP; isp++) {
            for (int jsp = 0; jsp < NSP; jsp++) {
                SChemJac = 0.0;
                // d_omega_s(i) / d_rho_s(j)
                //d[i][j] = Mw[i] * (sum of dRj with coefficients from omega_i)
                // air.nu[rx][sp] = (v'' - v') for species sp in reaxtion rx (related to omega[i])

                for (int krx = 0; krx < NSP; krx++) {
                    SChemJac += air.nu[krx][isp] * dRdrho[krx][jsp];//air.Mw[isp];
                }

                D[isp][jsp] += air.Mw[isp] * SChemJac;
            }
        }

        ///!!!!!!!Still Need d omega / d temperature Jacobian terms
    }

    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            if (_isnan(D[i][j])){
                printf("stop here");
            }
            ASSERT(!_isnan(D[i][j]), "Not a Number jacobian value")
        }
    }

}