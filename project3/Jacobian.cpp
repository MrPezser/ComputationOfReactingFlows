//
// Created by tskoepli on 3/18/2024.
//

#include "Jacobian.h"

void BuildJacobian(int isource, double dt, const double* unk, Chem &air, State& var, double** D) {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    //double D[NSP+3][NSP+3]{0.0};

    double dti = 1.0/dt;

    for (int i=0; i<NVAR; i++){
        for (int j=0; j<NVAR; j++){
            D[i][j] = 0.0;
        }
    }

    //Assuming that changes in density and energy from differing composition(within vapor) is negligible

    double Yv = unk[2];

    //      d__/d(yv1)
    D[0][0] =  dti * var.rho_mix*Yv;
    /// ADD d_rho/d_yv1 TERM

    // d__/d(yv2)
    D[1][1] =  dti * var.rho_mix*Yv;

    /// ADD d_rho/d_yv2 TERM

    //      d__/d(Yv)
    double drhodYv = -var.rho_mix*var.rho_mix * (1.0/var.rhov - 1.0/var.rhol);
    double dhdYv = var.hv - var.hl;

    D[0][2] = dti * (unk[0] * (drhodYv * Yv + var.rho_mix));
    D[1][2] = dti * (unk[1] * (drhodYv * Yv + var.rho_mix));
    D[2][2] = dti * (drhodYv * Yv + var.rho_mix);
    D[3][2] = dti * (drhodYv);
    D[4][2] = dti * (drhodYv * unk[4]);
    D[5][2] = dti * (drhodYv * var.h0 + var.rho_mix * dhdYv);
    D[6][2] = dti * (drhodYv * unk[6]);
    D[7][2] = dti * (drhodYv * unk[7]);

    //      d/d(p)
    D[0][3] = dti * (var.drhoP * Yv * unk[0]);
    D[1][3] = dti * (var.drhoP * Yv * unk[1]);
    D[2][3] = dti * (var.drhoP * Yv);
    D[3][3] = dti * (var.drhoP);
    D[4][3] = dti * (var.drhoP * unk[4]);
    D[5][3] = dti * (var.drhoP * var.h0 - 1.0);
    D[6][3] = dti * (var.drhoP * unk[6]);

    //      d/d(u)
    D[4][4] = dti * (var.rho_mix);
    D[5][4] = dti * (var.rho_mix*0.5*unk[4]);

    //      d/d(T)
    D[0][5] = dti * (var.drhoT * Yv * unk[0]);
    D[1][5] = dti * (var.drhoT * Yv * unk[1]);
    D[2][5] = dti * (var.drhoT * Yv);
    D[3][5] = dti * (var.drhoT);
    D[4][5] = dti * (var.drhoT * unk[4]);
    D[5][5] = dti * (var.drhoT * var.h0 + var.rho_mix*var.dhT - unk[3]/unk[5]); /// try adding this: dp/dT = P/T = u3/u5
    D[6][5] = dti * (var.drhoT * unk[6]);

    //      d/d(ntilde)
    D[6][6] = dti * var.rho_mix;

    // d / DelU
    D[7][7] = dti * var.rho_mix;

    for (int i=0; i<NVAR; i++){
        for (int j=0; j<NVAR; j++){
            ASSERT((!_isnan(D[i][j])), "NaN Jacobian")
        }
        /*
        if (D[i][i] <= 0.0){
            printf("negative jac diagonal\n");
        }*/
    }

}