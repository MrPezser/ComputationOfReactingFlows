//
// Created by tskoepli on 3/18/2024.
//

#include "VariableTransform.h"
#include "Indexing.h"


void BuildDudv() {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    double D[NSP+3][NSP+3]{};

    for (int i=0; i<NSP; i++){
        //Top left identity portion
        D[i][i] = 1.0;

        // momentum derivatives
        for (int j=0; j<NSP; j++){
            D[i][j] = u;
        }
        D[i][NSP] = rho_mix;
        D[i][NSP+1] = 0.0;
        D[i][NSP+2] = 0.0;

        //total energy derivatives

        //vibrational energy derivatives
    }

}