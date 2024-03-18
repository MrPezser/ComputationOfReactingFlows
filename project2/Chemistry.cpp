//
// Created by Tsail on 3/16/2024.
//

#include "Chemistry.h"


void Chem::LoadCurveFits(){
    //Function for reading in the therm.dat file and loading in thermo curve fits for 5 species air
    FILE* ftherm = fopen("../therm_airmix.dat","r");

    for (int isp=0; isp<NSP; isp++) {
        //get molecular weight from first line
        fscanf(ftherm, "%*18s %*6c %*1c%*f%*1c %*f %*f %*2f %*s %*f %*f %lf %*d", &(Mw[isp]));
        //get curve fits
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Ah[isp]),&(Bh[isp]),&(Ch[isp]),&(Dh[isp]),&(Eh[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Fh[isp]),&(Gh[isp]),&(Al[isp]),&(Bl[isp]),&(Cl[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%*15lE%*d",&(Dl[isp]),&(El[isp]),&(Fl[isp]),&(Gl[isp]));
    }
}

void Chem::Calc_h_Curve(int isp, double T, double* hs) {
    // Calculate species specific enthalpy
    // isp = index of species
    // T = temperature
    // Cp = output Cp at temperature

    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = pow(T,3);
    T4 = pow(T,4);
    T5 = pow(T,5);

    if (T > 1000) {
        hs[isp] = (Ruv/Mw[isp])*(Ah[isp]*T + Bh[isp]*T2 + Ch[isp]*T3 + Dh[isp]*T4 + Eh[isp]*T5 + Fh[isp]);
    } else {  //T<1000
        hs[isp] = (Ruv/Mw[isp])*(Al[isp]*T + Bl[isp]*T2 + Cl[isp]*T3 + Dl[isp]*T4 + El[isp]*T5 + Fl[isp]);
    }
}

double Chem::Calc_dhdT(int isp, double T) {
    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = pow(T,3);
    T4 = pow(T,4);

    if (T >= 1000) {
        return (Ruv/Mw[isp])*(Ah[isp] + Bh[isp]*T*2 + Ch[isp]*T2*3 + Dh[isp]*T3*4 + Eh[isp]*T4*5);
    } else {  //T<1000
        return (Ruv/Mw[isp])*(Al[isp] + Bl[isp]*T*2 + Cl[isp]*T2*3 + Dl[isp]*T3*4 + El[isp]*T4*5);
    }
}

double Chem::Calc_h_Mix(const double* unk) {
    //calculate the mixture enthalpy
    double hs, T, rho_h;
    T = unk[NSP+1];
    rho_h = 0.0;

    for (int isp=0; isp<NSP; isp++) {
        Calc_h_Curve(isp, T, &hs);
        rho_h += unk[isp]*hs;
    }

    return rho_h;
}