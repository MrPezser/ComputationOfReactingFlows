//
// Created by Tsail on 3/16/2024.
//

#include "Chemistry.h"


void Chem::LoadCurveFits(){
    //Function for reading in the therm.dat file and loading in thermo curve fits for 5 species air
    FILE* ftherm = fopen("../therm_ONRP.dat","r");
    ASSERT(ftherm, "Unable to open thermo file!") //if null, exit

    for (int isp=0; isp<NUMSP; isp++) {
        //get molecular weight from first line
        fscanf(ftherm, "%*18s %*6c %*1c%*f%*1c %*f %*f %*2f %*s %*f %*f %lf %*d", &(Mw[isp]));
        //get curve fits
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Ah[isp]),&(Bh[isp]),&(Ch[isp]),&(Dh[isp]),&(Eh[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Fh[isp]),&(Gh[isp]),&(Al[isp]),&(Bl[isp]),&(Cl[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%*15lE%*d",&(Dl[isp]),&(El[isp]),&(Fl[isp]),&(Gl[isp]));

        Rs[isp] = Ruv/Mw[isp];
    }

    fclose(ftherm);
}

double Chem::Calc_h_Curve(int isp, double T) {
    // Calculate species specific enthalpy according to curve fit alone
    // isp = index of species
    // T = temperature
    // Cp = output Cp at temperature

    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = T*T2;
    T4 = T*T3;
    T5 = T*T4;

    ASSERT(T > 200.0, "Temp below curve fit")
    ASSERT(T < 10000.0, "Temp way above curve fit")

    if (T > 1000.0) {
        return (Rs[isp])*(Ah[isp]*T + Bh[isp]*T2*0.5 + Ch[isp]*T3/3.0 + Dh[isp]*T4*0.25 + Eh[isp]*T5*0.2 + Fh[isp]);
    } else {  //T<100
        return (Rs[isp])*(Al[isp]*T + Bl[isp]*T2*0.5 + Cl[isp]*T3/3.0 + Dl[isp]*T4*0.25 + El[isp]*T5*0.2 + Fl[isp]);
    }
}
double Chem::Calc_cp_Curve(int isp, double T) {
    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = T*T2;
    T4 = T*T3;

    ASSERT(T > 200.0, "Temp below curve fit")
    ASSERT(T < 10000.0, "Temp way above curve fit")

    if (T >= 1000.0) {
        return (Rs[isp])*(Ah[isp] + Bh[isp]*T + Ch[isp]*T2 + Dh[isp]*T3 + Eh[isp]*T4);
    } else {  //T<1000
        return (Rs[isp])*(Al[isp] + Bl[isp]*T + Cl[isp]*T2 + Dl[isp]*T3 + El[isp]*T4);
    }
}

