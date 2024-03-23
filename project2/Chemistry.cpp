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

double Chem::Calc_h_Curve(int isp, double T) {
    // Calculate species specific enthalpy according to curve fit alone
    // isp = index of species
    // T = temperature
    // Cp = output Cp at temperature

    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = pow(T,3);
    T4 = pow(T,4);
    T5 = pow(T,5);


    if (T > 1000) {
        return (Ruv/Mw[isp])*(Ah[isp]*T + Bh[isp]*T2/2.0 + Ch[isp]*T3/3.0 + Dh[isp]*T4/4.0 + Eh[isp]*T5/5.0 + Fh[isp]);
    } else {  //T<100
        return (Ruv/Mw[isp])*(Al[isp]*T + Bl[isp]*T2/2.0 + Cl[isp]*T3/3.0 + Dl[isp]*T4/4.0 + El[isp]*T5/5.0 + Fl[isp]);
    }
}

double Chem::Calc_cp_Curve(int isp, double T) {
    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = pow(T,3);
    T4 = pow(T,4);

    if (T >= 1000) {
        return (Ruv/Mw[isp])*(Ah[isp] + Bh[isp]*T + Ch[isp]*T2 + Dh[isp]*T3 + Eh[isp]*T4);
    } else {  //T<1000
        return (Ruv/Mw[isp])*(Al[isp] + Bl[isp]*T + Cl[isp]*T2 + Dl[isp]*T3 + El[isp]*T4);
    }
}

//Translational Thermo
double Chem::Get_cptr(int isp, double T){
    if (T > 1000) {
        return (Ruv / Mw[isp]) * (Ah[isp]);
    }else{
        return (Ruv / Mw[isp]) * (Al[isp]);
    }
}
double Chem::Get_htr(int isp, double T){
    if (T > 1000) {
        return (Ruv / Mw[isp]) * (Ah[isp]*T + Fh[isp]);
    } else {
        return (Ruv / Mw[isp]) * (Al[isp]*T + Fl[isp]);
    }
}

//Vibrational Thermo
double Chem::Get_hv(int isp, double T, double Tv) {
    double cptr, hv;
    cptr = Get_cptr(isp, T);
    hv = Calc_h_Curve(isp, Tv);

    hv -= cptr*(Tv - Tref);
    return hv;
}
double Chem::Get_cpv(int isp, double T, double Tv) {
    double cptr, cpv;
    cptr = Get_cptr(isp, T);
    cpv = Calc_cp_Curve(isp, Tv);

    cpv -= cptr;
    return cpv;
}

//Total thermo
double Chem::Get_htotal(int isp, double T, double Tv){
    double hv, htr;
    htr = Get_htr(isp,T);
    hv = Get_hv(isp,T,Tv);

    return hv + htr;
}
double Chem::Get_cptotal(int isp,double T,double Tv){
    double cptr, cpv;
    cptr = Get_cptr(isp, T);
    cpv = Get_cpv(isp,T,Tv);

    return cptr + cpv;
}

double Chem::Calc_rho_htr_Mix(const double* unk) {
    //calculate the mixture enthalpy
    double hs[NSP], T, Tv, rho_h;
    T = unk[NSP+1];
    Tv = unk[NSP+2];
    rho_h = 0.0;

    for (int isp=0; isp<NSP; isp++) {
        hs[isp] = Get_htr(isp, T);
        rho_h += unk[isp]*hs[isp];
    }

    return rho_h;
}