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

    if (T >= 1000.0) {
        return (Rs[isp])*(Ah[isp] + Bh[isp]*T + Ch[isp]*T2 + Dh[isp]*T3 + Eh[isp]*T4);
    } else {  //T<1000
        return (Rs[isp])*(Al[isp] + Bl[isp]*T + Cl[isp]*T2 + Dl[isp]*T3 + El[isp]*T4);
    }
}

//Translational Thermo
double Chem::Get_cptr(int isp, double T){
    //if (isp>=3) return Calc_cp_Curve(isp,T);

    if (T > 1000.0) {
        return (Rs[isp]) * (Ah[isp]);
    }else{
        return (Rs[isp]) * (Al[isp]);
    }
}
double Chem::Get_htr(int isp, double T){
    //if (isp>=3) return Calc_h_Curve(isp,T);

    if (T > 1000.0) {
        return (Rs[isp]) * (Ah[isp]*(T) + Fh[isp]);
    } else {
        return (Rs[isp]) * (Al[isp]*(T) + Fl[isp]);
    }
}

//Vibrational Thermo
double Chem::Get_hv(int isp, double T, double Tv) {
    //if (isp>=3) return 0.0;

    double cptr, hv, h_f;
    cptr = Get_cptr(isp, T);
    hv = Calc_h_Curve(isp, Tv);

    if (T > 1000.0) {
        h_f = Fh[isp];
    } else {
        h_f = Fl[isp];
    }

    hv += -cptr*(Tv - Tref) - h_f;
    return hv;
}
double Chem::Get_cpv(int isp, double T, double Tv) {
    //if (isp>=3) return 0.0;

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

//Chemical Kinetics
void Chem::InitializeChemistry(){
    thetav[0] = 3395.0;
    thetav[1] = 2239.0;
    thetav[2] = 2817.0;
    thetav[3] = -1.0;
    thetav[4] = -1.0;

    //not given for N2, O2
    e0[2] = 2.996231e6;
    e0[3] = 3.362161e7;
    e0[4] = 1.543119e7;

    Cf[0] = 2.75e16;
    Cf[1] = 3.7e18;
    Cf[2] = 2.3e14;
    Cf[3] = 3.18e10;
    Cf[4] = 2.16e5;

    eta[0] = -1.0;
    eta[1] = -1.6;
    eta[2] = -0.5;
    eta[3] = 0.1;
    eta[4] = 1.29;

    thetad[0] = 59500.0;
    thetad[1] = 113200.0;
    thetad[2] = 75500.0;
    thetad[3] = 37700.0;
    thetad[4] = 19220.0;

    //3rd body efficiencies
    itb[0] = 1;
    itb[1] = 1;
    itb[2] = 1;
    itb[3] = 0;
    itb[4] = 0;

    //[reaction][species]
    tb[0][0] = 1.0;
    tb[0][1] = 1.0;
    tb[0][2] = 1.0;
    tb[0][3] = 3.0;
    tb[0][4] = 3.0;

    tb[1][0] = 1.0;
    tb[1][1] = 1.0;
    tb[1][2] = 1.0;
    tb[1][3] = 3.0;
    tb[1][4] = 3.0;

    tb[2][0] = 1.0;
    tb[2][1] = 1.0;
    tb[2][2] = 1.0;
    tb[2][3] = 2.0;
    tb[2][4] = 2.0;

    alpham[0] = 1000.0;
    alpham[1] = 1000.0;
    alpham[2] = 1000.0;
    alpham[3] = 1.0;
    alpham[4] = 1.0;

    //[reaction][coeff #]
    Arxn[0][0] = 1.335;
    Arxn[0][1] = -4.127;
    Arxn[0][2] = -0.616;
    Arxn[0][3] = 0.093;
    Arxn[0][4] = -0.005;

    Arxn[1][0] = 3.898;
    Arxn[1][1] = -12.611;
    Arxn[1][2] = 0.683;
    Arxn[1][3] = -0.118;
    Arxn[1][4] = 0.006;

    Arxn[2][0] = 1.549;
    Arxn[2][1] = -7.784;
    Arxn[2][2] = 0.228;
    Arxn[2][3] = -0.043;
    Arxn[2][4] = 0.002;

    Arxn[3][0] = 2.349;
    Arxn[3][1] = -4.828;
    Arxn[3][2] = 0.455;
    Arxn[3][3] = -0.075;
    Arxn[3][4] = 0.004;

    Arxn[4][0] = 0.215;
    Arxn[4][1] = -3.652;
    Arxn[4][2] = 0.843;
    Arxn[4][3] = -0.136;
    Arxn[4][4] = 0.007;

    //Differnece of species coefficients, v'' - v'  || [reaction][species]
    nu[0][1] = -1.0;
    nu[0][4] = 2.0;

    nu[1][0] = -1.0;
    nu[1][3] = 2.0;

    nu[2][2] = -1.0;
    nu[2][3] = 1.0;
    nu[2][4] = 1.0;

    nu[3][0] = -1.0;
    nu[3][2] = 1.0;
    nu[3][3] = 1.0;
    nu[3][4] = -1.0;

    nu[4][1] = 1.0;
    nu[4][2] = -1.0;
    nu[4][3] = 1.0;
    nu[4][4] = -1.0;

    //pointer from species to reaction (binary encoded, 1 2 4 8 16)
    sprxn[0] = 10;//N2 is in reaction 2 and 4, 2^1 + 2^3 = 2 + 8 = 10
    sprxn[1] = 17;//O2 is in reactions 1 and 5 = 1 + 16 = 17
    sprxn[2] = 28;//NO is in reactions 3, 4, and 5=  4 + 8 + 16 = 28
    sprxn[3] = 30;//N is in reactions 2, 3, 4, and 5= 30
    sprxn[4] = 29;//O is in reactions 1, 3, 4, 5 = 29
};
