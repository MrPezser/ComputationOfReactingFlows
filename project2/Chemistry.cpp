//
// Created by Tsail on 3/16/2024.
//

#include "Chemistry.h"
#include "Indexing.h"
#include <iostream>


void LoadCurveFits(){
    //Function for reading in the therm.dat file and loading in thermo curve fits for 5 species air
    double Al[NSP], Bl[NSP], Cl[NSP], Dl[NSP], El[NSP], Fl[NSP], Gl[NSP];
    double Ah[NSP], Bh[NSP], Ch[NSP], Dh[NSP], Eh[NSP], Fh[NSP], Gh[NSP];
    double Mw[NSP];
    FILE* ftherm = fopen("../therm_airmix.dat","r");

    char test[15];

    for (int isp=0; isp<NSP; isp++) {
        //get molecular weight from first line
        fscanf(ftherm, "%*18s %*6c %*1c%*f%*1c %*f %*f %*2f %*s %*f %*f %lf %*d", &(Mw[isp]));
        //get curve fits
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Ah[isp]),&(Bh[isp]),&(Ch[isp]),&(Dh[isp]),&(Eh[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Fh[isp]),&(Gh[isp]),&(Al[isp]),&(Bl[isp]),&(Cl[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%*15lE%*d",&(Dl[isp]),&(El[isp]),&(Fl[isp]),&(Gl[isp]));
    }
}