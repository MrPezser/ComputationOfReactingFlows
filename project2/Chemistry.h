//
// Created by Tsail on 3/16/2024.
//

#ifndef PROJECT2_CHEMISTRY_H
#define PROJECT2_CHEMISTRY_H

#define NRX 5

#include "Indexing.h"

#include <iostream>
#include <cmath>

class Chem {
private:

    //curve fit coefficients
    double Al[NSP]{}, Bl[NSP]{}, Cl[NSP]{}, Dl[NSP]{}, El[NSP]{}, Fl[NSP]{}, Gl[NSP]{};
    double Ah[NSP]{}, Bh[NSP]{}, Ch[NSP]{}, Dh[NSP]{}, Eh[NSP]{}, Fh[NSP]{}, Gh[NSP]{};
    void LoadCurveFits();

public:
    double e0[NSP]{}, Cf[NRX]{}, eta[NRX]{}, thetad[NRX]{}, tb[NRX][NSP]{}, Arxn[NRX][5]{}, alpham[NRX]{};
    double nu[NRX][NSP]{};
    int itb[NRX]{}, sprxn[NSP]{};

    //molecular weight
    double Mw[NSP]{};
    double Rs[NSP]{};
    double Ruv = 8314.34;
    double Tref = 298.15;
    double thetav[NSP]{};

    Chem() {
        LoadCurveFits();
        InitializeChemistry();
    }
    //thermo
    double Calc_h_Curve(int isp, double T);
    double Calc_rho_htr_Mix(const double* unk);
    double Calc_cp_Curve(int isp, double T);
    double Get_cptr(int isp, double T);
    double Get_htr(int isp, double T);
    double Get_hv(int isp, double T, double Tv);
    double Get_cpv(int isp, double T, double Tv);
    double Get_htotal(int isp, double T, double Tv);
    double Get_cptotal(int isp,double T,double Tv);

    //chem
    void InitializeChemistry();

    };

#endif //PROJECT2_CHEMISTRY_H
