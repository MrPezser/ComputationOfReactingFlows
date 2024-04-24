//
// Created by Tsail on 3/16/2024.
//

#ifndef PROJECT2_CHEMISTRY_H
#define PROJECT2_CHEMISTRY_H

#include "Indexing.h"

#include <iostream>
#include <cmath>

class Chem {
private:

    //curve fit coefficients
    double Al[NUMSP]{}, Bl[NUMSP]{}, Cl[NUMSP]{}, Dl[NUMSP]{}, El[NUMSP]{}, Fl[NUMSP]{}, Gl[NUMSP]{};
    double Ah[NUMSP]{}, Bh[NUMSP]{}, Ch[NUMSP]{}, Dh[NUMSP]{}, Eh[NUMSP]{}, Fh[NUMSP]{}, Gh[NUMSP]{};
    void LoadCurveFits();

public:

    //molecular weight
    double Mw[NUMSP]{};
    double Rs[NUMSP]{};
    double Ruv= 8314.34;
    double RHOLREF = 800.0;
    double PREF = 101325.0;
    double AREF = 1300.0;
    double HVAP = 3.5e5;




    Chem() {
        LoadCurveFits();
    }
    //thermo
    double Calc_h_Curve(int isp, double T);
    double Calc_cp_Curve(int isp, double T);

    };

#endif //PROJECT2_CHEMISTRY_H
