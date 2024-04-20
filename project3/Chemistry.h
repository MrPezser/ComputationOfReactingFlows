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
    double Al[NSP]{}, Bl[NSP]{}, Cl[NSP]{}, Dl[NSP]{}, El[NSP]{}, Fl[NSP]{}, Gl[NSP]{};
    double Ah[NSP]{}, Bh[NSP]{}, Ch[NSP]{}, Dh[NSP]{}, Eh[NSP]{}, Fh[NSP]{}, Gh[NSP]{};
    void LoadCurveFits();

public:

    //molecular weight
    double Mw[NSP]{};
    double Rs[NSP]{};
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
