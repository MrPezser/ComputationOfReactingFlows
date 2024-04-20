//
// Created by tskoepli on 3/27/2024.
//

#ifndef PROJECT2_STATEVARIABLES_H
#define PROJECT2_STATEVARIABLES_H

#include "Indexing.h"
#include "Chemistry.h"
// class for storing calculated properties in an element

class State {

private:
    double* unk{};

public:
    double cptr[NSP]{},rho_mix{},rhoR{},rhoCv{},a{},h0{}, rhol{}, rhov{}, hv{}, hl{};

    State() = default;

    void Initialize(double* u){
        unk = u;
    }

    void UpdateState(Chem& air, int ireact) {
        rhoR = 0.0;
        rhoCv = 0.0;
        rho_mix = 0.0;
        h0=0.0;

        double y[3]{};
        y[0] = unk[0];
        y[1] = unk[1];
        y[2] = 1.0 - y[0] - y[1];
        double Yv = unk[2];
        double p = unk[3];
        double u = unk[4];
        double T = unk[5];

        rhol = air.RHOLREF + (p - air.PREF)/(air.AREF*air.AREF);
        rhov = p / (T*(y[0]*air.Rs[0] + y[1]*air.Rs[1] + y[2]*air.Rs[2]));
        rho_mix = 1.0 / ((Yv/rhov) + ((1.0-Yv)/rhol));

        rhoR = p/T;
        rhoCv = rhov * (y[0]*(air.Calc_cp_Curve(0,T)-air.Rs[0]) +
                        y[1]*(air.Calc_cp_Curve(1,T)-air.Rs[1]) +
                        y[2]*(air.Calc_cp_Curve(2,T)-air.Rs[2]));

        hv = y[0]*air.Calc_h_Curve(0,T) + y[1]*air.Calc_h_Curve(1,T) + y[2]*air.Calc_h_Curve(2,T);
        hl = air.Calc_h_Curve(2,T) - air.HVAP;

        h0 = Yv*hv + (1.0-Yv)*hl + 0.5*u*u;                       //Mixture Total Enthaply
        a = sqrt((p/rho_mix)*(1.0 + rhoR/rhoCv));        //Wavespeed

        ASSERT(!_isnan(a), "Failed to calculate sound speed!")
    }

};

#endif //PROJECT2_STATEVARIABLES_H
