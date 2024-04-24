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
    double rho_mix{},a{},h0{}, rhol{}, rhov{}, hv{}, hl{}, dhT{}, dhp{}, drhoT{}, drhoP{};
    //double rhoR{},rhoCv{};

    State() = default;

    void Initialize(double* u){
        unk = u;
    }

    void UpdateState(Chem& air, int ireact) {
        //rhoR = 0.0;
        //rhoCv = 0.0;
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

        //rhoR = p/T;
        //rhoCv = rhov * (y[0]*(air.Calc_cp_Curve(0,T)-air.Rs[0]) +
        //                y[1]*(air.Calc_cp_Curve(1,T)-air.Rs[1]) +
        //                y[2]*(air.Calc_cp_Curve(2,T)-air.Rs[2]));

        hv = y[0]*air.Calc_h_Curve(0,T) + y[1]*air.Calc_h_Curve(1,T) + y[2]*air.Calc_h_Curve(2,T);
        hl = air.Calc_h_Curve(2,T) - air.HVAP;

        h0 = Yv*hv + (1.0-Yv)*hl + 0.5*u*u;                       //Mixture Total Enthaply

        hderivs(air);
        rhoderivs(air);

        a = sqrt(rho_mix*dhT / ( rho_mix*(dhT*drhoP-drhoT*dhp) + drhoT));    //Wavespeed



        ASSERT(!_isnan(a), "Failed to calculate sound speed!")
    }

    void hderivs(Chem& air){
        double dT = unk[5]*1e-8;
        //double dp = unk[3]*1e-8;
        double y[3]{};
        y[0] = unk[0];
        y[1] = unk[1];
        y[2] = 1.0 - y[0] - y[1];

        dhp = 0.0;

        double Tnew = unk[5] + dT;
        double hnew = y[0]*air.Calc_h_Curve(0,Tnew) + y[1]*air.Calc_h_Curve(1,Tnew) + (1.0 + y[2])*air.Calc_h_Curve(2,Tnew) - air.HVAP;
        dhT = (hnew - (hv+hl)) / dT;

        dhp = 0.0;

    }

    void rhoderivs(Chem& air){
        double dT = unk[5]*1e-8;
        double dp = unk[3]*1e-8;
        double Tnew = unk[5] + dT;
        double p = unk[3];

        double y[3]{};
        y[0] = unk[0];
        y[1] = unk[1];
        y[2] = 1.0 - y[0] - y[1];

        double Yv = unk[2];

        double rholnew = rhol;  //air.RHOLREF + (p - air.PREF)/(air.AREF*air.AREF);
        double rhovnew = rhov*unk[5] / (Tnew)      ;//   p/(Tnew* (y[0]*air.Rs[0] + y[1]*air.Rs[1] + y[2]*air.Rs[2])));
        double rho_mixnew = 1.0 / ((Yv/rhovnew) + ((1.0-Yv)/rholnew));

        drhoT = (rho_mixnew - rho_mix) / dT;

        double pnew = p+dp;
        rholnew = air.RHOLREF + (pnew - air.PREF)/(air.AREF*air.AREF);
        rhovnew = pnew / (unk[5]*(y[0]*air.Rs[0] + y[1]*air.Rs[1] + y[2]*air.Rs[2]));
        rho_mixnew = 1.0 / ((Yv/rhovnew) + ((1.0-Yv)/rholnew));
        drhoP = (rho_mixnew - rho_mix) / dp;


        /// INCORRECT drhoP = -rho_mix*rho_mix*((-Yv/rhov*rhov)*(rhov/p) + (1.0-Yv)/(-rhol*rhol*air.AREF*air.AREF) );
    }

};

#endif //PROJECT2_STATEVARIABLES_H
