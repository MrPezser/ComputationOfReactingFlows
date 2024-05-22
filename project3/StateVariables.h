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
    double rho_mix{},a{},h0{}, mu{}, rhol{}, rhov{}, hv{}, hl{}, dhT{}, dhp{},
    drhoT{}, drhoP{}, dp{}, drhoyvn2{}, drhoyvo2{}, dhyvo2{}, dhyvn2{};
    double rhoR{},rhoCv{};

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

        rhoR = p/T;
        rhoCv = rhov * (y[0]*(air.Calc_cp_Curve(0,T)-air.Rs[0]) +
                        y[1]*(air.Calc_cp_Curve(1,T)-air.Rs[1]) +
                        y[2]*(air.Calc_cp_Curve(2,T)-air.Rs[2]));

        hv = y[0]*air.Calc_h_Curve(0,T) + y[1]*air.Calc_h_Curve(1,T) + y[2]*air.Calc_h_Curve(2,T);
        hl = air.Calc_h_Curve(2,T) - air.HVAP;

        h0 = Yv*hv + (1.0-Yv)*hl + 0.5*u*u;                       //Mixture Total Enthaply

        hderivs(air);
        rhoderivs(air);

        a =sqrt(rho_mix*dhT / ( rho_mix*(dhT*drhoP-drhoT*dhp) + drhoT));    //Wavespeed
            // sqrt((p/rho_mix)*(1.0 + rhoR/rhoCv));

        dp = cbrt( 6.0*(1.0 - unk[2])/(M_PI*rhol*unk[6]) );
        dp = fmax(0.0,dp);


        mu = (1.716e-5) * pow(T/273.15, 3.0/2.0) * ((273.15 + 110.4)/(T + 110.4));

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

        double Yv = unk[2];
        double Yl = 1.0 - Yv;

        double Tnew = unk[5] + dT;
        double hnew = Yv*(y[0]*air.Calc_h_Curve(0,Tnew) + y[1]*air.Calc_h_Curve(1,Tnew) + y[2]*air.Calc_h_Curve(2,Tnew))
                        + Yl*(air.Calc_h_Curve(2,Tnew) - air.HVAP);
        dhT = (hnew - (Yv*hv+Yl*hl)) / dT;

        //Derivatives with respect to vapor species fractions
        double dyvo2, dyvn2, yvo2new, yvn2new;
        dyvo2 = unk[0]*1e-8;
        dyvn2 = unk[1]*1e-8;
        yvo2new = unk[0] + dyvo2;
        yvn2new = unk[1] + dyvn2;

        y[0] = yvo2new;
        y[1] = unk[1];
        y[2] = 1.0 - yvo2new - unk[1];
        hnew = Yv*(y[0]*air.Calc_h_Curve(0,unk[5]) + y[1]*air.Calc_h_Curve(1,unk[5]) + y[2]*air.Calc_h_Curve(2,unk[5]))
               + Yl*(air.Calc_h_Curve(2,unk[5]) - air.HVAP);
        dhyvo2 = (hnew - (Yv*hv+Yl*hl)) / dyvo2;

        y[0] = unk[0];
        y[1] = yvn2new;
        y[2] = 1.0 - yvn2new - unk[0];
        hnew = Yv*(y[0]*air.Calc_h_Curve(0,unk[5]) + y[1]*air.Calc_h_Curve(1,unk[5]) + y[2]*air.Calc_h_Curve(2,unk[5]))
               + Yl*(air.Calc_h_Curve(2,unk[5]) - air.HVAP);
        dhyvo2 = (hnew - (Yv*hv+Yl*hl)) / dyvn2;

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

        //Derivatives with respect to vapor species fractions
        double dyvo2, dyvn2, yvo2new, yvn2new;
        dyvo2 = unk[0]*1e-8;
        dyvn2 = unk[1]*1e-8;
        yvo2new = unk[0] + dyvo2;
        yvn2new = unk[1] + dyvn2;

        y[0] = yvo2new;
        y[1] = unk[1];
        y[2] = 1.0 - yvo2new - unk[1];
        rholnew = air.RHOLREF + (p - air.PREF)/(air.AREF*air.AREF);
        rhovnew = p / (unk[5]*(y[0]*air.Rs[0] + y[1]*air.Rs[1] + y[2]*air.Rs[2]));
        rho_mixnew = 1.0 / ((Yv/rhovnew) + ((1.0-Yv)/rholnew));
        drhoyvo2 = (rho_mixnew - rho_mix) / dyvo2;

        y[0] = unk[0];
        y[1] = yvn2new;
        y[2] = 1.0 - yvn2new - unk[0];
        rholnew = air.RHOLREF + (p - air.PREF)/(air.AREF*air.AREF);
        rhovnew = p / (unk[5]*(y[0]*air.Rs[0] + y[1]*air.Rs[1] + y[2]*air.Rs[2]));
        rho_mixnew = 1.0 / ((Yv/rhovnew) + ((1.0-Yv)/rholnew));
        drhoyvn2 = (rho_mixnew - rho_mix) / dyvn2;

    }

};

#endif //PROJECT2_STATEVARIABLES_H
