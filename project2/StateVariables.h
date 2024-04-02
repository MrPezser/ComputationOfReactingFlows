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
    double p{},cptr[NSP]{},rho_mix{},rhoR{},rhoCv{},rhoCp{},a{},h0{},hv{},taui[NSP]{},evT[NSP]{},evTv[NSP]{},etr[NSP]{},
            htr[NSP]{}, cpvT[NSP]{}, cpvTv[NSP]{}, Tstar{}, RRtb[NRX]{}, RRlma[NRX]{}, kf[NRX]{}, kb[NRX]{};

    State() = default;

    void Initialize(double* u){
        unk = u;
    }

    void UpdateState(Chem& air, int ireact) {
        p=0;


        rhoR = 0.0;
        rhoCv = 0.0;
        rhoCp = 0.0;
        rho_mix = 0.0;
        h0=0.0;
        hv=0.0;

        for (int isp=0.0; isp<NSP; isp++){
            rho_mix   += unk[isp];
        }

        if (ireact==1) {
            //full tc noneq
            for (int isp = 0; isp < NSP; isp++) {
                cptr[isp] = air.Get_cptr(isp, unk[NSP + 1]);

                rhoR += unk[isp] * air.Ruv / air.Mw[isp];
                rhoCv += unk[isp] * (cptr[isp] - air.Ruv / air.Mw[isp]);
                rhoCp += unk[isp] * (cptr[isp]);
                h0 += unk[isp] * air.Get_htotal(isp, unk[NSP + 1], unk[NSP + 2]) / rho_mix;
                hv += unk[isp] * air.Get_hv(isp, unk[NSP + 1], unk[NSP + 2]) / rho_mix;   //Mixture Vibrational Enthalpy
            }
            Tstar = sqrt(unk[NSP+1]*unk[NSP+2]);
        } else {
            //themochemical equilibrium
            for (int isp = 0; isp < NSP; isp++) {
                cptr[isp] = air.Calc_cp_Curve(isp, unk[NSP + 1]);

                rhoR += unk[isp] * air.Ruv / air.Mw[isp];
                rhoCv += unk[isp] * (cptr[isp] - air.Rs[isp]);
                rhoCp += unk[isp] * (cptr[isp]);
                h0 += unk[isp] * air.Calc_h_Curve(isp, unk[NSP + 1]) / rho_mix;
                hv += 0.0;
            }
            Tstar = unk[NSP+1];
        }


        h0 +=  0.5*unk[NSP]*unk[NSP];                       //Mixture Total Enthaply
        p = rhoR*unk[NSP+1];                                //Presure
        a = sqrt((p/rho_mix)*(1.0 + rhoR/rhoCv));        //Wavespeed

        ASSERT(!_isnan(a), "failed to calculate sound speed")

        CalcEnergies(air, ireact);

        if (ireact==1){CalcRlxTime(air);}
    }

    void CalcRlxTime(Chem& air){
        double T, Tv, fracsum{};
        T = unk[NSP + 1];
        Tv = unk[NSP + 2];

        for (int isp = 0; isp < NSP; isp++) { // NOLINT(modernize-loop-convert)
            double tsri{};
            fracsum = 0.0;

            if (air.thetav[isp] <= 0.0) {
                taui[isp] = -1.0;
                continue;
            }
            for (int jsp = 0; jsp < NSP; jsp++) {
                double x, Aij, Bij, muij;

                muij = air.Mw[isp] * air.Mw[jsp] / (air.Mw[isp] + air.Mw[jsp]);

                Aij = (1.16e-3) * sqrt(muij) * pow(air.thetav[isp], 4.0 / 3.0);
                Bij = 0.015 * pow(muij, 0.25);

                x = unk[jsp] / (air.Mw[jsp] * rho_mix);
                fracsum += x;
                tsri += (x * (p / 101325.0)) / (exp(Aij * (pow(T, -1.0 / 3.0) - Bij) - 18.42));
            }
            taui[isp] = tsri / fracsum;

            /*
            //Park vib relax time
            double sigma, c, n, NA;
            NA = 6.0221408e+23;
            sigma = 3.0e-21*(50000.0/T);
            c = sqrt( 8.0 * air.Ruv * T / (M_PI*air.Mw[isp]) );
            n = rho_mix*NA*fracsum;

            tau[isp] += 1.0 / (sigma * c * n);
            */


            }
    }

    void CalcEnergies(Chem& air, int ireact) {
        double T, Tv;
        if (ireact==1) {
            //Vibration energy at Ttr and Tve
            T = unk[NSP + 1];
            Tv = unk[NSP + 2];
            for (int isp = 0; isp < NSP; isp++) {
                evT[isp] = air.Get_hv(isp, T, T);
                evTv[isp] = air.Get_hv(isp, T, Tv);
                htr[isp] = air.Get_htr(isp, T);
                etr[isp] = htr[isp] - air.Rs[isp] * T;
                cpvT[isp] = air.Get_cpv(isp, T, T);
                cpvTv[isp] = air.Get_cpv(isp, T, Tv);
            }
            return;
        } else {
            //Vibration energy at Ttr and Tve
            T = unk[NSP + 1];
            for (int isp = 0; isp < NSP; isp++) {
                evT[isp] = 0.0;
                evTv[isp] = 0.0;
                htr[isp] = air.Calc_h_Curve(isp, T);
                etr[isp] = htr[isp] - air.Rs[isp] * T;
                cpvT[isp] = 0.0;
                cpvTv[isp] = 0.0;
            }
            return;
        }
    }

};

#endif //PROJECT2_STATEVARIABLES_H
