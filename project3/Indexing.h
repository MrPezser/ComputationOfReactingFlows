//
// Created by tskoepli on 3/8/2024.
//

#ifndef PROJECT2_INDEXING_H
#define PROJECT2_INDEXING_H

#include <cmath>

#define NDEGR (1)
#define NUMSP (3)
#define NVAR (8) //(7)
#define NALTF (4)

#define uIJK(ielem,jdegr,kvar)  (((ielem)*NDEGR + (jdegr))*(NVAR) + (kvar))
#define aIJ(ielem,jdegr) ((ielem)*NDEGR + (jdegr))
#define fIJ(ielem,kvar)  ((ielem)*(NVAR) + (kvar))
#define IJ(i,j,nj) ((i)*(nj) + (j))

//#define MXITER 1
#define MXITER (3e4)
#define RESTOL (1e-6) //residual drop to declare convergence
#define RXTOL (1e-3) //residual drop to activate thermochemical source terms
#define CFLTCNE 0.1 //Thermochem CFL factor, need lower timestep to accomodate nonequilibrium effects
#define IDETAIL (1) //printout more detailed residual

#define ASSERT(cond, msg) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n", __FILE__, __LINE__, #cond, msg); abort();}
#define ASSERTVEC(cond, msg, ind) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n Index: %d\n", __FILE__, __LINE__, #cond, msg, ind);abort();}


#endif //PROJECT2_INDEXING_H
