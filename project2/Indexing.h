//
// Created by tskoepli on 3/8/2024.
//

#ifndef PROJECT2_INDEXING_H
#define PROJECT2_INDEXING_H

#define NDEGR 1
#define NSP 5

#define uIJK(ielem,jdegr,kvar)  (((ielem)*NDEGR + (jdegr))*(NSP+3) + (kvar))
#define aIJ(ielem,jdegr) ((ielem)*NDEGR + (jdegr))
#define fIJ(ielem,kvar)  ((ielem)*(NSP+3) + (kvar))
#define IJ(i,j,nj) ((i)*(nj) + (j))

//#define MXITER 1
#define MXITER 1e7
#define RESTOL 1e-6 //residual drop to decare convergence
#define RXTOL 1e-3 //residual drop to activate thermochemical source terms
#define CFLTCNE 0.25 //Thermochem CFL factor, need lower timestep to accomodate nonequilibrium effects
#define IDETAIL 1 //printout more detailed residual

#define ASSERT(cond, msg) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n", __FILE__, __LINE__, #cond, msg); exit(0);}

#endif //PROJECT2_INDEXING_H
