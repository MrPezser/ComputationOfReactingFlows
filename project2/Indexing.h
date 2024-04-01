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
#define RXTOL 1e-4 //residual drop to activate thermochemical source terms

#define IDETAIL 1 //printout more detailed residual

#endif //PROJECT2_INDEXING_H
