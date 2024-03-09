//
// Created by tskoepli on 3/8/2024.
//

#ifndef PROJECT2_INDEXING_H
#define PROJECT2_INDEXING_H

#define NDEGR 1
#define NVAR 8

#define uIJK(ielemi,jdegr,kvar)  (((ielem)*NDEGR + (jdegr))*NVAR + (kvar))
#define aIJ(ielem,jdegr) ((ielem)*NDEGR + (jdegr))

#endif //PROJECT2_INDEXING_H
