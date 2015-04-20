#ifndef _MATRICE_COMPUTING_H
#define _MATRICE_COMPUTING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrices.h"


matrice*  solinf(matrice* L, matrice* b);
matrice*  solsup(matrice* L, matrice* b);
matrice** trigGauss(matrice* _A, matrice* _b);
matrice*  ResolutionGauss(matrice* A, matrice* b);
matrice** FactLU(matrice* A);
matrice*  resolutionLU(matrice* A, matrice* b);
matrice*  inverse(matrice *A);
matrice*  Cholesky(matrice* A);
matrice*  ResolutionCholesky(matrice* A,matrice* b);
double    norm1(matrice* A);
double    cond(matrice* A);
#endif
