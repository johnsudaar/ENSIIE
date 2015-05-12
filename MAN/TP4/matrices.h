#ifndef _MATRICES_H
#define _MATRICES_H

#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

struct matrice{
  double** values;
  int size_x;
  int size_y;
};

typedef struct matrice matrice;


matrice* constructMatrice(int s_x, int s_y, double def);
matrice* id(int k);
double get(int x, int y, matrice* m);
int set(int x, int y, double v, matrice* m);

double getMath(int y, int x, matrice* m);
int setMath(int y, int x, double v, matrice* m);

int isTriangInf(matrice* m);
int isTriangSup(matrice* m);
int isSquare(matrice* m);
matrice* transpose(matrice* A);
matrice* multiply(matrice* A, matrice* B);
matrice* sum(matrice* A, matrice* B);
matrice* diff(matrice* A, matrice* B);
matrice** decomp(matrice* A);

void printMatrice(matrice* m);
void resetMatrice(matrice *m, double val);
matrice* copy(matrice* base);
void freeMatrice(matrice* base);
matrice* hilbert(int n);

matrice* loadFromFile(char* file);
#endif
