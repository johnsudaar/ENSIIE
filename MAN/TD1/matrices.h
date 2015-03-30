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

void printMatrice(matrice* m);
void resetMatrice(matrice *m, double val);
matrice* copy(matrice* base);
void freeMatrice(matrice* base);

matrice* loadFromFile(char* file);
#endif
