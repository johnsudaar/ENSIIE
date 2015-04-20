#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"
#include "matrice_computing.h"

int main(){
  matrice* A;
  matrice* b;
  matrice* t;
  matrice* dA;
  matrice* db;
  matrice* r;

  int taille;

  printf("***** 1.a ****** \n");
  A = hilbert(3);
  printMatrice(A);
  freeMatrice(A);

  printf("***** 1.b ****** \n");

  for(taille = 1; taille <= 10; taille ++ ){
    A = hilbert(taille);
    printf("Conditionnement de Hibert(%2d) =  %f\n",taille,cond(A));
    freeMatrice(A);
  }

  printf("***** 2.a ****** \n");
  A = loadFromFile("Matrices/1.mat");
  b = loadFromFile("Matrices/1.vect");

  t = inverse(A);
  printMatrice(t);
  freeMatrice(t);

  printf("***** 2.b *****\n");

  t = resolutionLU(A,b);
  printMatrice(t);
  freeMatrice(t);

  printf("***** 2.c *****\n");
  dA = loadFromFile("Matrices/d1.mat");
  db = loadFromFile("Matrices/d1.vect");

  printf("Resolution de Ax = b + db : \n");
  r = sum(b,db);
  t = resolutionLU(A,r);
  printMatrice(t);
  freeMatrice(r);
  freeMatrice(t);

  printf("Resolution de (A+dA)x = b : \n");
  r = sum(A,dA);
  t = resolutionLU(r,b);
  printMatrice(t);
  freeMatrice(r);
  freeMatrice(t);

  return 1;
}
