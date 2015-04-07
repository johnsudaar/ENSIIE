#include <stdio.h>
#include <stdlib.h>
#include "matrices.h"
#include "matrice_computing.h"

int main(){

  // ******************************************************** PARTIE 1 A ***************************
  printf("**************** 1A ***************\n");
  matrice *L = loadFromFile("Matrices/1a.mat");
  matrice *b = loadFromFile("Matrices/1a.vect");

  printf("** L **\n");
  printMatrice(L);
  printf("** b **\n");
  printMatrice(b);
  matrice* x = solinf(L,b);
  printf("** x **\n");
  printMatrice(x);

  freeMatrice(L);
  freeMatrice(b);
  freeMatrice(x);

  // ******************************************************** PARTIE 1 B ***************************

  printf("**************** 1B ***************\n");

  L = loadFromFile("Matrices/1b.mat");
  b = loadFromFile("Matrices/1b.vect");

  printf("** U **\n");
  printMatrice(L);
  printf("** b **\n");
  printMatrice(b);
  x = solsup(L,b);
  printf("** x **\n");
  printMatrice(x);

  freeMatrice(L);
  freeMatrice(b);
  freeMatrice(x);

  // ******************************************************** PARTIE 2 A ***************************

  printf("**************** 2A ***************\n");

  L = loadFromFile("Matrices/2a.mat");
  b = loadFromFile("Matrices/2a.vect");

  printf("** A **\n");
  printMatrice(L);
  printf("** b **\n");
  printMatrice(b);
  matrice** res = trigGauss(L,b);
  printf("** At **\n");
  printMatrice(res[0]);
  printf("** bt **\n");
  printMatrice(res[1]);

  freeMatrice(L);
  freeMatrice(b);
  freeMatrice(res[0]);
  freeMatrice(res[1]);
  free(res);

  // ******************************************************** PARTIE 2 B ***************************

  printf("**************** 2B ***************\n");

  L = loadFromFile("Matrices/2b.mat");
  b = loadFromFile("Matrices/2b.vect");

  printf("** A **\n");
  printMatrice(L);
  printf("** b **\n");
  printMatrice(b);
  x = ResolutionGauss(L,b);
  printf("** x **\n");
  printMatrice(x);

  freeMatrice(L);
  freeMatrice(b);
  freeMatrice(x);

  // ******************************************************** PARTIE 3 A ***************************

  printf("**************** 3A ***************\n");

  L = loadFromFile("Matrices/3a.mat");

  printf("** A **\n");
  printMatrice(L);

  res = FactLU(L);
  printf("** L **\n");
  printMatrice(res[0]);
  printf("** U **\n");
  printMatrice(res[1]);

  freeMatrice(L);
  freeMatrice(res[0]);
  freeMatrice(res[1]);
  free(res);


  // ******************************************************** PARTIE 3 B ***************************

  printf("**************** 3B ***************\n");

  L = loadFromFile("Matrices/2b.mat");
  b = loadFromFile("Matrices/2b.vect");

  printf("** A **\n");
  printMatrice(L);
  printf("** b **\n");
  printMatrice(b);
  x = resolutionLU(L,b);
  printf("** x **\n");
  printMatrice(x);

  freeMatrice(L);
  freeMatrice(b);
  freeMatrice(x);

  // ******************************************************** PARTIE 4 ***************************

  printf("**************** 4  ***************\n");

  L = loadFromFile("Matrices/2b.mat");
  printf("** A **\n");
  printMatrice(L);
  x = inverse(L);
  printf("** B **\n");
  printMatrice(x);

  freeMatrice(L);
  freeMatrice(x);
  return 1;
}
