#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"
#include "matrice_computing.h"

int main(){
  matrice *A = loadFromFile("Matrices/A.mat");
  matrice *b = loadFromFile("Matrices/b.vect");

  printf("====== 1) b) ======\n");

  matrice *Id = id(4);
  matrice *IdA = diff(Id, A);
  printf("Id - A )\n");
  printMatrice(IdA);
  printf("b )\n");
  printMatrice(b);
  printf("Res )\n");
  matrice* res = iter(IdA,b,20,b);
  printMatrice(res);
  freeMatrice(res);

  printf("===== 2) b) ======\n");
  matrice** res2 = jac(A,b);
  printf("J )\n");
  printMatrice(res2[0]);
  printf("b )\n");
  printMatrice(res2[1]);

  freeMatrice(res2[0]);
  freeMatrice(res2[1]);
  free(res2);

  res2 = GS(A,b);
  printf("G )\n");
  printMatrice(res2[0]);
  printf("b )\n");
  printMatrice(res2[1]);

  freeMatrice(res2[0]);
  freeMatrice(res2[1]);
  free(res2);

  printf("==== 2) c) ======\n");
  int i;
  for(i=10; i<=30; i+=10){
    printf(" === < %d > === \n",i);
    res = iterjac(A,b,i);
    printf("Jacobi : \n");
    printMatrice(res);
    freeMatrice(res);
    res = iterGS(A,b,i);
    printf("Gauss Seidel : \n");
    printMatrice(res);
    freeMatrice(res);
  }

  freeMatrice(A);
  freeMatrice(b);
  freeMatrice(Id);
  freeMatrice(IdA);
  return 1;
}
