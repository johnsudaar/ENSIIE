#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"
#include "matrice_computing.h"

matrice* generateMatrice(int n){
  matrice *m = constructMatrice(n,n,0);
  int i;
  for(i=0;i<n; i++){
    set(i,i,3,m);
    if(i!= 0)
      set(i,i-1,1,m);
    if(i!= n-1)
      set(i,i+1,1,m);
  }
  return m;
}

void computeTimes(){
  int i;
  clock_t t;
  matrice* A;
  matrice* AA;
  matrice* b;
  matrice* temp;
  double time_taken;
  for(i=150; i<=250; i+=10){
    printf("I = %d\n",i);
    A = generateMatrice(i);
    AA = multiply(A,A);
    b = constructMatrice(1,i,1);
    t = clock();
    resolutionLU(AA,b);
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Temps de resolution A*A (LU) : %f\n",time_taken);

    t = clock();
    temp = resolutionLU(A,b);
    resolutionLU(A,temp);
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Temps de resolution A (LU) : %f\n",time_taken);
    freeMatrice(temp);

    t = clock();
    ResolutionCholesky(AA,b);
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Temps de resolution A (Cholesky) : %f\n",time_taken);

    t = clock();
    temp = ResolutionCholesky(A,b);
    ResolutionCholesky(A,temp);
    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Temps de resolution A (Cholesky) : %f\n",time_taken);
    freeMatrice(temp);
    freeMatrice(AA);
    freeMatrice(A);
    freeMatrice(b);
  }
}


int main(){

  // ******************************************************** PARTIE 1 A ***************************
  /*printf("**************** 32 ***************\n");
  matrice *A = loadFromFile("Matrices/32.mat");

  printf("** A **\n");
  printMatrice(A);
  matrice* B = Cholesky(A);
  printf("** B **\n");
  printMatrice(B);

  freeMatrice(B);
  printf("**************** 33 ***************\n");

  matrice* b = loadFromFile("Matrices/33.vect");

  printf("** A **\n");
  printMatrice(A);
  printf("** b **\n");
  printMatrice(b);
  B = ResolutionCholesky(A,b);
  printf("** x **\n");
  printMatrice(B);

  freeMatrice(b);
  freeMatrice(A);
  freeMatrice(B);

  printf("**************** 4 *****************\n");
  A = generateMatrice(3);
  printMatrice(A);

  B = multiply(A,A);
  printMatrice(B);*/

  computeTimes();

  return 1;
}
