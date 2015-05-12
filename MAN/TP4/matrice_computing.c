#include "matrice_computing.h"

matrice* solinf(matrice* L, matrice* b){
  if(! isTriangInf(L)){
    fprintf(stderr,"Il se peut que vous ayez fait une erreure, j'ai l'impression que la matrice n'est pas triangulaire inférieure ...\n");
    return NULL;
  }
  if(L->size_y != b->size_y){
    fprintf(stderr,"Il se peut que vous ayez fait une erreure, j'ai l'impression que le vecteur n'as pas la bonne taille ...\n");
  }

  double sum;
  int x,y;
  matrice* res = constructMatrice(1,L->size_y,0);
  for(y=0; y<L->size_y; y++){
    sum = 0;
    for(x= 0 ;x < y; x++){
      sum += L->values[x][y] * res->values[0][x];
    }
    res->values[0][y] = (b->values[0][y] - sum) / L->values[y][y];
  }
  return res;
}

matrice* solsup(matrice* L, matrice* b){
  if(! isTriangSup(L)){
    fprintf(stderr,"Il se peut que vous ayez fait une erreure, j'ai l'impression que la matrice n'est pas triangulaire supérieure...\n");
    return NULL;
  }
  if(L->size_y != b->size_y){
    fprintf(stderr,"Il se peut que vous ayez fait une erreure, j'ai l'impression que le vecteur n'as pas la bonne taille ...\n");
  }

  double sum;
  int x,y;
  matrice* res = constructMatrice(1,L->size_y,0);
  for(y=L->size_y - 1; y>=0; y--){
    sum = 0;
    for(x= y+1 ; x < L->size_x; x++){
      sum += L->values[x][y] * res->values[0][x];
    }
    res->values[0][y] = (b->values[0][y] - sum) / L->values[y][y];
  }
  return res;
}

matrice** trigGauss(matrice* _A, matrice* _b){
  if( ! isSquare(_A))
    fprintf(stderr,"Il se peut que vous ayez fait une erreure, j'ai l'impression que la matrice n'est pas carré...\n");
  if( _A->size_y != _b->size_y)
    fprintf(stderr,"Il se peut que vous ayez fait une erreure, j'ai l'impression que le vecteur n'as pas la bonne taille ...\n");

  matrice* A = copy(_A);
  matrice* b = copy(_b);
  double eps = 0.001;

  int k,i,j;

  double c;
  for(k = 0; k<A->size_x-1; k++){
    if(abs(A->values[k][k]) < eps){
      fprintf(stderr,"Pivot null...");
      return NULL;
    }else{
      for(i=k+1 ; i < A->size_x ; i++ ){
        c = A->values[k][i] / A->values[k][k];
        b->values[0][i] = b->values[0][i] - c * b->values[0][k];
        A->values[k][i] = 0;
        for(j =k+1; j<A->size_x; j++){
          A->values[j][i] = A->values[j][i] - c * A->values[j][k];
        }
      }
    }
  }
  matrice** res = (matrice**) malloc(2*sizeof(matrice*));
  res[0] = A;
  res[1] = b;
  return res;
}


matrice*  ResolutionGauss(matrice* A, matrice* b){
  matrice** trig = trigGauss(A,b);
  matrice * res = solsup(trig[0],trig[1]);
  freeMatrice(trig[0]);
  freeMatrice(trig[1]);
  return res;
}

matrice** FactLU(matrice* A){
  int n = A->size_x;
  matrice* L = id(n);
  matrice* U = constructMatrice(n,n,0);

  int i,j,k;
  double sum;
  for(i=1; i<= n -1; i++){
    for(j = i; j<= n ; j++){
      sum = 0;
      for(k = 1; k<= i -1; k++){
        sum += getMath(i,k,L)*getMath(k,j,U);
      }
      sum = getMath(i,j,A) - sum;
      setMath(i,j,sum,U);
    }
    for(j=i+1; j<=n; j++){
      sum = 0;
      for(k=1; k<= i-1; k++){
        sum += getMath(j,k,L) * getMath(k,i,U);
      }
      sum = getMath(j,i,A) - sum;
      setMath(j,i,sum/getMath(i,i,U),L);
    }
  }
  sum = 0;
  for(k=1; k<=n-1; k++){
    sum += getMath(n,k,L) * getMath(k,n,U);
  }

  sum = getMath(n,n,A) - sum;
  setMath(n,n,sum,U);

  matrice ** res = (matrice**) malloc(2* sizeof(matrice*));
  res[0] = L;
  res[1] = U;
  return res;
}

matrice*  resolutionLU(matrice* A, matrice* b){
  matrice **facto = FactLU(A);
  matrice *y = solinf(facto[0],b);
  matrice *res = solsup(facto[1],y);

  freeMatrice(facto[0]);
  freeMatrice(facto[1]);
  freeMatrice(y);
  return res;
}

matrice* resolutionLU2(matrice** facto, matrice* b){
  matrice* y = solinf(facto[0],b);
  matrice* res = solsup(facto[1],y);
  freeMatrice(y);
  return res;
}

matrice* inverse(matrice* A){
  if( ! isSquare(A)){
    fprintf(stderr, "Matrice pas carrée !");
    return NULL;
  }
  else{
    int n = A->size_x;
    int i,j;
    matrice** facto = FactLU(A);
    matrice* b = constructMatrice(1,A->size_x,0);
    matrice* res;
    matrice* final = constructMatrice(A->size_x, A->size_y, 0);
    for(i=0; i<n; i++){
      // --- Génération du vecteur pour la resolution ---
      resetMatrice(b,0);
      set(0,i,1,b);

      // --- Resolution du système ...
      res = resolutionLU2(facto,b);

      // --- Recopie dans la matrice ---

      for(j=0; j<A->size_y;j++){
        set(i,j,get(0,j,res),final);
      }
      freeMatrice(res);
    }
    freeMatrice(facto[0]);
    freeMatrice(facto[1]);
    free(facto);
    freeMatrice(b);
    return final;
  }
}

matrice* Cholesky(matrice* B){
  matrice* A = copy(B);
  int n = A->size_x;
  int j,k,i;
  double res;
  for(j = 1; j<=n; j++){
    for(k = 1; k<= j-1; k++){
      res = getMath(j,j,A) - pow(getMath(j,k,A),2);
      setMath(j,j,res,A);
    }
    //setMath(j,j,sqrt(getMath(j,j,A)),A);
    for(i=j+1; i<=n; i++){
      for(k = 1; k<= j - 1; k++){
        res = getMath(i,j,A) - getMath(i,k,A)*getMath(j,k,A);
        setMath(i,j,res,A);
      }
      res = getMath(i,j,A) / getMath(j,j,A);
      setMath(i,j,res,A);
    }
  }
  for(j = 0; j<n; j++){
    for(i = j+1; i<n; i++){
      set(i,j,0,A);
    }
  }
  return A;
}

matrice* ResolutionCholesky(matrice* A, matrice* b){
  matrice *B = Cholesky(A);
  matrice *Bt = transpose(B);

  matrice *y = solinf(B,b);
  matrice *res = solsup(Bt, y);
  freeMatrice(B);
  freeMatrice(Bt);
  freeMatrice(y);
  return res;
}

double norm1(matrice* A){
  int i,j;
  double max, cur;
  for(j = 1; j<= A->size_x; j++){
    cur = 0;
    for(i = 1; i<= A->size_y; i++){
      cur += absD(getMath(i,j,A));
    }
    if(j == 1 || cur > max){
      max = cur;
    }
  }
  return max;
}

double cond(matrice* A){
  matrice* B = inverse(A);
  double res = norm1(A) * norm1(B);
  freeMatrice(B);
  return res;
}


matrice* iter(matrice *B, matrice* b, int n, matrice* x0){
  // x(k+1) = B*x(k) +b
  if(n <= 0){
    return copy(x0);
  }else{
    matrice* mult = multiply(B,x0);
    matrice* summ = sum(mult,b);
    freeMatrice(mult);
    matrice *res = iter(B,b,n-1,summ);
    freeMatrice(summ);
    return res;
  }
}

matrice** jac(matrice* A, matrice* b){
  matrice** dec = decomp(A);
  matrice*  d1  = inverse(dec[0]);

  matrice*  p2 = sum(dec[1], dec[2]);
  matrice*  J  = multiply(d1,p2);
  matrice*  d  = multiply(d1,b);
  freeMatrice(dec[0]);
  freeMatrice(dec[1]);
  freeMatrice(dec[2]);
  free(dec);
  freeMatrice(p2);
  freeMatrice(d1);
  matrice **res = (matrice**)malloc(2 * sizeof(matrice*));
  res[0] = J;
  res[1] = d;
  return res;
}

matrice** GS(matrice* A, matrice* b){
  matrice** dec = decomp(A);
  matrice*  DE  = diff(dec[0],dec[1]);
  matrice*  DE1 = inverse(DE);
  matrice*  G   = multiply(DE1,dec[2]);
  matrice*  b1  = multiply(DE1,b);
  matrice** res = (matrice**)malloc(2*sizeof(matrice*));
  res[0] = G;
  res[1] = b1;
  return res;
}

matrice* iterjac(matrice* A, matrice*b , int n){
  matrice** jacs = jac(A,b);
  matrice* res  = iter(jacs[0],jacs[1],n,jacs[1]);

  freeMatrice(jacs[0]);
  freeMatrice(jacs[1]);
  free(jacs);

  return res;
}

matrice* iterGS(matrice* A, matrice*b , int n){
  matrice** gss = GS(A,b);

  matrice* res  = iter(gss[0],gss[1],n,gss[1]);

  freeMatrice(gss[0]);
  freeMatrice(gss[1]);
  free(gss);

  return res;
}
