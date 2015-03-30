#include "matrices.h"

matrice* constructMatrice(int s_x, int s_y, double def){
  int x,y;
  matrice* cur = (matrice*) malloc(sizeof(matrice));

  cur->values = (double**) malloc(s_x*sizeof(double*));

  for(x=0 ; x<s_x ; x++){
    cur->values[x] = (double*) malloc(sizeof(double));
    for(y=0; y<s_y; y++){
      cur->values[x][y] = def;
    }
  }
  cur->size_x = s_x;
  cur->size_y = s_y;
  return cur;
}

int isPosValide(int x, int y, matrice* m){
  return x < m->size_x && y < m->size_y;
}

double get(int x, int y, matrice* m){
  if(isPosValide(x,y,m)){
    return m->values[x][y];
  }else{
    fprintf(stderr, "[ERREUR] ACCES A UNE CASE INVALIDE X=%d Y=%d S_X=%d S_Y=%d\n",x,y,m->size_x,m->size_y);
    return 0;
  }
}

int set(int x, int y, double v, matrice* m){
  if(isPosValide(x,y,m)){
    m->values[x][y] = v;
    return 1;
  }else{
    fprintf(stderr, "[ERREUR] ACCES A UNE CASE INVALIDE X=%d Y=%d S_X=%d S_Y=%d\n",x,y,m->size_x,m->size_y);
    return 0;
  }
}

double getMath(int y, int x, matrice* m){
  return get(x-1,y-1,m);
}
int setMath(int y, int x, double v, matrice* m) {
  return set(x-1,y-1,v,m);
}

int isTriangInf(matrice* m){
  if(! isSquare(m))
    return 0;
  int x,y;
  for(x=0 ; x<m->size_x ; x++){
    for(y=0 ; y<m->size_y ; y++){
      if(x>y && get(x,y,m) != 0){
        return 0;
      }
    }
  }
  return 1;
}

int isTriangSup(matrice* m){
  if(! isSquare(m))
    return 0;
  int x,y;
  for(x=0 ; x<m->size_x ; x++){
    for(y=0 ; y<m->size_y ; y++){
      if(x<y && get(x,y,m) != 0){
        return 0;
      }
    }
  }
  return 1;
}


int isSquare(matrice* m){
  return m->size_x == m->size_y;
}

void printMatrice(matrice* m){
  printf("M(%d,%d) :\n",m->size_x, m->size_y);
  int x,y;
  for(y=0; y<m->size_y; y++){
    for(x=0; x<m->size_x; x++){
      printf(" %.3f ",m->values[x][y]);
      if(m->values[x][y] >= 0)
        printf(" ");
      if(x != m->size_x - 1)
        printf("|");
    }
    printf("\n");
    if(y != m->size_y - 1){
      for(x=0; x<8*m->size_x; x++){
        printf("-");
      }
    }
    printf("\n");
  }
}

void resetMatrice(matrice *m, double val){
  int x,y;
  for(x=0; x<m->size_x; x++){
    for(y=0; y<m->size_y;y++){
      m->values[x][y] = val;
    }
  }
}

matrice* copy(matrice* base){
  matrice* new = constructMatrice(base->size_x, base->size_y,0);
  int x,y;
  for(x = 0; x<base->size_x; x++){
    for(y = 0; y<base->size_y; y++){
      new->values[x][y] = base->values[x][y];
    }
  }
  return new;
}

void freeMatrice(matrice* base){
  int x;
  for(x = 0; x<base->size_x; x++){
    free(base->values[x]);
  }
  free(base->values);
  free(base);
}

matrice* id(int k){
  matrice* res = constructMatrice(k,k,0);
  int i;
  for(i = 0; i<k;i++){
    res->values[i][i] = 1;
  }
  return res;
}

matrice* loadFromFile(char* file){
  FILE* fichier = fopen(file,"r");
  if(fichier == NULL){
    perror(file);
    exit(EXIT_FAILURE);
  }
  int s_x, s_y,x,y;
  double cur;
  fscanf(fichier, "%d %d",&s_x, &s_y);
  getLine(fichier);
  matrice * res = constructMatrice(s_x,s_y,0);
  for(y =0; y<s_y; y++){
    for(x = 0; x<s_x; x++){
      fscanf(fichier, "%lf", &cur);
      set(x,y,cur,res);
    }
    getLine(fichier);
  }
  fclose(fichier);
  return res;
}
