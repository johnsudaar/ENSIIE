#include "utils.h"

void getLine(FILE * fichier){
  char cur;
  do{
    cur = fgetc(fichier);
  }while(cur != '\n' && cur != EOF);
}

double absD(double v){
  return v > 0? v : (-1.0)*v;
}
