#include "utils.h"

void getLine(FILE * fichier){
  char cur;
  do{
    cur = fgetc(fichier);
  }while(cur != '\n' && cur != EOF);
}
