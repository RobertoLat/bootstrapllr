#ifndef BOOSTRAP_H
#define BOOSTRAP_H
#include "datasample.h"
#include "ranlxs.h"
#include <iostream>
using namespace std;

class bootstrap {
public:
  int intervals;
  datasample* data;
  datasample* bs;
  long double * aval;

  bootstrap(int ivls,int seed) { 
    intervals = ivls; 
    data = new datasample[ivls]; 
    bs= new datasample;  
    aval = new long double[ivls];
    rlxs_init(2,seed);
  }
  ~bootstrap() { delete[] data; delete[] bs;  delete[] aval; }
  void purge() { for(int i=0;i<intervals;i++) data[i].clear();}
  void generate_sample(int intv, int elen);
  void generate_sample(int intv);

};

#endif
