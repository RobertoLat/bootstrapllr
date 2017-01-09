#include "bootstrap.h"
#include "ranlxs.h"

void bootstrap::generate_sample(int intv){
  this->generate_sample(intv,this->data[intv].size());
}

void bootstrap::generate_sample(int intv, int elen){
  const int rl=128;
  static float r[rl];

  this->bs->clear();


  for (int k=elen;k>0;k-=rl) {
    int n=rl; if(k<rl) n=k;
    ranlxs(r,rl);
    for (;n>0;) {
      --n;
      this->bs->push_back(this->data[intv][(int)(r[n]*float(this->data[intv].size()))]);
    }
  }

  this->aval[intv]=this->bs->avr().val;

}

