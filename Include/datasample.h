#ifndef DATASAMPLE_H
#define DATASAMPLE_H

#include <vector>
#include <iostream>


using std::vector;
using std::ostream;


struct estimate {
  long double val;
  long double err;

  friend ostream &operator<<(ostream &out, estimate e){
    return out<<e.val<<" +- "<<e.err;
  }
};



class datasample: public vector<long double> {
public:
  datasample(){}
  ~datasample(){}

  long double E0,dE,aval;

  estimate avr() const;
  estimate stderr() const;
  estimate autocorr() const;

  void ci(long double low, long double med, long double hi, long double *il, long double *im, long double *ih) const;
  void hist(long double *n, int npt, long double *x, long double w) const;
  void cumhist(long double *n, int npt, long double *x) const;

  estimate JKavr(const int blsize) const;
  estimate JKvar(const int blsize) const;

  friend estimate JKcov(const datasample &dt1, const datasample &dt2, const int blsize);

};




#endif //#ifndef DATASAMPLE_H
