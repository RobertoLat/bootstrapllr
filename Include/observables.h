#ifndef OBSERVABLES_H
#define OBSERVABLES_H
#include "datasample.h"
#include "stdafx.h" 
#include "interpolation.h"


class observables {
 public:
  observables(int invl){ 
    intervals = invl;
  };
  double operator()(double x){return ret(x);};
  virtual double ret(double x){cout<<"#[OBSERVABLES]: error operator () not defined"<<endl;return 0.;};
  void generate_functions(long double* beta,bootstrap* d);
  void generate_functions(bootstrap* d);
  alglib::spline1dinterpolant aint;
  double logrho(double x);
  long double beta;
  double Emax,Emin;
  int intervals;
  void setbeta(double x);
  virtual ~observables(){};
};


class Z : public observables
{
 public:
  Z(int invl)
    :observables(invl)
  {};
  //double operator()(double x);
  virtual double ret(double x);
};

class E:public observables
{
public:
  E(int invl)
    :observables(invl)
  {};
  //double operator()(double x);
  virtual double ret(double x);
};


class ES:public observables
{
public:
  ES(int invl)
    :observables(invl)
  {};
  //double operator()(double x);
  virtual double ret(double x);
};



struct params{
  long double beta;
  long double shift;
  observables * obs;
};


#endif
