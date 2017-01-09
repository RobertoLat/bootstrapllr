#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "bootstrap.h"
#include "observables.h"
#include "logint.cpp"
using namespace std;


void observables::generate_functions(long double *beta,bootstrap* boots){
  
  this->beta=*beta;

  alglib::real_1d_array AX, AY;
  AX.setlength(intervals);
  AY.setlength(intervals);
  for(int i=0; i < intervals; i++){
    AX[i]=boots->data[i].E0;
    AY[i]=boots->aval[i];
  }

  alglib::spline1dbuildcubic(AX,AY,aint);
  Emin=boots->data[0].E0;
  Emax=boots->data[intervals-1].E0;
  
  
};



void observables::generate_functions(bootstrap* boots){

  alglib::real_1d_array AX, AY;
  AX.setlength(intervals);
  AY.setlength(intervals);
  for(int i=0; i < intervals; i++){
    AX[i]=boots->data[i].E0;
    AY[i]=boots->aval[i];
  }

  alglib::spline1dbuildcubic(AX,AY,aint);
  Emin=boots->data[0].E0;
  Emax=boots->data[intervals-1].E0;
  
  
};

double observables::logrho(double x){return alglib::spline1dintegrate(this->aint,x);}
//double observables::operator()(double x) {return logrho(x);}
void observables::setbeta(double x){this->beta=x;}
//virtual double observables::operator()(double x){return 0;}


double Z::ret(double x){return logrho(x)+this->beta*x;}
//double Z::operator()(double x) {return this->ret(x);}

double E::ret(double x) {return logrho(x)+this->beta*x+log(x);}
//double E::operator()(double x) {return this->ret(x);}

double ES::ret(double x) {return logrho(x)+this->beta*x+2.*log(x);}
//double ES::operator()(double x) {return this->ret(x);}
/*
double Z::operator()(double x) {return logrho(x)+this->beta*x;}
double E::operator()(double x) {return logrho(x)+this->beta*x+log(x);}
double ES::operator()(double x) {return logrho(x)+this->beta*x+2.*log(x);}
*/
