#ifndef __logint
#define __logint
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

static inline double logsum(double a, double b){
  double res;
  
  if(a>b)
    res=a;
  else
    res=b;

  return res+log1p(exp(-std::abs(a-b)));
}

template <class F>
double logint(F &func,double a,double b,double precision=0.0000001) {
  double c0=log(3./8.);
  double c1=log(7./6.);
  double c2=log(23./24.);
  double res=0.;
  double res_old=1.;
  int ni=10;
  int np;
  double dh,ldh;
  std::vector<double> f;
  f.resize(1024);
  while (abs(res-res_old)>precision){
    res_old=res;
    dh=(b-a)/static_cast<double>(ni);
    ldh=log(dh);
    np=ni+1;
    if(np>f.size())
      f.resize(np*16);

    for(int i=0;i<np;i++){
      f[i]=func(i*dh+a)+ldh;
      //      cout<<"f("<<i*dh+a<<") ="<<f[i]<<endl;
    }
    
    res=logsum(c0+f[0],c1+f[1]);
    res=logsum(res,c2+f[2]);
    res=logsum(res,c0+f[np-1]);
    res=logsum(res,c1+f[np-2]);
    res=logsum(res,c2+f[np-3]);
    
    for(int i=3;i<np-3;i++){
      res=logsum(res,f[i]);
    }


    ni*=2;
  }

  return res;
}


template <class F>
double findmax(F &func,double a,double b,double precision=0.0000001) {
  long double xa=a,xb=b,xm=(xa+xb)/2.,xl,xr;
  long double Fa=func(xa),Fb=func(xb),Fm=func(xm),Fl,Fr;
 
  while(abs(xa-xb)>precision){
    xl=0.5*(xa+xm);
    xr=0.5*(xm+xb);
    Fl=func(xl);
    Fr=func(xr);
    std::vector<long double> v;
    
    v.push_back(Fa);
    v.push_back(Fl);
    v.push_back(Fm);
    v.push_back(Fr);
    v.push_back(Fb);

    std::vector<long double>::iterator result;
    result = std::max_element(v.begin(), v.end());
    
    int res=std::distance(v.begin(),result);
    
    if(res==0 or res==1){
      xb=xm;
      xm=xl;
      Fb=Fm;
      Fm=Fl;
    }else if(res==2){
      xa=xl;
      xb=xr;
      Fa=Fl;
      Fb=Fr;
      
    }else{
      xa=xm;
      xm=xr;
      Fa=Fm;
      Fm=Fr;

    }
    
  }
  return xb;
}


template <class F>
double bisection(F &func,double off, double a0, double a1, double tol ){
  double middle;
  for (;;){
    middle = (a0 + a1)/2.0;
    if (fabs(middle - a0) < tol)
      return middle;
    else if ((func(middle)-off)*(func(a0)-off) < 0.0)
      a1 = middle;
    else
      a0 = middle;
  }
}

#endif






