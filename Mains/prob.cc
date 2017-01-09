#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <functional>
#include "bootstrap.h"
#include "datareader.h"
#include "observables.h"
#include "logint.cpp"
#include "integration.h"
#include "optimization.h"
using namespace std;

class cvmax{
public:
  cvmax(int n):e(n),z(n),es(n){}
  E e;
  Z z;
  ES es;
  double temp,temp2;
  void generate_functions(bootstrap *mybs){
    z.generate_functions(mybs);
    e.generate_functions(mybs);
    es.generate_functions(mybs);

  }
  double operator()(double x){
    e.setbeta(x);
    z.setbeta(x);
    es.setbeta(x);
    temp=logint(z,z.Emin,z.Emax);
    temp2=exp(logint(e,e.Emin,e.Emax)-temp);
    return exp(logint(es,es.Emin,es.Emax)-temp)-temp2*temp2;

  }

};

void check(observables *y){
  double temp;
  temp=findmax(*(y),y->Emin,y->Emax);
  cout<<"#check at Emin the integrand is suppresed "<<exp((*y)(y->Emin)-(*y)(temp))<<" times with respect to the maximum"<<endl;
  cout<<"#check at Emax the integrand is suppresed "<<exp((*y)(y->Emax)-(*y)(temp))<<" times with respect to the maximum"<<endl;
}

void energywindow(observables *y){
  double max;
  max=findmax(*(y),y->Emin,y->Emax);
  cout<<"maximum probability is at "<<max<<endl;
  double emin=bisection(*(y),(*y)(max)-30.,y->Emin,max,0.0000001);
  cout<<"emin is at "<<emin<<endl;
  double emax=bisection(*(y),(*y)(max)-30.,max,y->Emax,0.0000001);
  cout<<"emax is at "<<emax<<endl;
}


int main(int argc, char* argv[]){

  inputdata mydata;

  mydata.read_cmdline(argc, argv);

  bootstrap * mybs=mydata.populate_bs_sample(mydata.binsize);


  cvmax cvm(mybs->intervals);  
  Z myz(mybs->intervals);
  E mye(mybs->intervals);
  //double betac=0.,dbetac=0.;
  

  double *a;
  double *da;
 
  a=new double[mybs->intervals];
  da=new double[mybs->intervals];

  for( int i=0;i<mybs->intervals; i++){
    a[i]=0;
    da[i]=0;

  }

  long double beta= mydata.beta;
  long double L= mydata.L;  

  cout << setprecision(10)<< "# beta " << beta<< " L " << L << endl; 


  for( int j=0; j< mydata.nboots ; j++) {
    
    for( int i=0;i<mybs->intervals; i++)
      mybs->generate_sample(i);
    double temp1=0;
    double lognorm=0;
    double logprob;
    for( int i=0;i<mybs->intervals; i++){
      logprob=(temp1+beta*mybs->data[i].E0);
      lognorm=logsum(lognorm,logprob+log(mybs->data[1].E0-mybs->data[0].E0));
      temp1+=mybs->aval[i]*(mybs->data[1].E0-mybs->data[0].E0);
    }

    temp1=0;

    for( int i=0;i<mybs->intervals; i++){
      a[i]+=(temp1+beta*mybs->data[i].E0-lognorm)/static_cast<double>(mydata.nboots);
      da[i]+=(temp1+beta*mybs->data[i].E0-lognorm)*(temp1+beta*mybs->data[i].E0-lognorm)/static_cast<double>(mydata.nboots);
      temp1+=mybs->aval[i]*(mybs->data[1].E0-mybs->data[0].E0);
	}
  }

    for( int i=0;i<mybs->intervals; i++){
      cout<<mybs->data[i].E0<<" "<<a[i]<<" "<<sqrt(da[i]-a[i]*a[i])<<endl;
  
    }
  //cout<<"betac = "<<betac<<" "<<sqrt(dbetac-betac*betac)<<endl;

  
  return 0;
}
