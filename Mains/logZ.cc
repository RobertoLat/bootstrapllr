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
  cout<<"#energy window maximum probability is at "<<max<<endl;
  double emin=bisection(*(y),(*y)(max)-30.,y->Emin,max,0.0000001);
  cout<<"#energy window emin is at "<<emin<<endl;
  double emax=bisection(*(y),(*y)(max)-30.,max,y->Emax,0.0000001);
  cout<<"#energy window emax is at "<<emax<<endl;
}


int main(int argc, char* argv[]){

  inputdata mydata;

  mydata.read_cmdline(argc, argv);

  bootstrap * mybs=mydata.populate_bs_sample(mydata.binsize);


  cvmax cvm(mybs->intervals);  
  Z myz(mybs->intervals);
  E mye(mybs->intervals);
  //double betac=0.,dbetac=0.;
  


  long double beta= mydata.beta;
  long double L= mydata.L;

  double lZ=0.,dlZ=0.;
  double lE=0.,dlE=0.;
  double cv=0.,dcv=0.;
  

  cout << setprecision(10)<< "# beta " << beta<< " L " << L << endl; 


  for( int j=0; j< mydata.nboots ; j++) {
    
    for( int i=0;i<mybs->intervals; i++)
      mybs->generate_sample(i);
    
    
    myz.generate_functions(&beta,mybs);
    mye.generate_functions(&beta,mybs);
    cvm.generate_functions(mybs);
    
    double temp;
    if(j==0){
      observables *ch=&myz;
      check(ch);
      energywindow(ch);
      ch=NULL;
    }
    temp=logint(myz,myz.Emin,myz.Emax);
    lZ+=temp/static_cast<double>(mydata.nboots);
    dlZ+=temp*temp/static_cast<double>(mydata.nboots);
    double temp2=exp(logint(mye,mye.Emin,mye.Emax)-temp);
    lE+=temp2/static_cast<double>(mydata.nboots);
    dlE+=temp2*temp2/static_cast<double>(mydata.nboots);

    
    
    
    temp=cvm(beta);
    cv+=temp/static_cast<double>(mydata.nboots);
    dcv+=temp*temp/static_cast<double>(mydata.nboots);
    
    //temp=findmax(cvm,0.2,2.5);
    //betac+=temp/static_cast<double>(mydata.nboots);
    //dbetac+=temp*temp/static_cast<double>(mydata.nboots);



  }

  cout<<"logZ = "<<lZ<<" "<<sqrt(dlZ-lZ*lZ)<<endl;
  cout<<"E = "<<lE/(static_cast<double>(L)*static_cast<double>(L)*static_cast<double>(L)*24.)<<" "<<sqrt(dlE-lE*lE)/(static_cast<double>(L)*static_cast<double>(L)*static_cast<double>(L)*24.)<<endl;
  cout<<"cv = "<<cv/(static_cast<double>(L)*static_cast<double>(L)*static_cast<double>(L)*24.)<<" "<<sqrt(dcv-cv*cv)/(static_cast<double>(L)*static_cast<double>(L)*static_cast<double>(L)*24.)<<endl;



  
  return 0;
}
