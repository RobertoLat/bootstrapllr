#include "datasample.h"
#include <cmath>
#include <algorithm>
#include <iostream>
using namespace std;

estimate datasample::avr() const {
  const int len = size();
  estimate res = { 0. , 0. };

  if(len==0) return res; //Check for zero len

  for(int i=0; i<len; i++) {
    res.val += (*this)[i];
    res.err += (*this)[i] * (*this)[i];
  }
  res.val /= static_cast<long double>(len);
  res.err -= res.val * res.val * static_cast<long double>(len);
  res.err /= static_cast<long double>(len)*static_cast<long double>(len-1);
  res.err = sqrt(res.err);

  return res;
}

estimate datasample::stderr() const {
  const int len = size();
  estimate res = { 0. , 0. };

  if(len==0) return res; //Check for zero len

  for(int i=0; i<len; i++) {
    res.val += (*this)[i];
    res.err += (*this)[i] * (*this)[i];
  }
  res.val /= static_cast<long double>(len);
  res.val = res.err - res.val * res.val * static_cast<long double>(len);
  res.val /= static_cast<long double>(len-1);
  res.val = sqrt(res.val);
  res.err = res.val/sqrt(2.*len);

  return res;
  
}

/* compute the confidence interval corresponding to the low and hi quantile */
void datasample::ci(long double low, long double med, long double hi, long double *il, long double *im, long double *ih) const {
  datasample tmp=(*this);
  int len = size();

  sort( tmp.begin(), tmp.end() );

  int l=int((long double)(len)*(low<0.?0.:low>1.?1.:low));
  int m=int((long double)(len)*(med<0.?0.:med>1.?1.:med));
  int h=int((long double)(len)*(hi<0.?0.:hi>1.?1.:hi));

  *il=tmp[l];
  *im=tmp[m];
  *ih=tmp[h];

}

/* histograms */
void datasample::hist(long double *n, int npt, long double *x, long double w) const {
  datasample tmp=(*this);
  int len = size();
  long double w2=0.5*w;

  sort( tmp.begin(), tmp.end() );
  
  for (int i=0;i<npt; ++i) {
    long double lo=x[i]-w2;
    long double hi=x[i]+w2;
    n[i]=0;
    for (int j=0;j<len; ++j) {
      if (tmp[j]>hi) break;
      if (tmp[j]>lo) n[i]+=1.;
    }
  }
}

void datasample::cumhist(long double *n, int npt, long double *x) const {
  datasample tmp=(*this);
  int len = size();

  sort( tmp.begin(), tmp.end() );
  
  for (int i=0;i<npt; ++i) {
    n[i]=0;
    for (int j=0;j<len; ++j) {
      if (tmp[j]>x[i]) break;
      else n[i]+=1.;
    }
  }

}

estimate datasample::autocorr() const {

  estimate res = { 0. , 0. };
  int len = size();
  if(len==0)  //Empty data set;
    return res;

  long double C0, Ct, rho, tint;
  long double avr=0.0; // Media f(x)
  long double f2=0.0;  // Somma di f(x)f(x+t)
  long double f=0.0;   // Somma di f(x)+f(x+t)

  //t=0
  for (int i=0; i<len; i++) {
    avr+=(*this)[i];
    f2+=(*this)[i]*(*this)[i];
  }
  f = 2.0 * avr;
  avr/=static_cast<long double>(len);
  C0 = (f2/static_cast<long double>(len))-avr*avr;
  tint = 0.5;

  bool valid=false; int M;
  for (M=1; M<len; M++) {
    f2=f=avr=0.0;
    for (int i=0; i<len-M; i++) {
      f2+=(*this)[i]*(*this)[i+M];
      f+=(*this)[i]+(*this)[i+M];
      avr+=(*this)[i];
    }
    avr/=static_cast<long double>(len-M);
    Ct = (f2/static_cast<long double>(len-M))+avr*(avr-(f/static_cast<long double>(len-M)));
    rho = Ct/C0;
    tint += rho;
    //Check for end condition M/4>=tint
    if(M>6.0*tint) {valid=true; break;}
  }

  if(valid) {
    res.val=tint;
    res.err=sqrt(2.0*(2*M+1)/len)*tint;
  }

  return res;

}

estimate datasample::JKavr(const int blsize) const {

  estimate res = { 0. , 0. };
  int sz = size();
  if(sz==0 || blsize==0)  //Empty data set;
    return res;
  int nblocks = sz/blsize; //Numero di blocchi

  if(nblocks<2) {
    cerr<<"datasamples contains too few data!\n";
    return res;
  }

  long double psum[nblocks]; //Somme parziali del blocco
  long double avp[nblocks]; //Medie togliendo un blocco
  long double tsum=0.0; //Somma totale

  for(int j=0; j<nblocks; j++) {
    psum[j]=0.0;
    for(int i=j*blsize; i<blsize*(j+1); i++) {
      psum[j]+=(*this)[i];
    }
    tsum+=psum[j];
  }

  long double avr, err;
  avr = err = 0.;
  for(int j=0; j<nblocks; j++) {
    avp[j]=(tsum-psum[j])/(blsize*(nblocks-1));
    avr += avp[j];
    err += avp[j]*avp[j];
  }

  avr /= static_cast<long double>(nblocks);
  err -= avr * avr * nblocks;
  err *= static_cast<long double>(nblocks-1) / static_cast<long double>(nblocks);
  err = sqrt(err);

  for(int j=blsize*nblocks; j<sz; j++)
    tsum+=(*this)[j];
  avr = tsum/static_cast<long double>(sz);

  res.val = avr;
  res.err = err;
  
  return res;

}

estimate JKcov(const datasample &dt1, const datasample &dt2, const int blsize) {
  estimate res = { 0. , 0. };
  int sz1 = dt1.size(); int sz2 = dt2.size();
  int sz = (sz1>sz2)?sz2:sz1;
  if(sz==0 || blsize==0)  //Empty data set;
    return res;
  int nblocks = sz/blsize;

  //  cerr<<"blocks:"<<nblocks<<"\tblsize:"<<blsize<<endl<<flush;
  if(nblocks<2) {
    cerr<<"datasamples contains too few data!\n";
    return res;
  }

  long double xpsum[nblocks], ypsum[nblocks], zpsum[nblocks]; //Somme parziali del blocco
  long double avp[nblocks]; //Medie togliendo un blocco
  long double xtsum=0.0, ytsum=0.0, ztsum=0.0; //Somme totali
  
  for(int j=0; j<nblocks; j++) {
    xpsum[j]=ypsum[j]=zpsum[j]=0.0;
    for(int i=j*blsize; i<blsize*(j+1); i++) {
      xpsum[j]+=dt1[i]*dt2[i];
      ypsum[j]+=dt1[i];
      zpsum[j]+=dt2[i];
    }
    xtsum+=xpsum[j];
    ytsum+=ypsum[j];
    ztsum+=zpsum[j];
  }

  long double susc, err;
  susc = err = 0.;
  for(int j=0; j<nblocks; j++) {
    avp[j]=((xtsum-xpsum[j])-(ytsum-ypsum[j])*(ztsum-zpsum[j])/(blsize*(nblocks-1)))/(blsize*(nblocks-1));
    susc += avp[j];
    err += avp[j]*avp[j];
  }

  susc /= static_cast<long double>(nblocks);
  err -= susc * susc * nblocks;
  err *= static_cast<long double>(nblocks-1) / static_cast<long double>(nblocks);
  err = sqrt(err);

  for(int j=blsize*nblocks; j<sz; j++) {
    xtsum+=dt1[j]*dt2[j];
    ytsum+=dt1[j];
    ztsum+=dt2[j];
  }
  susc = (xtsum-(ytsum*ztsum)/sz)/sz;

  res.val = susc;
  res.err = err;
  
  return res;


}


estimate datasample::JKvar(const int blsize) const {
  return JKcov(*this, *this, blsize);
}

