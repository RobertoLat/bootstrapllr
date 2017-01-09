#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <list>
#include <iomanip>
#include "datareader.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

struct mypair{
  string dE;
  int n;
};

static bool compare (string first, string second)
{
  if (atof(first.c_str())<atof(second.c_str())) return true;
  else return false;
}

bootstrap * inputdata::populate_bs_sample(int i){
	if(i==1)
		return this->populate_bs_sample();
	
	string E0,E0old;
	string dE,dEold,a;
	stringstream line;
	string sb;
	
	
	map<string,mypair> mymap;
	map<string,mypair>::iterator it;
	
	list<string> mylist;
	list<string>::iterator itl;
	
	
	ifstream fsinput(this->inputfilename.c_str());
	
	while(1) {
		getline(fsinput,sb);
		if(!fsinput.good()) break;
		line.str(sb);
		if(line.good()){
			line>>E0;
			line>>dE;
			
			it=mymap.find(E0);
			if (it==mymap.end() ) {
				mymap[E0].dE=dE;
				mylist.push_back(E0);
			}
			else if( it->second.dE!= dE ){
				cerr << "#[inputdata::populate_bs_sample] Error, in file "<<this->inputfilename<<" there are two E0 with different dE"<<endl;
				exit(1);
			}
		}
	}
	
	mylist.sort(compare);
	
	int Ecounter=0;
	for (itl=mylist.begin(); itl!=mylist.end(); ++itl){
		mymap[*itl].n=Ecounter;
		Ecounter++;
	}
	
	cout << "#[inputdata::populate_bs_sample] found "<<(int) mymap.size()<<" different energy intervals "<< endl;
	
	bootstrap * bs = new bootstrap((int) mymap.size(), this->seed);
	
	fsinput.clear();
	
	fsinput.seekg(0,ios_base::beg);
	long double amean=0.;
	int count=0;
	while(1) {
		getline(fsinput,sb);
		if(!fsinput.good())  break;
		line.str(sb);
		if(line.good()){
			line>>E0;
			line>>dE;
			line >>a;
			Ecounter = mymap[E0].n;
			//cout << Ecounter << " " << E0 << " " << dE << " " <<a << endl;
			count++;
			amean+=atof(a.c_str());
			if(count%this->binsize==0){
			bs->data[Ecounter].E0=atof(E0.c_str());
			bs->data[Ecounter].dE=atof(dE.c_str())-bs->data[Ecounter].E0;
			bs->data[Ecounter].push_back(amean/static_cast<double>(binsize));
			amean=0.;
			}
		}
		line.clear();
	}
	fsinput.close();
	/*if(count%this->binsize!=0){
		cerr<<"#[inputdata::populate_bs_sample] Error: The number of measurements is not a multiple of binsize "<<endl;
		exit(1);
	}*/
		
	cout << "#[inputdata::populate_bs_sample] Each intervall has "<<(int)  bs->data[0].size()<<" measurements in it "<< endl;
	double dEvalue;
	bool samedE=true;
	dEvalue= bs->data[0].dE;
	for(int i=1;i<bs->intervals;i++)
		if(abs(dEvalue-bs->data[i].dE)/dEvalue > 10.E-8) samedE=false;
	if(samedE) 
		cout << "#[inputdata::populate_bs_sample] All intervals have the same dE "<<dEvalue << endl;
	
	
	cout << "#[inputdata::populate_bs_sample] bin size has been set to "<<(int) this->binsize << endl;
	
	// cout <<"# E0      dE     <a_i>       d<a_i>" << endl;
	// for(int i=0;i<bs->intervals;i++)
	//   {
	//     cout<<setprecision(10) <<bs->data[i].E0<<" "<<bs->data[i].dE<<" "<< bs->data[i].avr()<<endl;
	//     //cout<< "#"<<setprecision(10) <<bs->data[i].E0+bs->data[i].dE/2.0<<" "<< bs->data[i].avr()<<endl;
	//   }   
	return bs;
}



bootstrap * inputdata::populate_bs_sample(){
  string E0,E0old;
  string dE,dEold,a;
  stringstream line;
  string sb;


  map<string,mypair> mymap;
  map<string,mypair>::iterator it;
 
  list<string> mylist;
  list<string>::iterator itl;
  

  ifstream fsinput(this->inputfilename.c_str());
  
  while(1) {
    getline(fsinput,sb);
    if(!fsinput.good()) break;
    line.str(sb);
    if(line.good()){
      line>>E0;
      line>>dE;

      it=mymap.find(E0);
      if (it==mymap.end() ) {
	mymap[E0].dE=dE;
	mylist.push_back(E0);
      }
      else if( it->second.dE!= dE ){
	cerr << "#[inputdata::populate_bs_sample] Error, in file "<<this->inputfilename<<" there are two E0 with different dE"<<endl;
	exit(1);
      }
    }
  }
  
  mylist.sort(compare);

  int Ecounter=0;
  for (itl=mylist.begin(); itl!=mylist.end(); ++itl){
    mymap[*itl].n=Ecounter;
    Ecounter++;
  }
  
  cout << "#[inputdata::populate_bs_sample] found "<<(int) mymap.size()<<" different energy intervals "<< endl;

  bootstrap * bs = new bootstrap((int) mymap.size(), this->seed);

  fsinput.clear();

  fsinput.seekg(0,ios_base::beg);

  while(1) {
    getline(fsinput,sb);
    if(!fsinput.good())  break;
    line.str(sb);
    if(line.good()){
      line>>E0;
      line>>dE;
      line >>a;
      Ecounter = mymap[E0].n;
      //cout << Ecounter << " " << E0 << " " << dE << " " <<a << endl;
      bs->data[Ecounter].E0=atof(E0.c_str());
      bs->data[Ecounter].dE=atof(dE.c_str());
      bs->data[Ecounter].push_back(atof(a.c_str()));
    }
    line.clear();
  }
  fsinput.close();

  cout << "#[inputdata::populate_bs_sample] Each intervall has "<<(int)  bs->data[0].size()<<" measurements in it "<< endl;
  double dEvalue;
  bool samedE=true;
  dEvalue= bs->data[0].dE;
  for(int i=1;i<bs->intervals;i++)
    if(abs(dEvalue-bs->data[i].dE)/dEvalue > 10.E-8) samedE=false;
  if(samedE) 
    cout << "#[inputdata::populate_bs_sample] All intervals have the same dE "<<dEvalue << endl;
  

  cout << "#[inputdata::populate_bs_sample] bin size has been set to "<<(int) this->binsize << endl;
 
  // cout <<"# E0      dE     <a_i>       d<a_i>" << endl;
  // for(int i=0;i<bs->intervals;i++)
  //   {
  //     cout<<setprecision(10) <<bs->data[i].E0<<" "<<bs->data[i].dE<<" "<< bs->data[i].avr()<<endl;
  //     //cout<< "#"<<setprecision(10) <<bs->data[i].E0+bs->data[i].dE/2.0<<" "<< bs->data[i].avr()<<endl;
  //   }   
  return bs;
}

void inputdata::read_cmdline(int argc, char* argv[]){
  if(argc!=6 && argc!=7){
    cerr<<"[read_cmdline] Missing Parameter\n\tUsage "<<argv[0]<<" <inputfile> <blocksize> <beta> <lattice size> <seed>  [<n bootstrap>]"<<endl;
    exit(1);
  }

  this->inputfilename=argv[1];
  ifstream input(this->inputfilename.c_str());
  if(!input){
    cerr<<"[read_cmdline] Missing File "<<inputfilename<<endl;
    exit(1);
  }
  

  this->binsize=atoi(argv[2]);
  if(this->binsize<=0) {
    cerr<<"[read_cmdline] Error blocksize must be greater than zero"<<endl;
    exit(1);
  }
  
  this->beta=atof(argv[3]);
  if(this->beta<=0) {
    cerr<<"[read_cmdline] Error beta must be greater than zero"<<endl;
    exit(1);
  }
  
  this->L=atoi(argv[4]);
  if(this->L<=0) {
    cerr<<"[read_cmdline] Error lattice size must be greater than zero"<<endl;
    exit(1);
  }
  
  this->seed=atoi(argv[5]);

  if(argc==7){
    this->nboots=atoi(argv[6]);
    if(this->nboots<=0) {
      cerr<<"[read_cmdline] Error n bootstrap must be greater than zero"<<endl;
      exit(1);
    }
  }
  
}
