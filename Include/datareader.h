#ifndef DATAREADER_H
#define DATAREADER_H
#include <string>
#include "bootstrap.h"


class inputdata {
 public:
  std::string inputfilename;
  int nboots;
  int binsize;
  inputdata(){nboots=1000;};
  ~inputdata(){};
  long double beta;
  int L;
  int seed;

  void read_cmdline(int argc, char* argv[]);
  bootstrap* populate_bs_sample();
  bootstrap* populate_bs_sample(int i);

};

#endif //#ifndef DATAREADER_H
