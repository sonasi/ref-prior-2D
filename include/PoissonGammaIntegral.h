#ifndef POISSONGAMMA_H
#define POISSONGAMMA_H

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>

typedef std::vector<double>  vdouble;
typedef std::vector<std::vector<double> > vvdouble;

#ifdef __WITH_CINT__
#include "TObject.h" 
#endif

using namespace std;


namespace pg {
  double 
  poissongamma(vector<double>     D,         // Observed counts
               vector<double>	    p,         // Weights "p_j" 
               vector< vector<double> > 	A,         // Counts  "A_ji" for up to 10 sources
               vector< vector<double> > 	f,  // Scale factors associated with counts
			   double mu = 1.0,
               bool returnlog=false, // return log(P) if true
               bool scale=true);     // Scale p_j if true  
}

#endif