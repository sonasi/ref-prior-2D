#ifndef IntrinsicPrior_H
#define IntrinsicPrior_H

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



namespace iprior {
//  vdouble   vNULL2;
//  vvdouble vvNULL2;

  double 
  intrinsicprior(vector<double>	    p,         // Weights "p_j" 
               vector< vector<double> > 	A,         // Counts  "A_ji" for up to 10 sources
               vector< vector<double> > 	f,  // Scale factors associated with counts
               bool returnlog=false, // return log(P) if true
               bool scale=true,
               double conv_int = 1000);     // Scale p_j if true  
}

#endif