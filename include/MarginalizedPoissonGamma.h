#ifndef MARGINALIZEDPOISSONGAMMA_H
#define MARGINALIZEDPOISSONGAMMA_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: MarginalizedPoissonGamma.h
// Description: Implements the Poisson-Gamma model 
//              that is marginalized (that is, 
//              integrated over the nuisance parameters).
// 
// Created: June 11, 2010
// Modifications: 
//
//--------------------------------------------------------------
#include <vector>
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
#include "PDFunction.h"

namespace rfp
{

  /** Implements the Poisson-Gamma model that is marginalized 
      (that is, integrated over the nuisance parameters).
   */
  class MarginalizedPoissonGamma : public rfp::PDFunction
  {
  public:
    
    ///
    MarginalizedPoissonGamma();
     
 	/** Constructor:  
        @param yy - count for background prior
        @param xx - count for signal prior
        @param b - scale for background prior
        @param a - scale for signal prior
        @maxcount - maximum allowed value of observed count 
        (this sets the sizes of some internal arrays).
        <p>
        Here, mean of Poisson is defined by:   
        \f$\sigma*\epsilon + \mu,\f$ where 
        <p>
        \f$\sigma\f$ is the "signal cross section" and 
        parameter of interest,  
        <p> 
        \f$\hat\epsilon \sim (xx + 0.5) / a \f$, and 
        <p> 
        \f$\hat\mu \sim (yy + 0.5) / b \f$.        
    */
    MarginalizedPoissonGamma(std::vector<std::vector<double> >& yy, 
                             std::vector<double>& xx,
						     std::vector<std::vector<double> >& eyy,
						   	 std::vector<double>& exx,
                             int maxcount=50000);
  	/** Constructor:  
        @param yy - count for background prior
        @param xx - count for signal prior
        @param b - scale for background prior
        @param a - scale for signal prior
        @maxcount - maximum allowed value of observed count 
        (this sets the sizes of some internal arrays).
        <p>
        Here, mean of Poisson is defined by:   
        \f$\sigma*\epsilon + \mu,\f$ where 
        <p>
        \f$\sigma\f$ is the "signal cross section" and 
        parameter of interest,  
        <p> 
        \f$\hat\epsilon \sim (xx + 0.5) / a \f$, and 
        <p> 
        \f$\hat\mu \sim (yy + 0.5) / b \f$.        
    */ 
     
    ///
    ~MarginalizedPoissonGamma();
	     
    /** Generate and cache data for one experiment.
        @param sigma - value of parameter of interest
    */
    std::vector<double>& generate(double sigma);
	     
    /** Computes PDF.
        @param xdata - observed data (could be binned)
        @param sigma - value of parameter of interest 
    */

    long double operator() (std::vector<double>& xdata, double sigma);
	     
  private:	     
    std::vector<std::vector<double> > _yy;
    std::vector<double> _xx;
    std::vector<std::vector<double> > _eyy;
    std::vector<double> _exx;
    std::vector<double> _b;
    std::vector<double> _a;
    int _maxcount;
    std::vector<double> _xdata;
    ROOT::Math::Random<ROOT::Math::GSLRngMT> _gslRan;
    int _nbins;
    int _nbkgs;
  };
}

#endif

