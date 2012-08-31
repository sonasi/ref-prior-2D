#ifndef MarginalizedPoissonGammaInt_H
#define MarginalizedPoissonGammaInt_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: MarginalizedPoissonGammaInt.h
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
#include "Math/Interpolator.h"


namespace rfp2d
{

  /** Implements the Poisson-Gamma model that is marginalized 
      (that is, integrated over the nuisance parameters).
   */
  class MarginalizedPoissonGammaInt : public rfp::PDFunction
  {
  public:
    
    ///
    MarginalizedPoissonGammaInt();
     
 	/** Constructor:  
        @param yy - count for background prior
        @param xx - count for signal prior
        @param b - scale for background prior
        @param a - scale for signal prior
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


    MarginalizedPoissonGammaInt(
    						 std::vector<std::vector<double> >& yy,
                             std::vector<std::vector<double> >& xx, 
							 std::vector<std::vector<double> >& eyy,
                             std::vector<std::vector<double> >& exx, 
                             double b, 
                             double a,
							 double maxsigma,
                             int nsigma,
                             int nmass,
                             std::vector<double>& signals,
							 std::vector<std::vector<double> >& firstprior
							  );
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
    ~MarginalizedPoissonGammaInt();
	     
    /** Generate and cache data for one experiment.
        @param sigma - value of parameter of interest
    */
    std::vector<double>& generate(double thismass);
	double randomSample(double thismass);

    /** Computes PDF.
        @param xdata - observed data (could be binned)
        @param sigma - value of parameter of interest 
    */

    long double operator() (std::vector<double>& xdata, double thismass);

  private:

	int _intpoints;
	int _nsigma;
	int _nmass;
	int _nbkgs;
	const int _maxsigma;
	double _priormax;
    std::vector<std::vector<double> > _firstprior;    
    std::vector<std::vector<double> > _xx;
    std::vector<std::vector<double> > _yy;
    std::vector<std::vector<double> > _eyy;
    std::vector<std::vector<double> > _exx;
    std::vector<double> _b;
    std::vector<double> _a;
    std::vector<double> _signals;
    ROOT::Math::Interpolator itp;
    std::vector<double> _xdata;
    ROOT::Math::Random<ROOT::Math::GSLRngMT> _gslRan;
    int _nbins;
  };
}

#endif

