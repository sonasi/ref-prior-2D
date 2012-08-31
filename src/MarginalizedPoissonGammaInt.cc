//--------------------------------------------------------------
//
// Author: Joe Bochenek (bochenek.joe@gmail.com)
// 
// 
// File: MarginalizedPoissonGammaInt.cpp
// Description: Prior probability on Higgs mass from poisson 
// 				gamma model and Jeffrey's prior
// 
// Created: June 04, 2012
// Modifications: 
//
//--------------------------------------------------------------
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <cmath>
#include  <TMath.h>
#include "refpriors/MarginalizedPoissonGammaInt.h"
#include "Math/Interpolator.h"
#include "refpriors/JeffreysPrior.h"
#include "MarginalizedPoissonGamma.h"

using namespace std;
using namespace rfp2d;
using namespace rfp;


MarginalizedPoissonGammaInt::MarginalizedPoissonGammaInt(
												   vector<vector<double> >& yy, 
                                                   vector<vector<double> >& xx,
                     						       vector<vector<double> >& eyy,
					  							   vector<vector<double> >& exx,
                                                   double b, 
                                                   double a,
												   double maxsigma,
  												   int nsigma,
												   int nmass,
					                               vector<double>& signals,
						                           std::vector<std::vector<double> >& firstprior
                                                   )
  : PDFunction(),
    _xdata(vector<double>()),
    _yy(yy),
    _xx(xx),
    _eyy(eyy),
    _exx(exx),
    _maxsigma(maxsigma),
    _nsigma(nsigma),
    _intpoints(nsigma),
    _nmass(nmass),
    _gslRan(ROOT::Math::Random< ROOT::Math::GSLRngMT >()),
    _signals(signals),
    _firstprior(firstprior),
    _b(b),
    _a(a)
{
  _nbins = (int)xx[0].size();
  _nbkgs = yy.size();
  _b.clear();
  _a.clear();
  for(int ibin=0;ibin<_nbins; ++ibin)
	{
      _b.push_back(b);
      _a.push_back(a);
	}

	double max = 0;	
//	vector<ROOT::Math::Interpolator> itp;
//	cout << "inside: " << _nmass << "\t" << _intpoints << endl;
//	cout << "vectors: " << _firstprior.size() << "\t" << _firstprior[0].size() << endl;
	
	for(int imass = 0; imass < _nmass; imass++)
	{
//		vector <double> splinedomain;
//		vector <double> splinevalues;
		for(int isigint = 0; isigint < _intpoints; isigint++)
		{
			double thissigma = (double(isigint) * double(_maxsigma)) / double(_intpoints);
			if(_firstprior[imass][isigint] > max) max = _firstprior[imass][isigint];
//			cout << imass << "\t" << isigint << "\t" << max <<endl;
//			splinevalues.push_back(_firstprior[imass][isigint]);
//			splinedomain.push_back(thissigma);
//			cout << "max: " << max << "\t " << _firstprior[imass][isigint]  << endl;
		}
//		ROOT::Math::Interpolator itp_i( splinedomain, splinevalues, ROOT::Math::Interpolation::kCSPLINE );
//	    itp.push_back(itp_i);
	}
	_priormax = max;
}


MarginalizedPoissonGammaInt::~MarginalizedPoissonGammaInt() 
{}


vector<double>& MarginalizedPoissonGammaInt::generate(double thismass)
{

  if(_nbins == 0)
	{
      Error("MarginalizedPoissonGammaInt", "required input vectors not supplied by user.");
      exit(0);
	}
	double epsilon;
	double mean;

    _xdata.clear();
    
    double sigma =  randomSample(thismass);

//	cout << "sigma: " << sigma << endl;
//	cout << "nbins: " << _nbins << endl;

	for(int ibin=0; ibin<_nbins; ++ibin)
	{
		int massind;
		for(int i = 0; i < _nmass; i++) if(thismass==_signals[i]) massind=i;
		epsilon= _gslRan.Gamma(_exx[massind][ibin]*_xx[massind][ibin], 1.0)/_exx[massind][ibin];
		double bmean = 0;
		for(int ibkg = 0; ibkg<_nbkgs; ibkg++)
		{
		bmean += _gslRan.Gamma(_eyy[ibkg][ibin]*_yy[ibkg][ibin], 1.0)/_eyy[ibkg][ibin];
		}
		mean = epsilon*sigma + bmean;
		double sample = _gslRan.Poisson(mean);
		_xdata.push_back(sample);
	}
	return _xdata;
}



// This samples directly from the prior distribution at a particular mass.
// Uses a rejection sampling method.

double MarginalizedPoissonGammaInt::randomSample(double thismass)
{
	int randx = 1;
	double x;
	double y;
	int massind;
	for(int i = 0; i < _nmass; i++) if(thismass==_signals[i]) massind=i;
	
	do
	{
	x = _gslRan.Uniform(_maxsigma);
	randx = int(x*_nsigma/_maxsigma);
	y = _gslRan.Uniform(_priormax);
	// cout << "priormax: " << _priormax  << "\tthismass: " << thismass << "\t randx: " << randx << "\tx: " << x << "\ty: " << y << "\t prior: " << _firstprior[massind][randx]  << endl;
	}
	while(_firstprior[massind][randx] < y);
	
	return x; 
}


long double MarginalizedPoissonGammaInt::operator() (std::vector<double>& xdata, double thismass)
{
	long double prob = 1.0;	
	// Parameter for numerical integration in Jeffreys prior
	int Mj = 10; // Number of Monte Carlo integration points

	double lo = 0.0;   // lower bound of model density
	double up = 500.0; // upper bound of model density

//	cout << "xdata: " << xdata.size() << endl;

	vector <double> splinevalues;
	// Interpolate over available mass points
	for(int imass = 0; imass < _nmass; imass++)
	{
		MarginalizedPoissonGamma pgamma(_yy, _xx[imass], _eyy, _exx[imass]);
		double priormass = 0; 
		
		// Integrate over sigma
		for(int isigint = 1; isigint < _intpoints+1; isigint++)
		{
			double thissigma = (double(isigint) * double(_maxsigma)) / double(_intpoints);
			vector <double> splinevalues;
			priormass += _firstprior[imass][isigint] * exp(pgamma(xdata, thissigma)) / _intpoints;
//			cout << "mass: " << imass <<  "\t pi(mu|m): " << _firstprior[imass][isigint] << "\t p(D|mu, m): " << pgamma(xdata, thissigma) << "\t priormass: " << priormass << endl;
		}
		cout << "mass: " << imass <<  "\t priormass: " << priormass << "\tlog(prior): " << log(priormass) << endl;

//		priormass *= (double(_maxsigma)/double(_intpoints));
//		splinevalues.push_back(priormass);
		splinevalues.push_back(log(priormass));
	}

    ROOT::Math::Interpolator itp( _signals, splinevalues, ROOT::Math::Interpolation::kCSPLINE);
//	prob = itp.Eval(thismass);
	prob = itp.Deriv2(thismass);
	cout << prob << endl;
	return prob;
}