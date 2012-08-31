//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: MarginalizedPoissonGamma.cpp
// Description: Implements the Poisson-Gamma model 
//              that is marginalized (that is, 
//              integrated over the nuisance parameters).
// 
// Created: June 11, 2010
// Modifications: June 04, 2012 Joe Bochenek - modified to include gamma functions 
// 				  instead of factorials
//
//--------------------------------------------------------------
#include <iostream>
#include <cmath>
#include  <TMath.h>
#include "MarginalizedPoissonGamma.h"
#include "PoissonGammaIntegral.h"
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace rfp;


MarginalizedPoissonGamma::MarginalizedPoissonGamma(
												   vector<vector<double> >& yy, 
                                                   vector<double>& xx,
                     						       vector<vector<double> >& eyy,
					  							   vector<double>& exx,
                                                   int maxcount)
  : PDFunction(),
    _yy(yy),
    _xx(xx),
    _eyy(eyy),
    _exx(exx),
    _gslRan(ROOT::Math::Random< ROOT::Math::GSLRngMT >()),
    _xdata(vector<double>()),
    _maxcount(maxcount)
{
  _nbkgs = yy.size();
  _nbins = (int)xx.size();
  if((int)_yy[0].size() != _nbins)
	{
      Error("MarginalizedPoissonGamma", "input vectors have different sizes.");
      exit(0);
	}
}



MarginalizedPoissonGamma::~MarginalizedPoissonGamma() 
{}

vector<double>&  
MarginalizedPoissonGamma::generate(double sigma)
{
	
	if(_nbins == 0)
	{
	  Error("MarginalizedPoissonGamma", "required input vectors not supplied by user.");
	  exit(0);
	}
	double epsilon;
	double mean;
	
	_xdata.clear();
	
//	cout << "sigma: " << sigma << endl;

//    cout << "generate" << endl;
	
	for(int ibin=0; ibin<_nbins; ++ibin)
	{
	  
	  epsilon= _gslRan.Gamma(_exx[ibin]*_xx[ibin], 1.0)/_exx[ibin];

	  double bmean = 0;
	  for(int ibkg = 0; ibkg<_nbkgs; ibkg++)
	  {
		bmean += _gslRan.Gamma(_eyy[ibkg][ibin]*_yy[ibkg][ibin], 1.0)/_eyy[ibkg][ibin];
//	    cout << "\tb[" << ibin << "]=" << _yy[ibkg][ibin] << "\tbmean = " << bmean << endl;
	  }

//		cout << "\ts " << _xx[ibin] << "\tepsilon: " << epsilon << endl;
	  mean = epsilon*sigma + bmean;
	  double sample = _gslRan.Poisson(mean);
//      if(sample > 0) cout << "sample: " << epsilon << endl;
	  _xdata.push_back(sample);
	
//	  cout<<"_xdata["<<ibin<<"] = "<<_xdata[ibin]<<endl;
	
	}
	
	return _xdata;
    
}

long double MarginalizedPoissonGamma::operator() (std::vector<double>& xdata, double sigma)
{
	long double C1[_maxcount+1];
	long double C2[_maxcount+1];
	double A1;  // signal count
	double A2;  // background count
	double nn;  // observed count
	double p1;
	double p2;  
	
	// loop over bins
//	cout << "operator" << endl;
	int nbins = (int)xdata.size();
    
	// Scale factor
	vector< vector<double> > A;
	vector< vector<double> > f;
	vector<double> p;
	vector<double> mva_src_s;
	vector<double> mva_error_s;
	// SIGNAL
	for(int i; i < _nbins; i++)
	{
		double bincontent = sigma * _xx[i];
		double binerror = sigma * _exx[i];
		double error = 0. ;
		if(bincontent < 0) bincontent = 0;
		//				if(bincontent > 0) error = double( binerror*binerror/bincontent );		
		if(bincontent > 0) error = 50 / bincontent;
		mva_src_s.push_back(bincontent);
		mva_error_s.push_back(error);
	}

	A.push_back(mva_src_s);
	f.push_back(mva_error_s);
	p.push_back( 1.0 );

	// BACKGROUND
	for(int j; j < _nbkgs; j++)
	{
	A.push_back(_yy[j]);
	f.push_back(_eyy[j]);	
	p.push_back( 1.0 );
	}

	double prob = pg::poissongamma(xdata, p, A, f, 1.0, true, false);

//	cout << "prob(" << sigma << "): " << prob << endl;
	return prob;
}