//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: JeffreysPrior.cpp
// Description: Implementation of reference prior using 
//              Jeffreys' rule on Fisher information.  
// 
// Created: June 11, 2010
// Modifications:
//
//--------------------------------------------------------------
#include  <cmath>
#include <iostream>
#include <cstdio>
#include "refpriors/JeffreysPrior.h"


using namespace rfp;
using namespace std;

JeffreysPrior::JeffreysPrior() 
  : PriorFunction(),
    _pdf(0),
    _xmin(0),
    _xmax(1000),
    _m(1000),
    _secondDef(false),
    _debug(false),
    _h(1.e-4),
    _xval(vector<double>()),
    _pval(vector<double>()),
    _interpolator(0)
{}

JeffreysPrior::JeffreysPrior(PDFunction& pdf, 
                             double xmin, 
                             double xmax,
                             int M,
                             bool secondDef,
                             int debug,
                             double h)
  : PriorFunction(),
    _pdf(&pdf),
    _xmin(xmin),
    _xmax(xmax),
    _m(M),
    _secondDef(secondDef),
    _debug(debug),
    _h(h),
    _xval(vector<double>()),
    _pval(vector<double>()),
    _interpolator(0)
{
//  cout << endl;
//  cout << "\t **** Instantiate JeffreysPrior for on-the-fly calculations ****" << endl;
//  cout << endl;
}

JeffreysPrior::JeffreysPrior(PDFunction& pdf, 
                             std::vector<double>& xval,
                             int M,
                             bool secondDef,
                             int debug,
                             double h)
  : PriorFunction(),
    _pdf(&pdf),
    _xmin(xval.front()),
    _xmax(xval.back()),
    _m(M),
    _secondDef(secondDef),
    _debug(debug),
    _h(h),
    _xval(xval),
    _pval(vector<double>(xval.size(), 0)),
    _interpolator(0)
{
//  cout << endl;
//  cout << "\t **** Instantiate JeffreysPrior for calculation based on given interpolation points ****" << endl;
//  cout << endl;

  // Move slightly away from the user-defined boundaries

  int N = xval.size();
  double minStep = min(_xval[1] - _xval[0], 
                       _xval[N-1] - _xval[N-2]);
  if ( 4*_h > minStep ) _h = minStep / 10; 
  
  _xval[0] = _xmin + 4*_h;
  _xval[N-1] =  _xmax - 4*_h;

  cout << endl;
  cout << "\t **** NOTE: Moving slightly away from user-defined boundaries:" 
       << endl;
  cout << "\t **** Bounds of probability density function are set to:" 
       << endl;
  cout << "\t **** xmin = " << _xmin 
       << ", xmax = " << _xmax << "," << endl;
  cout << "\t **** and the first and last interpolation points are shifted to:"
       << endl;
  cout << "\t **** " << _xval[0] << " and " << _xval[N-1] 
       << ", respectively." << endl;
  cout << endl;

  // Compute reference prior over requested range

  for(int i=0; i < N; ++i)
    {
      _pval[i] = (*this)(_xval[i]);
      if ( _debug > 0)
        {
          char record[255];
          sprintf(record, "xval[%3.3d]\t= %12.5f\tpval[%3.3d]\t= %12.5f", i, _xval[i], i, _pval[i]);
          cout << record << endl;
        }
    }

  // Set up interpolator

  _interpolator = new ROOT::Math::
    Interpolator(_xval, _pval,
                 ROOT::Math::Interpolation::kLINEAR);
}

JeffreysPrior::~JeffreysPrior()
{
  if ( _interpolator != 0 ) delete _interpolator;
}

//-------------------------------------------------------------
// Compute reference prior at specified value of parameter 
// of interest
//-------------------------------------------------------------
double JeffreysPrior::operator() (double xVal)
{

  if ( _interpolator != 0 )
	{
      double prior = _interpolator->Eval(xVal);
      return  prior;
	}
        
  //-------------------------------------
  // first, check for bounds of function
  //-------------------------------------
  
  if(xVal-2*_h < _xmin)
    {
      cout << "WARNING: xVal (= " << xVal 
           << ") - 2*" << _h 
           << " is less than lower bound of prob. density function (= " 
           << _xmin << ");" << endl;
      cout << "\t hence, shifting xVal to " << _xmin + 4*_h << endl;

      xVal =  _xmin + 4*_h;
    }
  if(xVal+2*_h > _xmax)
    {
      cout << "WARNING: xVal (= " << xVal 
           << ") + 2*" << _h 
           << " is greater than upper bound of prob. density function (= "
           << _xmax << ");" << endl;
      cout << "\t hence, shifting xVal to " << _xmax - 4*_h << endl;

      xVal =  _xmax - 4*_h;   
    }
     
  PDFunction& pdf = *_pdf; // Create a reference to the PDFfunction

  double prior;

  if(!_secondDef)
    {
    
      //----------------------------------------
      // loop over _M for numerical integration
      //----------------------------------------
      long double sumP = 0.0; 
      long double sumN = 0.0; 
      double t1 = 4.0/3;
      double t2 = 1.0/12;
      for (int iM=0; iM < _m; ++iM)
        {
          //----------------------------------------
          // generate fake-data
          //---------------------------------------- 
		  cout << "numint: " << iM << "\txVal: " << xVal << "\t" << endl;

          vector<double>& xdata = pdf.generate(xVal); 
          // Compute 2nd derivative of ln(P) to O(h^4)
		
          double fm1 = pdf(xdata, xVal - _h) ;
          double fm2 = pdf(xdata, xVal - 2*_h) ;
          double f   = pdf(xdata, xVal) ;
          double fp1 = pdf(xdata, xVal + _h) ;
          double fp2 = pdf(xdata, xVal + 2*_h) ;

//		  exit(0);
          double term1 = t1 * ( fp1 - 2*f + fm1 );
          double term2 =-t2 * ( fp2 - 2*f + fm2 );
      
          // Sum all positive and all negative terms separately
          if ( term1 < 0 )
            sumN +=-term1;
          else
            sumP += term1;

          if ( term2 < 0 )
            sumN +=-term2;
          else
            sumP += term2;
        } // iM    

      prior = sqrt( abs(sumN - sumP)/(_h * _h * _m) );
	  cout << prior << endl;
    }
  else 
    {
      // fake data sets (first dim: entries; second dim: bins)
      vector<vector<double> > xdata;
      //----------------------------------------
      // loop over _M for numerical integration
      //----------------------------------------
      double sum = 0.0; 

      for (int iM=0; iM < _m; ++iM)
        {
          //----------------------------------------
          // generate fake-data
          //---------------------------------------- 

          vector<double>& xdata = pdf.generate(xVal);
          // Get 2nd Derivative directly from pdf class (which uses spline interpolation or something)
		  sum += pdf(xdata, xVal);
    	}
    
      prior = sqrt(-sum / _m);
    }

  return prior;

}
