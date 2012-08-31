#ifndef JEFFREYSPRIOR_H
#define JEFFREYSPRIOR_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: JeffreysPrior.h
// Description: Implementation of reference prior using 
//              Jeffreys' rule on Fisher information.  
// 
// Created: June 11, 2010
// Modifications:
//
//--------------------------------------------------------------
#include <vector>
#include "PDFunction.h"
#include "PriorFunction.h"
#include "Math/Interpolator.h"

namespace rfp
{
  /** Conpute Jeffreys prior using two different methods:
      <p> 1. By setting <i>secondDef</i> = false,  
      \f$ 
      \pi(\theta) = \sqrt{
      E\left[ -\frac{d^2 \ln p(x|\theta)}{d\theta^2} \right]}
      \f$
      <p> 2. By setting <i>secondDef</i> = true,
      \f$ 
      \pi(\theta) = \sqrt{
      E\left[ \left( \frac{d \ln p(x|\theta)}{d\theta} \right)^2\right] }
      \f$
   */
  class JeffreysPrior : public rfp::PriorFunction
  {
  public:
	
    ///
	JeffreysPrior();

	/** Constructor for on-the-fly calculations.
        @param pdf - probability model
        @param xmin - lower bound of parameter space
        @param xmax - upper bound of parameter space
        @param M - number of MC integration points
        @param secondDef - if true, compute prior using formula 2 above
        @param debug - 0: least printout; 1: first level of printout
        @param h - step size for \f$O(h^4)\f$ 2nd derivative calculation
     */
    JeffreysPrior(rfp::PDFunction& pdf, 
                  double xmin=0.0, 
                  double xmax=1000.0,
                  int M=1000,
                  bool secondDef=false,
                  int debug=0,
                  double h=1.e-4		
                  );
	   
	/** Constructor for calculations using cached interpolation points.
        @param pdf - probability model
        @param xval - interpolation points spanning parameter space
        @param M - number of MC integration points
        @param secondDef - compute using alternate formula (Eq., arXiv:)
        @param debug - 0: least printout; 1: first level of printout
        @param h - step size for \f$O(h^4)\f$ 2nd derivative calculation
     */
    JeffreysPrior(rfp::PDFunction& pdf, 
                  std::vector<double>& xval,
                  int M=1000,
                  bool secondDef=false, 
                  int debug=0,
                  double h=1.e-4
                  );
    /// 
    ~JeffreysPrior();
	    
    /** Compute prior: 
        either on-the-fly or using cached interpolation points 
        depending on which
        constructor was called.
    */
    double operator() (double xVal);
    double operator() (double xVal, double yVal);
	    
     /** Return interpolation points if corresponding constructor 
         called before.
    */
    std::vector<double>& getX() {return _xval;}

    /** Return prior computed at interpolation points if 
        corresponding constructor called before.
    */
    std::vector<double>& getP() {return _pval;}

  private:
	    
    rfp::PDFunction* _pdf;
    double _xmin;
    double _xmax;
    int _m;
    bool _secondDef;
    int _debug;
    double _h;	    
    std::vector<double> _xval;
    std::vector<double> _pval;
    ROOT::Math::Interpolator* _interpolator;
      
  };
}

#endif


