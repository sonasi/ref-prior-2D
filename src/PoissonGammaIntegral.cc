///////////////////////////////////////////////////////////////////////////////
// File:	poissongamma.cc
// Description:	x = poissongamma(D,p,A,f)
// 
//              compute the marginal likelihood:
//              a product of Poisson distributions
//              and a prior that is a product of gamma distributions.
//
//              Model:
//
//              d_i = Sum p_j * a_ji,  i=1..M,   j=1..N
// 
//              D_i  is count corresponding to true mean d_i
//              A_ji is count corresponding to true mean a_ji
//              f_ji is an optional scale factor associated with A_ji
//
//              D is a vector of M observed counts (that is, over M bins)
//              p is a vector of N parameters (that is, over N sources)
//              A is a vector of vectors of size N x M counts.
//
//              Simple 1-D Application:
//              
//              d =  xsec * e_a * a + e_b * b 
//
//              where p_1 = xsec * e_a
//                    p_2 = e_b
//                    a_11= a
//                    a_21= b
//
//              and e_a, e_b are scale factors defined by
//              e = sigma^2 / estimate
//
// WARNING: The number of terms grows exponentially with the number of 
//          sources and the maximum bin count. Therefore, this routine
//          really only makes sense for rather few sources with modest
//          bin counts. 
//
// Created:     20-Dec-2001 Harrison B. Prosper
//              11-Mar-2004 HBP add additional interfaces
//              07-Jun-2005 HBP increase number of sources to 10! See warning
//                          above.
//              19-Feb-2011 HBP generalize to make use of user-supplied
//                          scale factors. This is needed for weighted 
//                          histograms. Reduce maximum number of sources to 6
//              27-Apr-2011 HBP increase number of sources to 8
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include "PoissonGammaIntegral.h"
#include <gsl/gsl_sf_gamma.h>

using namespace std;

// Global variables to avoid memory allocation.
namespace pg {
  const int MAXSRCS=8; // Maximum number of sources

  const int MAXSRC = MAXSRCS;     // Maximum number of sources
  const int MAXBUF = 50000; // Maximum count per bin

  double 
    poissongamma(vector<double>	D,   // Counts  "D_i" for data.
                 vector<double>	p,   // Weights "p_j" 
                 vector< vector<double> > 	A,   // Counts  "A_ji" for up to 8 sources
                 vector< vector<double> > 	f,   // scale factor for  "A_ji"
				 double mu,
                 bool returnlog, // return log(P) if true
                 bool scale)     // Scale p_j if true  
  {
	long double c[MAXSRC][MAXBUF];
	double ns[MAXSRC];
	double s [MAXSRC];
	double x [MAXSRC];
	double y [MAXSRC];

    int N = p.size(); // Number of sources (N)
    int M = D.size(); // Number of bins    (M)


    // Check inputs
    if ( A.size() != (unsigned int)N )
      {
        std::cout << "**Error - poissongamma - mis-match in number of sources"
                  << endl
                  << "size(p): " << N << " differs from size(A) = " << A.size()
                  << std::endl;
        exit(0);
      }

    if ( A[0].size() != (unsigned int)M )
      {
        std::cout << "**Error - poissongamma - mis-match in number of binss\n"
                  << "size(D): " << M << " differs from size(A[0]) = " 
                  << A[0].size()
                  << std::endl;
        exit(0);
      }

    if ( M < 1 ) return -1.0;


    if ( ( N < 1 ) || ( N > MAXSRC ) ) return -2.0;
  
    // Get total counts per source



    for (int j = 0; j < N; ++j)
      {
        ns[j] = 0;
        for (int i = 0; i < M; ++i)
          ns[j] += A[j][i];
      }
      

    // loop over the M terms of the product,
    // corresponding to the M bins

    double prob;
    if ( returnlog )
      prob = 0.0;
    else
      prob = 1.0;

    for (int i = 0; i < M; ++i)
      {
        int Di = (int)D[i]; // data count for bin i

//		cout << "Bin " << i << endl;
//		std::cout << "N: " << N << ", M: " << M << "\t";
//		cout << "D_" << i << " \t" << D[i] << "\t";
        // compute terms of sum from zero to D

        // first do zero...      
        double total = 0;
        for (int j = 0; j < N; ++j)
          {
//          	cout << "\tA_" << j << "," <<  i << " = " << p[j]*A[j][i]/ns[j]  << " +/- " << sqrt(A[j][i]/f[j][i]) << endl;
//          	cout << "f_" << j << "," <<  i << " = " << f[j][i] << endl;
//         	cout << "p_" << j << "," <<  i << " = " << p[j] << endl;

			total += p[j]*A[j][i]/ns[j];

            // Normalize sources to unit area so that
            // x[i] becomes the actual source count.
            if ( scale )
        		x[j] = p[j];
        	else
              x[j] = p[j];
            s[j] = A[j][i];

			// scale the signal by a scale factor
			s[0] *= mu;
			f[0][i] /= mu;
			

//			f[j][i] = 1.;
            // Apply user supplied scale factor
            // This is needed to take account of weighted histograms
            if ( (int)f[j].size() == M )
              {
                c[j][0] = 0;
                if ( f[j][i] > 0 )
                  {
                    x[j] *= f[j][i];
                    s[j] *= f[j][i];
                  }
                else 
                  {
                    x[j] = 100000.;
                    s[j] = 0.;                    
                  }
				  
              }
	    	  
//	    	  c[j][0] = exp( (s[j] + 0.5) * log(x[j]) - (s[j] - 0.5) * log(1 + x[j]) );
//	    	  c[j][0] = pow( (1 + x[j]), -(s[j]-0.5) );
//			  c[j][0] = pow(x[j], (s[j]+0.5) ) / pow( (1 + x[j]), (s[j]+0.5) );
			  c[j][0] = exp(  (s[j]-0.5) * log(x[j])  - (s[j]-0.5) * log(x[j]+1) );

//			  cout << "\t\t " << "f: " <<  f[j][i] << "\t x: " <<  x[j] <<  "\t s: " <<  s[j] << endl;	
//			  cout << "\t\t c[" << j << "][0] = pow(" << x[j] << ", (" << s[j] << " - 0.5) ) * pow( (1 + " << x[j] << "), -( " << s[j] << " - 0.5) ) = " << c[j][0] << endl;
	

//			  cout << c[j][0] << "\t" ;
          }


//        cout << "Total\t" << total << endl;
//        cout << "Diff\t" << fabs(D[i] - total) << endl;
        // ...then 1 to D
        if ( Di > 0 )
          {
//            cout << "Di: " << Di << endl;
            for (int k = 1; k < Di+1; ++k)
              for (int j = 0; j < N; ++j)
               {
                  c[j][k] = 
//                  pow(x[j], (s[j]+0.5) ) * pow( (1 + x[j]), -(k + s[j]+0.5) ) *
               	  exp(  (s[j]-0.5) * log(x[j])  - (k + s[j]-0.5) * log(x[j]+1) ) *
                  gsl_sf_gamma(s[j] + k + 0.5) / (  gsl_sf_gamma(k+1) * gsl_sf_gamma( s[j] + 0.5 )  ) ;
 //                	pow(x[j], k ) * pow( (1 + x[j]), -(k + s[j] + 1) ) *
 //                 	gsl_sf_gamma(s[j] + k + 1.0) / (  gsl_sf_gamma( k + 1 ) * gsl_sf_gamma( s[j] + 1.0 )  ) ;
//                   if ( s[j] > 100 ) cout << "s[ " << j << "]: " <<  s[j] + k + 0.5  <<  "\t" <<  k+1  << "\t" <<  s[j] + 0.5   << endl;

  //                cout << "\t\t s[" << j << "], " << k << ": " << s[j] << endl;

                  if(x[j] == 0) c[j][k] = 0;
				}
          }

        // compute sum

        double sum = 0.0;
        switch (N)
          {
          case 1:
                sum += c[0][Di];
//              sum += gsl_sf_gamma(s[0] + Di) / (gsl_sf_gamma(Di+1)*gsl_sf_gamma(s[0])) * pow(1 - x[0], (Di)) * pow(x[0], s[0]);
            break;
          case 2:
            for (int j = 0; j < Di+1; ++j){
              sum += 
                c[0][j] *
                c[1][Di-j];
//				cout << "\t Di: " << Di << "\tk: " << j << " prob: " << sum << endl;
			}
//              sum += gsl_sf_gamma(s[0] + Di) / (gsl_sf_gamma(Di+1)*gsl_sf_gamma(s[0])) * pow(1 - x[0], (Di)) * pow(x[0], s[0]);

            break;

          case 3:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                sum += 
                  c[0][j] *
                  c[1][k] *
                  c[2][Di-j-k];
            break;

          case 4:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][Di-j-k-l];
            break;
          
          case 5:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    sum += 
                      c[0][j] *
                      c[1][k] *
                      c[2][l] *
                      c[3][m] *
                      c[4][Di-j-k-l-m];
            break;
          
          case 6:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    for (int n = 0; n < Di+1-j-k-l-m; ++n)
                      sum += 
                        c[0][j] *
                        c[1][k] *
                        c[2][l] *
                        c[3][m] *
                        c[4][n] *
                        c[5][Di-j-k-l-m-n];
            break;
          
          case 7:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    for (int n = 0; n < Di+1-j-k-l-m; ++n)
                      for (int jj = 0; jj < Di+1-j-k-l-m-n; ++jj)
                        sum += 
                          c[0][j] *
                          c[1][k] *
                          c[2][l] *
                          c[3][m] *
                          c[4][n] *
                          c[5][jj] *
                          c[6][Di-j-k-l-m-n-jj];
            break;

          case 8:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    for (int n = 0; n < Di+1-j-k-l-m; ++n)
                      for (int jj = 0; jj < Di+1-j-k-l-m-n; ++jj)
                        for (int kk = 0; kk < Di+1-j-k-l-m-n-jj; ++kk)
                          sum += 
                            c[0][j] *
                            c[1][k] *
                            c[2][l] *
                            c[3][m] *
                            c[4][n] *
                            c[5][jj] *
                            c[6][kk] *
                            c[7][Di-j-k-l-m-n-jj-kk];
            break;
          };
//		cout << "\tD[" << i << "] = " << Di << "\t" << "A[" << 0 << "] = " << s[0] << "\t" << "p[" << 0 << "] = " << x[0] << "\tf[" << i << "] = " << f[0][i] << "\t Sum:" << log(sum) << endl;


//		Purely Poisson based likelihood		
//		double binsum = 0;
//		for (int j = 0; j < N; ++j)
//		{binsum += s[j];}
//		double fact;
//		if(Di > 0) { fact = gsl_sf_gamma(Di+1);} else {fact = 1.;}
//		sum = pow(binsum, Di) * exp(-binsum) / fact;


        if ( returnlog )
          prob += log(sum);
        else
          prob *= sum;

//        cout << "\tp final prob:\t" <<  prob << endl;

//          cout << "sum: " << sum << "\t" << endl;
//          cout << "p_" << i << ":\t" <<  prob << endl;

//		cout << endl;
      }

//	cout << "Probability: " << prob << endl;
    return prob;
  }
}

 
