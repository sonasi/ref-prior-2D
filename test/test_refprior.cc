//--------------------------------------------------------------
//
// Author: Joe Bochenek (bochenek.joe@gmail.com)
// 
// 
// File: testPGammaMethod2.cpp
// Description: Test ReferencePrior for a PoissonGamma distribution
//                     that is marginalized (integrated over 
//                                           nuisance parameters).
// 
//              The mean of the Poisson distribution is:
//              mean = sigma*epsilon + mu; 
//                     Here, sigma is the parameter of interest, 
//                     sometimes called "x" below, and epsilon 
//                     and mu pertain to the signal and background, 
//                     and are defined by Gamma priors with 
//                     parameters (xx, a) and (yy, b) respectively.  
//                     (See note in doc/doc.pdf)      
//
//
// This program creates two files:
//  1. PGamma.eps, which compares the numerical Reference prior 
//     and the numerical Jeffreys' prior calculations with
//     an analytical calculation possible for a single bin.
//  2. PGammaJeffreysPrior.eps, which compares the numerical 
//     Jeffreys' prior calculation through two different approaches.
//
//
// Created: June 06, 2012
// Modifications: 
//
//--------------------------------------------------------------

#include "refpriors/JeffreysPrior.h"
#include "refpriors/MarginalizedPoissonGammaInt.h"
#include "MarginalizedPoissonGamma.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TMath.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TROOT.h>

#include "samples.h"
#include "Slurper.h"


//-------------------------------------------------------------
using namespace std;
using namespace rfp;
using namespace rfp2d;

//-------------------------------------------------------------
int main()
{


	//---------------------------------------------------------------------------
	// Make Vectors/Histograms
	//---------------------------------------------------------------------------

	const int xsec_plots = 60;
	double lumi_scale = 1.0;
	//---------------------------------------------------------------------------
	// RETRIEVE HISTOGRAMS
	//---------------------------------------------------------------------------
	TH2F* histos_sig[maxsig];

	TFile file("output/hist_ensemble_noiso.root");
	file.cd();

	// ZZ background
	TH2F *histos_bkg=(TH2F*) gDirectory->Get("bkg2d");
	histos_bkg->SetDirectory(0);
	histos_bkg->Sumw2();
	histos_bkg->GetXaxis()->SetTitle("BNN(x)");
	histos_bkg->GetYaxis()->SetTitle("m_{4l}");
	double signal_weights[maxsig] = {0};

	// Data driven reducible background
	TH2F *histos_zjets=(TH2F*) gDirectory->Get("zjets2d");
	histos_zjets->SetDirectory(0);
	histos_zjets->Sumw2();
	histos_zjets->GetXaxis()->SetTitle("BNN(x)");
	histos_zjets->GetYaxis()->SetTitle("m_{4l}");

	// Data
	TH2F *histos_dat=(TH2F*) gDirectory->Get("dat2d");
	histos_dat->SetDirectory(0);	
	histos_dat->Sumw2();
	histos_dat->GetXaxis()->SetTitle("BNN(x)");
	histos_dat->GetYaxis()->SetTitle("m_{4l}");

	// Print some histograms for debugging

	TFile f5("output/hist_ensemble_morph.root");
	f5.cd();

	//---------------------------------------------------------------------------
	// FILL SIGNAL VECTOR
	//---------------------------------------------------------------------------

	for(int k=0; k < maxsig; k++) 
	{
		sigNames[k] = ymin + k;
	}


	for(Int_t k=0; k < maxsig; k++) 
	{
		std::stringstream _name;
//		sigNames[k] = double(2*k + 116);
		_name << "histos_sig_" << sigNames[k]; 
		TString name = _name.str();
		cout << k << "\t" << sigNames[k] << "\t" << name << endl;
		histos_sig[k]= (TH2F*)gDirectory->Get(name);
		histos_sig[k]->Scale(lumi_scale);
		signal_weights[k] = histos_sig[k]->Integral();
		cout << sigNames[k] << ": " << signal_weights[k] << endl;
		cout << "signal[" << k << "]: " << sigNames[k] << "\t" << signal_weights[k] << endl;
		histos_sig[k]->Sumw2();
	}



	//---------------------------------------------------------------------------
	// FILL DATA VECTOR
	//---------------------------------------------------------------------------
	
	vector<double> mva_src_d;
	vector<double> mva_error_d;			

	// BACKGROUND
	for(int i=1; i < xbins+1; i++)
	{
	for(int j=1; j < ybins+1; j++)
	{
		double bincontent = histos_dat->GetBinContent(i, j);
		mva_src_d.push_back(bincontent);
		//				cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " exp error: " << exp_error << " exp error: " << tot_error  << endl;
	}
	}


	//---------------------------------------------------------------------------
	// FILL BACKGROUND VECTORS
	//---------------------------------------------------------------------------
	vector< vector<double> > B;
	vector< vector<double> > fB;

	histos_bkg->Scale(lumi_scale);			
	// Calculate errors in bins for the backgrounds

	cout << xbins << " : " << ybins << endl;

	vector<double> mva_src_b;
	vector<double> mva_error_b;			
	// BACKGROUND
	// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
	for(int i=1; i < xbins+1; i++)
	{
	for(int j=1; j < ybins+1; j++)
	{
		double bincontent = histos_bkg->GetBinContent(i, j);
		double binerror = histos_bkg->GetBinError(i, j);
		double error = 0;
//		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_b.push_back(bincontent);
		mva_error_b.push_back(error);
//		cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " eff error: " << bincontent*error << " rel error: " << binerror/bincontent << endl;
		cout << "bin " << i << " bin " << j << " count " << mva_src_b.size() << endl;
		//				cout << hmlp[j]->GetBinContent(i) << " +/-" << hmlp[j]->GetBinContent(i)/sqrt(hmlp[j]->GetBinError(i))  << endl;
	}
	}
	
	B.push_back(mva_src_b);
	fB.push_back(mva_error_b);
	
	cout << "mva_src_b: " << mva_src_b.size() << endl;

	histos_zjets->Scale(0.0484985);
	histos_zjets->Scale(lumi_scale);
	histos_zjets->Smooth();
	histos_zjets->Smooth();
	histos_zjets->Smooth();
	histos_zjets->Smooth();
			
	// Calculate errors in bins for the backgrounds
	vector<double> mva_src_zjets;
	vector<double> mva_error_zjets;			
	// BACKGROUND
	for(int i=1; i < xbins+1; i++)
	{
	for(int j=1; j < ybins+1; j++)
	{
		double bincontent = histos_zjets->GetBinContent(i, j);
		double binerror = histos_zjets->GetBinError(i, j);
		double error = 0;
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_zjets.push_back(bincontent);
		mva_error_zjets.push_back(error);
	}
	}

	B.push_back(mva_src_zjets);
	fB.push_back(mva_error_zjets);
	
	
	cout << "Z+jets Integral: " << histos_zjets->Integral() << endl;
	cout << "ZZ Integral: " << histos_bkg->Integral() << endl;
	cout << "Signal: " << histos_sig[12]->Integral() << endl;
	cout << "Data Integral: " << histos_dat->Integral() << endl;



	//---------------------------------------------------------------------------
	// Remake  A, f vectors for IntrinsicPrior algorithm
	//---------------------------------------------------------------------------

	for(Int_t k=0; k < maxsig; k++) 
	{
		std::stringstream _name;
		_name << "histos_sig_" << sigNames[k]; 
		TString name = _name.str();
		histos_sig[k]= (TH2F*)gDirectory->Get(name);
		histos_sig[k]->Scale(lumi_scale);
		signal_weights[k] = histos_sig[k]->Integral();
		cout << sigNames[k] << ": " << signal_weights[k] << endl;
	}
	
	vector<double> signals;
	// Loop through m_H and \mu to calculate 2D prior and store the grid in a TGraph2D
	double xsec_factors[xsec_plots]; 
	double count = 0;


    TGraph2D *lh_tgraph = new TGraph2D();
    TGraph2D *lh_tgraph_norm = new TGraph2D();

	clock_t start, end;
	start = clock();





	//---------------------------------------------------------------------------
	// Remake  A, f vectors for PoissonGammaFit2
	//---------------------------------------------------------------------------
	vector< vector<double> > A2;
	vector< vector<double> > f2;
	vector<double> p2;

	// SIGNAL
	for(int l=0; l<maxsig; l++) {
		// Scale factor
		vector<double> mva_src_s;
		vector<double> mva_error_s;
		signals.push_back(sigNames[l]);		
		// SIGNAL
		for(int i=1; i < xbins+1; i++)
		{
		for(int j=1; j < ybins+1; j++)
		{

				double bincontent = histos_sig[l]->GetBinContent(i,j);
				double binerror = histos_sig[l]->GetBinError(i,j);				
				double error = 0;
				if(bincontent > 0) error = double(binerror*binerror/bincontent);		
				mva_src_s.push_back(bincontent);
				mva_error_s.push_back(error);
		}
		}

		A2.push_back(mva_src_s);
		f2.push_back(mva_error_s);
		p2.push_back(1.0);
	}

	cout << "A2: " << A2.size() << endl;
	cout << "Hellow orlds" << endl;

	int Mj = 10; // Number of Monte Carlo integration points
	double lo = 0.0;   // lower bound of model density
	double up = 500.0; // upper bound of model density
	double xref = 1.00;

    TGraph2D *dt1 = new TGraph2D();
	
    TGraph2D *dt = new TGraph2D();
    
//	double prior[maxsig][xsec_plots];
	double function[maxsig][xsec_plots];

//	MarginalizedPoissonGamma pgamma(mva_src_b, A2[l], 1.0, 1.0);
 //   vector<double>& xdata = pgamma.generate(2.0); 
    vector<vector<double> > prior;

	for(int l=0; l<maxsig-1; l++) {
		MarginalizedPoissonGamma pgamma(B, A2[l], fB, f2[l]);
		JeffreysPrior jeffprior(pgamma, lo, up, Mj);
		vector<double> thisrow;
		cout << "Mass: " << l << endl;

		for( int j = 0; j < xsec_plots; j++)
		{
			double xsecfactor =  j*0.1;
			xsec_factors[j] = xsecfactor;
			double point = jeffprior(xsecfactor);
//			double point = 1.;
			cout << "mass: " << l <<  "\txsec[" << j << "]\t" << point << endl; 
//			splinevalues.push_back(point);
//			prior[l][j] = point;
			thisrow.push_back(point);
//			cout << point << "\t";
		}
		prior.push_back(thisrow);
	}
	
	
	for(Int_t m=0; m < maxsig-1; m++)for(Int_t j=0; j < xsec_plots; j++)  dt1->SetPoint( (j+(m)*xsec_plots), sigNames[m], xsec_factors[j], prior[m][j]);


	TCanvas *A4 = new TCanvas("A3","Plotting Canvas",150,10,990,660);
	dt1->Draw("colz");
	A4->SaveAs("output/muprior.gif");
	
	TFile outfile("output/refprior.root", "RECREATE");
	outfile.cd();
	dt1->SetName("prior");
	dt1->Write();
	outfile.Close();

	exit(0);

	MarginalizedPoissonGammaInt pgammaint(B, A2, fB, f2,  1.0, 1.0, 6.0, xsec_plots, maxsig, signals, prior);
	cout << "between: " << maxsig << "\t" << xsec_plots << endl;

	JeffreysPrior jeffprior2d(pgammaint, lo, up, Mj, 1);

	for(int l=0; l<maxsig-1; l++) {
		cout << sigNames[l] << endl;
		double point = jeffprior2d(sigNames[l]);
		cout << sigNames[l] << ": " << point << endl;
	}

}