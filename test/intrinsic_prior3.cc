//-----------------------------------------------------------------------------
// File:        data_likelihood.C
// Description: Generate likelihood function in the form of a TGraph using
//		the data histogram.
//				discriminant histogram
//-----------------------------------------------------------------------------
// Created:     11-04-2012  Joe Bochenek

//-----------------------------------------------------------------------------

#include <time.h>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TMath.h>
#include <TTree.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <TGaxis.h>
#include <gsl/gsl_randist.h>
#include <TGraph2D.h>
#include <gsl/gsl_sf_gamma.h>
#include <TRandom.h>
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"

using namespace std;
#include "samples.h"

//#include "samples.h"
#include "IntrinsicPrior.h"
#include "PoissonGammaIntegral.h"



double Min( double x, double y )
{
    return ( x > y ) ? y : x;
}


int combinations_i(std::vector< std::vector<int> > &data, std::vector<int> & histogram, int total, int bin, int totalbins, int bincap, int count)
{
		if( bin > (totalbins - 2) ) {
			if(total < bincap + 1) {
				histogram[bin] = total;
				data.push_back(histogram);
				count++;
				if(!(count%1000000)) cout << count << endl;
//				cout << "bin: " << bin << "/" << totalbins << "\t total: " << total << "\tcount: " << count << endl;
			}
			return count;
		} else
		{
			for(int thisbincount = Min(bincap, total); thisbincount >= 0; thisbincount--){
				histogram[bin] = thisbincount;
//				cout << "\t" <<  thisbincount;
				count = combinations_i(data, histogram, (total-thisbincount), (bin+1), totalbins, bincap, count);
			}
			return count;
		}
}


int combinations(std::vector< std::vector<int> > &data, int max, int totalbins, int bincap, int count)
{
	for(int total = 0; total < max + 1; total++)
	{
//		cout << "total: " << total << endl;
		std::vector<int> histogram(totalbins);
		count = combinations_i(data, histogram, total, 0, totalbins, bincap, count);		
//		cout << "===========" << endl;
	}
	return count;
}




int main(){
	//---------------------------------------------------------------------------
	// Make Vectors/Histograms
	//---------------------------------------------------------------------------

	const int xsec_plots = 50;
	double lumi_scale = 1.0;
	double this_xsec_factor = 1.0;
	int thismass = 13;
    ROOT::Math::Random<ROOT::Math::GSLRngMT> _gslRan;

	TString scripttag = "input1D";	
	style();
	
	//---------------------------------------------------------------------------
	// RETRIEVE HISTOGRAMS
	//---------------------------------------------------------------------------
	TH2F* histos_sig[maxsig];

	TFile file("output/hist_iso.root");
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

	// Print some histograms for debugging
	

	char filename[256];
	char charname [] = "bnnvsmt_bkg";
	sprintf(filename, "output/all/%s_data.pdf", charname);
	TCanvas *c = new TCanvas(charname,"Average Likelihood vs. Mass",200,10,700,500);
	histos_bkg->Draw("colz");
	c->Update();
	c->Modified();
	cout << filename << endl;
	c->SaveAs(filename);
	delete c;






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
//		histos_sig[k]->Rebin(25);
		histos_sig[k]->Scale(lumi_scale);
		signal_weights[k] = histos_sig[k]->Integral();
		cout << sigNames[k] << ": " << signal_weights[k] << endl;
		cout << "signal[" << k << "]: " << sigNames[k] << "\t" << signal_weights[k] << endl;
		histos_sig[k]->Sumw2();
	}


	//---------------------------------------------------------------------------
	// FILL BACKGROUND VECTORS
	//---------------------------------------------------------------------------

	histos_bkg->Scale(lumi_scale);			
	// Calculate errors in bins for the backgrounds
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
		//				cout << hmlp[j]->GetBinContent(i) << " +/-" << hmlp[j]->GetBinContent(i)/sqrt(hmlp[j]->GetBinError(i))  << endl;
	}
	}


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
	
	cout << "Z+jets Integral: " << histos_zjets->Integral() << endl;
	cout << "ZZ Integral: " << histos_bkg->Integral() << endl;
	cout << "Signal: " << histos_sig[thismass]->Integral() << endl;
	


	//---------------------------------------------------------------------------
	// Remake  A, f vectors for IntrinsicPrior algorithm
	//---------------------------------------------------------------------------

	vector<double> signals;		

	
	// Loop through m_H and \mu to calculate 2D prior and store the grid in a TGraph2D
	double xsec_factors[xsec_plots]; 
	double count = 0;

    TGraph2D *lh_tgraph = new TGraph2D();
    TGraph2D *lh_tgraph_norm = new TGraph2D();

	clock_t start, end;
	start = clock();

	char charname2 [] = "bnnvsmt_bkg";
	sprintf(filename, "output/all/%s_data.pdf", charname2);
	c = new TCanvas(charname2,"Average Likelihood vs. Mass",200,10,700,500);
	histos_bkg->Draw("colz");
	c->Update();
	c->Modified();
	cout << filename << endl;
	c->SaveAs(filename);
	delete c;

	char filename_s[256];
	char charname_s[] = "bnnvsmt_sig";
	sprintf(filename_s, "output/all/%s_data.pdf", charname_s);
	c = new TCanvas(charname_s,"Average Likelihood vs. Mass",200,10,700,500);
	histos_sig[thismass]->Draw("colz");
	c->Update();
	c->Modified();
	cout << filename_s << endl;
	c->SaveAs(filename_s);
	delete c;
	


//	int bins = 5;
	vector< vector<int> > data;
//	combinations(data, 30, bins, 6, 0);
	int size = data.size();
	cout << "Size: " << size << "\t COunt: "<< endl;
	double likelihood_ave[maxsig][xsec_plots];
	for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++)  likelihood_ave[l][j] = 0;
	double final_integral = 0;
	double prob_null_integral = 0;


	int iter = 100;

	// Load Reference prior
	const int Npoints = 1500;
	TFile f2("output/refprior_1000.root");
	TGraph2D *refprior = (TGraph2D*) gDirectory->Get("prior");
	Double_t rppoints[Npoints];
	Double_t *rpPoints = refprior->GetZ();
	for(Int_t l=0; l < maxsig; l++) {
	  cout << "Mass: " << sigNames[l] << "\t";
	  for(Int_t j=0; j < xsec_plots; j++)  cout << " : " << rpPoints[(j+l*xsec_plots)]<< "\t";
	  cout << endl;
	}


	gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);

	// Loop through data permutations
	for(Int_t i_d=0; i_d<iter; i_d++)
	{
		double likelihood_sig[maxsig][xsec_plots]; 
		double likelihood_null[maxsig][xsec_plots]; 

		double l_integral_sig = 0;
		vector<double> D;
		D.clear();
		int this_size = 0;
		
		double cdf = 1.0;

		int _nbins = mva_src_zjets.size();
		for(int ibin=0; ibin<_nbins; ++ibin)
		{
			double bmean = 0;
			bmean += _gslRan.Gamma(mva_error_b[ibin]*mva_src_b[ibin], 1.0)/mva_error_b[ibin];
			bmean += _gslRan.Gamma(mva_error_zjets[ibin]*mva_src_zjets[ibin], 1.0)/mva_error_zjets[ibin];
			double sample = _gslRan.Poisson(bmean);
			D.push_back(sample);
			cdf *= gsl_ran_negative_binomial_pdf(D[ibin], 0.5, bmean);
//			std::cout << D[ibin] << endl;
		}

		prob_null_integral += cdf;
		cout << i_d << "\t" << cdf << endl;


		// First loop through signal mass points
		for(Int_t l=0; l<maxsig; l++) {
//			cout << endl;					
			// Next loop through signal strength modifier points
			for( int j = 0; j < xsec_plots; j++)
			{						
				double xsecfactor =  (j)*0.1;
				xsec_factors[j] = xsecfactor;

				vector< vector<double> > A;
				vector< vector<double> > f;
				vector<double> p;		
				vector<double> mva_src_s;
				vector<double> mva_error_s;

//				cout << "i_d: " << i_d << "/" << size << "\tSignal: m_H = " << sigNames_mc[l] << " GeV,\t mu = " << xsecfactor;

				// Build signal vectors to pass to intrinsicprior at this mass point
				// SIGNAL
				for(int i=1; i < xbins+1; i++)
				{
				for(int j=1; j < ybins+1; j++)
				{
					double bincontent = xsecfactor * histos_sig[l]->GetBinContent(i,j);
					double binerror = xsecfactor * histos_sig[l]->GetBinError(i,j);				
					double error = 0. ;
					if(bincontent < 0) bincontent = 0;
	//				if(bincontent > 0) error = double( binerror*binerror/bincontent );		
					if(bincontent > 0) error = 50 / bincontent;
	
	//				cout << sigNames[l] << "\txsecfactor: "  << xsecfactor << "\tbincontent: " << bincontent << "\t eff events: " << bincontent * 25000 << "\t error: " << error << "\t rel error: " << 1/sqrt(25000) << endl;				
	
					mva_src_s.push_back(bincontent);
					mva_error_s.push_back(error);
				}
				}
				
				A.push_back(mva_src_s);
				f.push_back(mva_error_s);
	
				A.push_back(mva_src_b);
				f.push_back(mva_error_b);	
	
				A.push_back(mva_src_zjets);
				f.push_back(mva_error_zjets);	
	
				p.push_back( 1.0 );
				p.push_back( 1.0 );
				p.push_back( 1.0 );

				double prob_sig = exp(pg::poissongamma(D, p, A, f, 1.0, true, false));
//				double prob_sig = exp(pg::poissongamma(D, p, A, f, true, true));
//				cout << (j+l*xsec_plots) << " : " <<  rpPoints[((j-1)+l*xsec_plots)] << endl;
				likelihood_sig[l][j] = prob_sig * rpPoints[(j+l*xsec_plots)];
				l_integral_sig += rpPoints[(j+l*xsec_plots)] * prob_sig * 0.2;
				lh_tgraph->SetPoint( (j+l*xsec_plots), sigNames[l], xsecfactor, prob_sig );
				p.clear();
				f.clear();		
				A.clear();
				mva_src_s.clear();
				mva_error_s.clear();
			}
		}

		// Do NULL (mu = 0) likelihood
		vector< vector<double> > A2;
		vector<double> p2;
		vector< vector<double> > f2;

		for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++) likelihood_sig[l][j] /= l_integral_sig;
		for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++) likelihood_ave[l][j] += likelihood_sig[l][j];
		
		if(!(i_d%100)) cout << "i_d: " << i_d << "/" << iter << "\t sig lh int: " << prob_null_integral << endl;
		D.clear();
	}
	
	

	for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++) likelihood_ave[l][j]/=(iter) ;
	for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++) lh_tgraph_norm->SetPoint( ((j-1)+l*xsec_plots), sigNames[l], xsec_factors[j], likelihood_ave[l][j] );
//	for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++) cout << ((j-1)+l*xsec_plots) << "\t" << sigNames_mc[l] << "\t" << xsec_factors[j] << "\t" << likelihood_ave[l][j] <<  endl;
	for(Int_t l=0; l < maxsig; l++) for(Int_t j=0; j < xsec_plots; j++) final_integral += likelihood_ave[l][j] * 0.2;
	cout << "null dist integral: " << prob_null_integral << "\tfinal integral: " << final_integral << endl;

	count+= 1.0;
	end = clock();
	
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";





	TFile outfile("output/int_prior.root", "RECREATE");
	outfile.cd();
	lh_tgraph_norm->SetName("intrinsic_prior");	
	lh_tgraph_norm->Write();
	outfile.Close();	


	//---------------------------------------------------------------------------
	// Posterior Plot
	//---------------------------------------------------------------------------


	std::stringstream name2;
	name2 << "intrinsic_prior_norm3";

	std::stringstream plotname2;
	plotname2 << "output/" << name2.str() << ".png";

	std::stringstream title;
	title << "Intrinsic Prior, m_{H} vs #mu";

	// Make plot of likelihood
	TCanvas *c4 = new TCanvas(name2.str().c_str(), title.str().c_str(), 200,10,700,500);
	TPad *pad5 = new TPad("pad5","",0,0,1.0,1);
	pad5->Draw();
	pad5->cd();

//	TGaxis: :SetMaxDigits(2);
	c4->SetRightMargin(0.16); 
	pad5->SetRightMargin(0.16); 


	lh_tgraph_norm->SetTitle(title.str().c_str());
    lh_tgraph_norm->GetXaxis()->SetTitleFont(42);
	lh_tgraph_norm->GetXaxis()->SetTitle("m_{H} (GeV)");
    lh_tgraph_norm->GetZaxis()->SetTitleFont(42);
	lh_tgraph_norm->GetZaxis()->SetTitle("#pi_{0}(m_{H},#mu)");
    lh_tgraph_norm->GetYaxis()->SetTitleFont(42);
	lh_tgraph_norm->GetYaxis()->SetTitle("#mu");

	c4->SetRightMargin(1.3);
	gStyle->SetPalette(1);

	lh_tgraph_norm->Draw("colz");

	c4->Update();
	c4->Modified();
	c4->Write();
	c4->SaveAs( plotname2.str().c_str() );

	name2.str("");
	plotname2.str("");
	title.str("");
}
