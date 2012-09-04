#ifndef SAMPLES_H
#define SAMPLES_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TROOT.h>
#include <algorithm>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TH2.h>

using namespace std;


//const TString dirIN_tmva =      "/home/jbochenek/data/HWWSep28-2011/tmva/";
const TString dirIN =		"/home/jbochenek/data/HZZ4lNov6-2011/hzz/";
//const TString dirIN =		"/storage/5/home/jbochenek/data/HZZ4l/";
//const TString dirIN =           "/storage/5/home/jbochenek/data/hzz_nomasscut/";
//const TString dirIN =           "/storage/5/home/jbochenek/data/hzz/";


//-----------------------------------------------------------------------------
// Utilities
//-----------------------------------------------------------------------------


// Poissonian distribution without factorial
float poissonnofact(float mean,float nobs)
{
float p=0;
p=exp(-mean)* pow(mean,nobs);
return p;
}


// likelihood ration Q=L(s+b)/L(b)
// hyphotesis used: nobs=signal+bkg instead of toy experiment
float likelihoodRatio(float signal,float bkg,float nobs)
{
float q=0;
//if (poisson(bkg,nobs)!=0) q= poisson(signal+bkg,nobs)/poisson(bkg,nobs);
if (poissonnofact(bkg,nobs)!=0) q= poissonnofact(signal+bkg,nobs)/poissonnofact(bkg,nobs);
return q;
}


// Significance Sl=sqrt (2lnQ)
// hyphotesis used: nobs=signal+bkg instead of toy experiment
float significance(float signal,float bkg,float nobs){
float sl=0;
sl=sqrt(2*log(likelihoodRatio(signal,bkg,nobs)));
return sl;
}


void style()
{	
gROOT->SetStyle("Pub");
TStyle* style = gROOT->GetStyle("Pub");
style->SetFrameBorderMode(0);
style->SetCanvasBorderSize(0);
style->SetCanvasBorderMode(0);
style->SetCanvasColor(0);
style->SetPadBorderMode(0);
style->SetPadColor(0);
style->SetHistLineWidth(0.5);



// For the Global title:
style->SetOptTitle(1);
style->SetTitleFont(42);
style->SetTitleColor(1);
style->SetTitleTextColor(1);
style->SetTitleFillColor(10);
style->SetTitleFontSize(0.05);

// For the axis titles:
style->SetTitleColor(1, "XYZ");
style->SetTitleFont(42, "XYZ");
style->SetTitleSize(0.05, "XYZ");

//gStyle->SetOptStat(0);

// For the axis labels:
style->SetLabelColor(1, "XYZ");
style->SetLabelFont(42, "XYZ");
style->SetLabelSize(0.05, "XYZ");

// For the axis:
style->SetAxisColor(1, "XYZ");
style->SetStripDecimals(kTRUE);
style->SetTickLength(0.03, "XYZ");
style->SetNdivisions(505, "XYZ");
style->SetPadTickX(1);  
style->SetPadTickY(1);

//gROOT->ForceStyle();
}


double modul(double a, double b)
{
int result = static_cast<int>( a / b );
return a - static_cast<double>( result ) * b;
}

void formatHist(TH1* h, int color)
{
	//int NDIVX =-505;	
	h->SetLineColor(color);
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitle("efficiency (BKG)");
	h->GetXaxis()->SetTitleOffset(1.2);
	//h->GetXaxis()->SetNdivisions(NDIVX);
	
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitle("efficiency (SIG)");
	h->GetYaxis()->SetTitleOffset(1.8);
}

std::vector<double> cdf(TH1* hist)
{
std::vector<double> c(hist->GetNbinsX());
c[0] = hist->GetBinContent(1);
for(int i=1; i < hist->GetNbinsX(); i++)
c[i] = c[i-1]+hist->GetBinContent(i+1);
return c;
}

std::vector<double> rcdf(TH1* hist)
{
std::vector<double> c(hist->GetNbinsX());
c[hist->GetNbinsX()-1] = hist->GetBinContent(hist->GetNbinsX());
for(int i=hist->GetNbinsX()-1; i > 0; i--)
c[i] = c[i]+hist->GetBinContent(i);
return c;
}



const int N = 43;
double histo_min[N] = {
0,	   //mu1_pt
0,	   //mu1_eta
0,	   //mu1_phi
0,	   //mu1_charge
0,	   //mu1_trackIso
0,	   //mu1_EcalIso
0,	   //mu1_HcalIso
0,	   //mu1_X
0,	   //mu1_SIP
0,	   //mu2_pt
0,	   //mu2_eta
0,	   //mu2_phi
0,	   //mu2_charge
0,	   //mu2_trackIso
0,	   //mu2_EcalIso
0,	   //mu2_HcalIso
0,	   //mu2_X
0,	   //mu2_SIP
0,	   //mu3_pt
0,	   //mu3_eta
0,	   //mu3_phi
0,	   //mu3_charge
0,	   //mu3_trackIso
0,	   //mu3_EcalIso
0,	   //mu3_HcalIso
0,	   //mu3_X
0,	   //mu3_SIP
0,	   //mu4_pt
0,	   //mu4_eta
0,	   //mu4_phi
0,	   //mu4_charge
0,	   //mu4_trackIso
0,	   //mu4_EcalIso
0,	   //mu4_HcalIso
0,	   //mu4_X
0,	   //mu4_SIP
0,	   //worst_iso_X_mu
0,	   //second_worst_iso_X_mu
0,	   //worst_vertex
0,	   //second_worst_vertex
40,		//mZ
0,	   //mZstar
110,	 //mbestH
};


string varnames[N] = {
	"p_{T,l1}",
	"#eta_{l1}",
	"#phi_{l1}",
	"charge_{l1}",
	"Iso_{track,l1}",
	"Iso_{ECAL,l1}",
	"Iso_{HCAL,l1}",
	"X_{l1}",
	"SIP_{l1}",
	"p_{T,l2}",
	"#eta_{l2}",
	"#phi_{l2}",
	"charge_{l2}",
	"Iso_{track,l2}",
	"Iso_{ECAL,l2}",
	"Iso_{HCAL,l2}",
	"X_{l2}",
	"SIP_{l2}",
	"p_{T,l3}",
	"#eta_{l3}",
	"#phi_{l3}",
	"charge_{l3}",
	"Iso_{track,l3}",
	"Iso_{ECAL,l3}",
	"Iso_{HCAL,l3}",
	"X_{l3}",
	"SIP_{l4}",
	"p_{T,l4}",
	"#eta_{l4}",
	"#phi_{l4}",
	"charge_{l4}",
	"Iso_{track,l4}",
	"Iso_{ECAL,l4}",
	"Iso_{HCAL,l4}",
	"X_{l4}",
	"SIP_{l4}",
	"Worst Isolated lepton X",
	"2nd Worst Iso lept X",
	"Worst Vertex",
	"2nd Worst Vertex",
	"m_{Z}",
	"m_{Z*}",
	"M_{4l}",
};

	

double histo_max[N] = {
150,	 //mu1_pt
5,	   //mu1_eta
5,	   //mu1_phi
2,	   //mu1_charge
0.1,	 //mu1_trackIso
0.1,	 //mu1_EcalIso
0.1,	 //mu1_HcalIso
0.1,	 //mu1_X
5,	   //mu1_SIP
150,	 //mu1_pt
3.5,	 //mu1_eta
3.5,	 //mu1_phi
2,	   //mu1_charge
0.1,	 //mu1_trackIso
0.1,	 //mu1_EcalIso
0.1,	 //mu1_HcalIso
0.1,	 //mu1_X
5,	   //mu1_SIP
150,	 //mu1_pt
3.5,	 //mu1_eta
3.5,	 //mu1_phi
2,	   //mu1_charge
0.1,	 //mu1_trackIso
0.1,	 //mu1_EcalIso
0.1,	 //mu1_HcalIso
0.1,	 //mu1_X
5,	   //mu1_SIP
100,	 //mu1_pt
3.5,	 //mu1_eta
3.5,	 //mu1_phi
2,	   //mu1_charge
0.1,	 //mu1_trackIso
0.1,	 //mu1_EcalIso
0.1,	 //mu1_HcalIso
0.1,	 //mu1_X
5,	   //mu1_SIP
0.2,	 //worst_iso_X_mu
0.1,	 //second_worst_iso_X_mu
30,		//worst_vertex
10,		//second_worst_vertex
200,	 //mZ
220,	 //mZstar
165,	 //mbestH  
};


double logscale[N] = {
1,	   //mu1_pt
0,	   //mu1_eta
0,	   //mu1_phi
0,	   //mu1_charge
0,	   //mu1_trackIso
0,	   //mu1_EcalIso
0,	   //mu1_HcalIso
0,	   //mu1_X
0,	   //mu1_SIP
1,	   //mu2_pt
0,	   //mu2_eta
0,	   //mu2_phi
0,	   //mu2_charge
0,	   //mu2_trackIso
0,	   //mu2_EcalIso
0,	   //mu2_HcalIso
0,	   //mu2_X
0,	   //mu2_SIP
1,	   //mu3_pt
0,	   //mu3_eta
0,	   //mu3_phi
0,	   //mu3_charge
0,	   //mu3_trackIso
0,	   //mu3_EcalIso
0,	   //mu3_HcalIso
0,	   //mu3_X
0,	   //mu3_SIP
1,	   //mu4_pt
0,	   //mu4_eta
0,	   //mu4_phi
0,	   //mu4_charge
0,	   //mu4_trackIso
0,	   //mu4_EcalIso
0,	   //mu4_HcalIso
0,	   //mu4_X
0,	   //mu4_SIP
1,	   //worst_iso_X_mu
1,	   //second_worst_iso_X_mu
1,	   //worst_vertex
1,	   //second_worst_vertex
1,	   //mZ
1,	   //mZstar
1,	   //mbestH
};


double binning[N] = {
25,		//mu1_pt
25,		//mu1_eta
25,		//mu1_phi
25,		//mu1_charge
25,		//mu1_trackIso
25,		//mu1_EcalIso
25,		//mu1_HcalIso
25,		//mu1_X
25,		//mu1_SIP
25,		//mu2_pt
25,		//mu2_eta
25,		//mu2_phi
25,		//mu2_charge
25,		//mu2_trackIso
25,		//mu2_EcalIso
25,		//mu2_HcalIso
25,		//mu2_X
25,		//mu2_SIP
25,		//mu3_pt
25,		//mu3_eta
25,		//mu3_phi
25,		//mu3_charge
25,		//mu3_trackIso
25,		//mu3_EcalIso
25,		//mu3_HcalIso
25,		//mu3_X
25,		//mu3_SIP
25,		//mu4_pt
25,		//mu4_eta
25,		//mu4_phi
25,		//mu4_charge
25,		//mu4_trackIso
25,		//mu4_EcalIso
25,		//mu4_HcalIso
25,		//mu4_X
25,		//mu4_SIP
25,		//worst_iso_X_mu
25,		//second_worst_iso_X_mu
25,		//worst_vertex
25,		//second_worst_vertex
25,		//mZ
25,		//mZstar
25,		//mbestH
};




// -------------  -------------  -------------  -------------  ------------- 
// SIGNAL LOOP
// -------------  -------------  -------------  -------------  ------------- 

string thischannel = "H#rightarrow ZZ#rightarrow 4l";

const int maxNum = 3;

double mvacut = 0.5;

int nbins = 15;
const int maxdat = 1;
const int maxbkg = 7;
const int maxsig = 50;
const int maxsig_plots = 6;


double lumifactor = 4.71;
int xbins = 50;
int ybins = 50;
double xmin = 0.0;
double xmax = 1.0;
double ymin = 112;
double ymax = 162;

int ybins_morph = 280;
int xbins_morph = 100;
double ymin_morph = 110;
double ymax_morph = 166;

double otherfactor =  1;
double  lumi_factor = 1;
double  scale_factor = 1.41/0.94;

TString datfiles[maxdat] = {
"all_data.dat",
};


TString bkgfiles[maxbkg] = {
"all_TT.dat",
"all_WZ.dat",
"all_ZBB.dat",
"all_ZZ.dat",
"DY_2e2mu.dat",
};

TString bkgNames[maxbkg] = {
"TT",
"WZ",
"Z+bb",
"ZZ",
"Z+cc",
"Z+light",
};



string sampleNames_sig = "sig_ggH";

string sampleNames[7] = {
"NOSYST",
"NOSYST",
"bkg_zjets",
"bkg_qqzz",
"bkg_zjets",
"bkg_zjets",
"bkg_ggzz"
};

TString sigfiles[maxsig] = {
//"all_H110.dat",
"all_H115.dat",
"all_H120.dat",
"all_H130.dat",
"all_H140.dat",
"all_H150.dat",
"all_H160.dat",
//"all_H170.dat",
//"all_H180.dat",
//"all_H190.dat",
//"all_H200.dat"
};

Double_t sigNames[maxsig];


double sigNames_mc[maxsig_plots] = {
115,
120,
130,
140,
150,
160
};


Double_t colors[maxbkg] = {
38,
5,
3,
4,
7,
6
};



Double_t bkgweight[maxbkg] = {
lumi_factor * scale_factor,		// qq WW
lumi_factor * scale_factor,		// gg WW
lumi_factor * scale_factor,		// DY mumu
lumi_factor * scale_factor,		// DY ee
lumi_factor * scale_factor * 3 / 45	//  DY MC was set with weights of 1.  Adjusting to the AN  
};


TString signalTag[maxsig] = {
//"H#rightarrow ZZ#rightarrow 4l (M_{H} = 110)",
"H#rightarrow ZZ#rightarrow 4l (M_{H} = 115)",
"H#rightarrow ZZ#rightarrow 4l (M_{H} = 120)",
"H#rightarrow ZZ#rightarrow 4l (M_{H} = 130)",
"H#rightarrow ZZ#rightarrow 4l (M_{H} = 140)",
"H#rightarrow ZZ#rightarrow 4l (M_{H} = 150)",
"H#rightarrow ZZ#rightarrow 4l (M_{H} = 160)",
//"H#rightarrow ZZ#rightarrow 4l (M_{H} = 170)",
//"H#rightarrow ZZ#rightarrow 4l (M_{H} = 180)",
//"H#rightarrow ZZ#rightarrow 4l (M_{H} = 190)",
//"H#rightarrow ZZ#rightarrow 4l (M_{H} = 200)",
};

int signalFill[maxsig] = {
//3006,
3004,
3005,
3002,
3003,
3001,
3004,
//3005,
//3002,
//3003,
//3001,
};

int signalColor[maxsig] = {
//2,
1,
2,
13,
40,
49,
6,
//7,
//8,
//9,
//12
};

double sigfactor = 1;

Double_t signalweight[maxsig] = {
//10,
1.5,
1.5,
1.5,
1.5,
1.5,
1.5,
//1.5,
//1.5,
//1.5,
//1.5
};

//---------------------------------------------------------------------------
// PARAMETERS
//---------------------------------------------------------------------------




#endif
