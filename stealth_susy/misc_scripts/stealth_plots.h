//
//  stealth_plots.h
//  
//
//  Created by Benjamin Carlson on 4/10/14.
//
//

#include "/uscms/home/btcarlso/stealth_code/root_headers.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "Math/QuantFuncMathCore.h"
#include "TVirtualFitter.h"

#ifndef __CINT__
#include "Riostream.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#endif

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooExtendPdf.h"
#include "RooVoigtian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooRandomizeParamMCSModule.h"
#include "RooChi2Var.h"
#include "RooKeysPdf.h"
#include "RooGExpModel.h"
#include "RooNovosibirsk.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooMinuit.h"
#include "RooChebychev.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooRealSumPdf.h"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#include "RooCFunction2Binding.h"
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h"



using namespace RooFit ;

std::map<TString, TH1F*> hName;
std::map<TString, TH2F*> hName2D;
std::map<TString,THStack*> stackName;
std::map<TString, TGraphAsymmErrors*> grName;
std::map<TString,TFile*> fName;
std::map<TString, TCanvas*> CName;

bool overlay_stlabels=false;
bool plot_systematics=false;
bool expected_bkg=false;
bool scaleUpDown=false; 
bool notblind=true;
bool drawQCD=false;
bool ErrorX=false; 

TString cms_pre = "CMS Preliminary";
TString cms_sim = "CMS Simulation";
TString lumi = "#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
TLatex L1;
TLatex L2;
TLatex L1sim;

//stuff for lumi text

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = false;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );
//end lumi stuff 

void stealth_plots();
int nJetmax=7;

void open_file(TString name);
void combine_histograms(TString variable_name,std::vector<TString> names, TString outname);
void clone_histogram(TString name1, TString clone_name);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, Double_t xLow, Double_t xUp,
					 Int_t nBinsY,Double_t yLow, Double_t yUp);
void CreateStack(TString Name,TString Title);
void CreateCanvas(TString Name,TString Title, int x, int y );
void read_histograms(TString file);
void draw_correction(int nbsig,int ist);
void draw_correction(int nbsig);
double poisson(double n);

void compute_DY(int nbtag);
void compute_QCD(int nbtag);
int category_bin(int nJ, int ist);
void model_independent_datacard(TString variable, int nJet, int nbtag);
void compute_systematic(std::vector<TString> names);
void compute_systematic_MC(TString variable,std::vector<TString> names);
void normalize_systematic_categories(TString variable_name, TString sys); 
void blind_categoryplots(int nb, TString scale);
void fill_stack_eventcategories(TString variable_name, TString data_name, std::vector<TString> signal_names, std::vector<TString> names, std::vector<TString> legend_names, TString MC_name);
void table_simulation(int nbtag, int ist);
void table(int nbtag,int ist);
void shape_difference(int nbtag);
void plot_zbi(int nbsig, int nJ, TString stType, TString signal);
void summarize_zbi(int nbsig,TString stType);
double getZBI(double sig, double bkg, double bkge);
void plot_MC_exp(int nbtag);
void acceptance_table();
TString print_zbi(double mass, TString variable_name, TString signal, double st_used);
void print_loose_events();
void datacards(TString variable_name, TString data_name, TString signal_names, TString MC_name);
void datacards_new(TString variable_name, TString variable_name_cntrl, TString signal_names);
void sum_histograms();
void compute_ratio(TString variable_name, TString MC_name);
void compute_ratio(TString variable_name, TString MC_name, int N, double *x);
void CreateGraph(TH1F *h,TString graphName);
void draw_header();
void plot_categories();
void plot_distributions();
void plot_jet_distributions();
double compute_R0(double Ndata, double *N);
void normalization(int nbsig, TString tag, TString scale);
void predict_background(int nbsig,TString tag, TString scale);
void predict_background(int nbsig, int nJ,TString stType);
void plot_categories_expected();
void prediction_systematic(int nbsig, TString tag);
void predict_systematic(int nbsig, int nJ,TString stType);
void plot_mass();
void print_fractions();

