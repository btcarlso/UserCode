/*
 *  ;
 *  
 *
 *  Created by Benjamin Carlson on 11/15/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

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

std::map<TString, TH1F*> hName;
std::map<TString, TGraphAsymmErrors*> grName; 
std::map<TString, TH2F*> hName2D;
std::map<TString, TCanvas*> CName; 
std::map<TString,TFile*> fName;
std::map<TString,THStack*> stackName;
std::map<TString, float> Ngen; 
std::map<TString, float> MC_corr; 
std::map<TString, float> Nevent; 


TString cms_pre = "CMS Preliminary"; 
TString cms_sim = "CMS Simulation"; 
TString lumi = "#sqrt{s} = 8 TeV, L_{int} = 19.6 fb^{-1}"; 
TLatex L1; 
TLatex L2; 
TLatex L1sim; 
int Cx=800,Cy=600; 

void susy_plots();
void open_files();
void write();
void setcolor(TGraphAsymmErrors *gr, int num);
void compare_ratios(int num, int den,int mu);
void compute_ratios();

TString st_name(int nJ, int Mu, TString weight, TString file);
TString bkg_systematic(int nJ, int bkg_nJ, int mu, TString mode, TString file);
void fractional_uncertainty(int bkg_nJ, int mu, TString file);
TString ratio_name(int num, int den, int mu, TString file);

TString sig_bkg(int nJ,int mu,TString file);
void draw_frac(TString variable, TString canvasname);
void compute_sensitivity( int nJ,int mu, TString signal);

void plot_ratios(int den, int mu);

TString expected_bkg(int nJ, int bkg_nJ, int mu, TString mode, TString file);
void bkg_est(int nJ, int bkg_nJ, int mu, TString file);
void overlay_bkg_expectations(int nJ, int mu, TString file);
void plot_expected(int nJ, int bkg_nJ, int mu, TString file);

TString acc_name(int nJ, int mu);
void compute_acceptance(int nJ,int mu, int irpv);

void frac_bkg(TString variable, TString file);

double get_N(double st1, double st2, TH1F *h);
void compute_correction(int nJ,double st1,double st2);
void get_column(TString name, TString file, double st1, double st2, double *N);
void corrections_table();
void print_numbers();

void quantile_ratio(TH1F *h_num, TH1F *h_den, TH1F *h_num_uW, TH1F *h_den_uW, TString ratio_name);
void combine_histograms();
void fill_stack(TString variable, bool blind);
void CreateStack(TString Name,TString Title);
void CreateCanvas(TString Name,TString Title, int x, int y );
void open_file(TString name);
void GetHistograms(TString file);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void draw_headersim();
void draw_header();