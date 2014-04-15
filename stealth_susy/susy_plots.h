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

std::map<TString, TF1*> f1Name;

std::map<TString, float> Nevent; 
std::map<TString, float> NeventE; 


TString cms_pre = "CMS Preliminary"; 
TString cms_sim = "CMS Simulation"; 
TString lumi = "#sqrt{s} = 8 TeV, L_{int} = 19.6 fb^{-1}"; 
TLatex L1; 
TLatex L2; 
TLatex L1sim; 
int Cx=800,Cy=600; 

int nJetsmax=6;

void susy_plots();
void open_files();
void close();
void write();
void setcolor(TGraphAsymmErrors *gr, int num);
void compare_ratios(int num, int den,int mu);
void compute_ratios();

TString st_name(int nJ, int Mu, TString weight, TString file);

TString ratio_name(int num, int den, int mu, TString file);


void draw_frac(TString variable, TString canvasname);
void compute_sensitivity( int nJ,int mu, TString signal);
void compute_expected_obs(int nJ, int mu, int iRPV);
void compute_expected_obs(int iRPV);

void divide_hist(TH1F *data, TH1F *mc, double *x, int N);
void plot_ratios(int den, int mu);

double compute_contamination(int mu, int nJ, TString signal, float st1, float st2);
TString acc_name(int nJ, int mu);
void compute_acceptance(int nJ,int mu, int irpv);
void table_contamination(int mu);
void table_obs(int nJ, int mu);

void signal_distributions();

float correction_bins[5];
int NbinsCor=5; 

float correction_bins_stealth[5];
int NbinsCor_stealth=3; 

void add_corrected_bkg(int nJ, int mu);
void relative_systematic(int nJ, int mu); 

void initialze_correction_bins();
void bookHisto();
void fit_corrections();
void draw_corrections();

void SF_systematic(TString variable, TString file);
void frac_bkg(TString variable, TString file);
void contamination_sensitivity(TString variable, TString signal);
void draw_nJet_contamination_plots(TString signal);

double get_N(double st1, double st2, TH1F *h);
double get_NE(double st1, double st2, TH1F *h);

void plot_process_ratios(int den, int mu);
void plot_compare_muon(int num, int den);

void compute_correction(int nJ,double st1,double st2);
void compute_correction_stealth(int nJ, double st1, double st2);
void get_column(TString name, TString file, double st1, double st2, double *N, double *NE);
void get_column_stealth(TString name, TString file, double st1, double st2, double *N, double *NE);
void transfer_factor();
void corrections_table();
void print_numbers(int nJ, float st1, float st2);
void rescaleMC_nJ(TString variable);
void rescaleMC(int nJ, int mu);
double scaleStBins(TString hist_name,double st1, double st2, double Corr);

void quantile_ratio(TH1F *h_num, TH1F *h_den, TH1F *h_num_uW, TH1F *h_den_uW, TString ratio_name);
void combine_histograms();
void R_iterative(TString on_shell, TString off_shell, TString signal,int nJ);
void print_backgrounds(TString variable,bool blind);
void fill_stack(TString variable, bool blind);
void CreateStack(TString Name,TString Title);
void CreateCanvas(TString Name,TString Title, int x, int y );

void jet_meanpt();


void Poisson(double N, double &mu_p,double &mu_m);
void setHistogramErrors(TH1F *h);

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