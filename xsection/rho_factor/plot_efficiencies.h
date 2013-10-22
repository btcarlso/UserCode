/*
 *  plot_efficiencies.h
 *  
 *
 *  Created by Benjamin Carlson on 7/12/13.
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

std::map<TString, TH1F*> hName;
std::map<TString, TProfile*> prName;
std::map<TString, TH2F*> hName2D;
std::map<TString, TCanvas*> CName;           // Map for histograms

int UPS; 

TString cms_pre;
TLatex L1; 

TFile *input_file=new TFile("/uscms_data/d3/btcarlso/Ups1S_Analysis/onia1S_Acceptance_tree.root","READ");
							///uscms/home/btcarlso/CMSSW_4_2_9_HLT1_hltpatch1/src/onia2MuMu_tree.root","READ");
							///uscms_data/d3/btcarlso/Ups1S_Analysis/onia_Acceptance_tree.root","READ");
							

TFile *output_canvas;// = new TFile("output_file_efficiencies1S.root","RECREATE");
void draw_header();
void CreateCanvas(TString Name,TString Title, int x, int y );
void plot_efficiencies();
void plot_2Deta();
void plot_vertex_efficiency();
void remove_initial_points(TGraph *gr);
void Rebin(TString name);
void plot_polarization(); 
void compute_acceptance(string method, int iy);
void compute_acceptance(string method);
void plot_acceptance(int iy);
void rho_table(TGraphAsymmErrors *gr);
void draw_histograms(TString name, int dim);
void reweight_rho(TH2F *rho, string method);
void divide_histograms(TString name, int dim);
void GetProfile(TString name);
void GetHistogram(TString name);
void GetHistogram2D(TString name);
void CreateCanvas(TString Name,TString Title, int x, int y );
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void projectionY(TString name, float y1, float y2);
void load_histograms();