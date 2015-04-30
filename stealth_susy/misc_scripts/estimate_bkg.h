/*
 *  estimate_bkg.h
 *  
 *
 *  Created by Benjamin Carlson on 10/30/13.
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
std::map<TString, TGraphAsymmErrors*> grName; 
std::map<TString, TH2F*> hName2D;
std::map<TString, TCanvas*> CName; 
std::map<TString,TFile*> fName;
std::map<TString,THStack*> stackName;

int Cx=700; 
int Cy=600;

void estimate_bkg();
void fit(int mu);
void bkg_est(int nJ, int mu,string file);
void open_files();
void open_file(TString name);
void get_ratios();
void get_hists();
void GetGraph(TString graph_name, TString file);
void GetHistogram(TString hist_name, TString file);

void CreateCanvas(TString Name,TString Title, int x, int y );
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins);
