/*
 *  set_limit.h
 *  
 *
 *  Created by Benjamin Carlson on 1/15/14.
 *  Copyright 2014 Carnegie Mellon University. All rights reserved.
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
void set_limit();
std::map<TString, TH1F*> hName;
std::map<TString, TCanvas*> CName;

std::map<TString, TGraphAsymmErrors*> grName; 
void create_theory();

void read_data();
void CreateGraph(TString name);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void bookGraphs();
void GetHistograms();

void CreateCanvas(TString Name,TString Title, int x, int y );
void draw_plots(int nJ, int mu);
void draw_plots();
double lumi=19.7;