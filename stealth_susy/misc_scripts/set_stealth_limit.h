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



void set_stealth_limit();
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
void draw_model_independent();
void model_Ind_Limit();
void draw_sensitivity(int mass);
void draw_header();

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

void CMS_lumi( TPad* pad, int iPeriod, int iPosX );

void CreateCanvas(TString Name,TString Title, int x, int y );
void draw_plots();
double lumi=19.7;