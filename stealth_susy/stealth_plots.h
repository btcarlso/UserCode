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

std::map<TString, TH1F*> hName;
std::map<TString, TH2F*> hName2D;
std::map<TString,THStack*> stackName;

std::map<TString,TFile*> fName;
std::map<TString, TCanvas*> CName;

void stealth_plots();

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

void blind_categoryplots(int nb, TString scale);
void fill_stack_eventcategories(TString variable_name, TString data_name, std::vector<TString> signal_names, std::vector<TString> names);

void compute_ratio(TString variable_name, TString MC_name);




