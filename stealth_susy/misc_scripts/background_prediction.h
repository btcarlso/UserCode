//
//  background_prediction.c
//  
//
//  Created by Benjamin Carlson on 4/4/14.
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

std::map<TString,TFile*> fName;
std::map<TString, TCanvas*> CName;



void background_prediction();
void open_file(TString name);
void read_histograms(TString file);
void clone_histogram(TString name1, TString clone_name);
void combine_jetbins(TString file);
TString histname(int nbtag, TString file);
TString histnameP(int nbtag, TString file);
TString histnameM(int nbtag, TString file);
TString histname(int nbtag,int nJ, TString file);
TString histnameP(int nbtag,int nJ, TString file);
TString histnameM(int nbtag,int nJ, TString file);
double transfer_cor(double *x, double *par);
void transfer_factor(int nb_in, int nb_out, int nj_in, int nj_out, TString file);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, Double_t xLow, Double_t xUp,
					 Int_t nBinsY,Double_t yLow, Double_t yUp);
void CreateCanvas(TString Name,TString Title, int x, int y );
double transfer_cor_nJ(double *x, double *par);
void transfer_factor_st(int nb_in, int nb_out, int nj_in, int nj_out, TString file);
void subtract_smallbkg();
int iBin=1;