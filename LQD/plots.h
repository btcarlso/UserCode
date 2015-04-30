//
//  plots.c
//  
//
//  Created by Benjamin Carlson on 11/24/14.
//
//

#include <stdio.h>


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

int M=800;

bool drawLogo      = false;
bool plotSystematic=false;
bool right_=true;
void plots();
void open_file(TString name);
void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );
void fill_stack(TString variable_name, TString data_name, std::vector<TString> signal_names, std::vector<TString> names, std::vector<TString> legend_names, TString MC_name);
void rescale();
void ratio();
void combine_histograms(TString variable_name,std::vector<TString> names, TString outname);
void clone_histogram(TString name1, TString clone_name);
void combine_all();
void CreateGraph(TH1F *h,TString graphName);
double poisson(double n);
void plot_distribution();
void CreateStack(TString Name,TString Title);
void CreateCanvas(TString Name,TString Title, int x, int y );
void read_histograms(TString file);

