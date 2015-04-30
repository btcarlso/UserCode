/*
 *  fit_z.h
 *  
 *
 *  Created by Benjamin Carlson on 2/25/14.
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
std::map<TString,TFile*> fName;
std::map<TString,TCanvas*> CName;

ofstream output;
ofstream outputR;
ofstream output_injection;
ofstream output_mc;

TString cms_pre = "CMS Preliminary";
TString cms_sim = "CMS Simulation";
TString lumi = "#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}";
TLatex L1;
TLatex L2;
TLatex L1sim;

void fit_z();
void draw_header();

double integralE(TH1F *h,double x1,double x2);
double integral(TH1F *h,double x1,double x2);
TString hist_name(TString sample, int nJets, int btags);
void open_file(TString name);
void GetHistograms();
double fit_fraction(int nJet, int nbtags, double &unc);
void CreateCanvas(TString Name,TString Title, int x, int y );


