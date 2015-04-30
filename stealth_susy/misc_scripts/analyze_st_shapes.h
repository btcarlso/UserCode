/*
 *  analyze_st_shapes.h
 *  
 *
 *  Created by Benjamin Carlson on 7/5/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "root_headers.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "Math/QuantFuncMathCore.h"
#include "TVirtualFitter.h"

TFile *input_file; 

//photonHad_All_1mu

TString type="ttSemiLept";

TFile *output_file = new TFile("../Figures/canvas_plots_"+type+".root","RECREATE"); 

std::map<TString, TCanvas*> CName;           // Map for histograms


double pol_scaled(double *x, double *par);
void plot_ratios();
double power_law(double *x, double *par); 
void fit(); 
void fit_ratio_pTmin();
void plot_ratio_pTmin();
void draw_muon_pt();
void ht_st();
void draw_jet_pt();
void st_peak();
void st_overlay();
void energy_fractions();
void st_pTmin();
void shift_histogram(TH1F *h,float deltaX);
void shape_comparison();
double shifted_shape(double *x, double *par);
void CreateCanvas(TString Name,TString Title, int x, int y );
TH1F *hist; 