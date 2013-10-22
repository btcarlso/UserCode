/*
 *  Generate_LS.h
 *  
 *
 *  Created by Benjamin Carlson on 7/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */
#include "root_headers.h"
//#include "/Users/carlsonbt1/stealth_susy/code/root_headers.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "Math/QuantFuncMathCore.h"
#include "TVirtualFitter.h"
#include "smearing_parameters.c"
#include "bins_final.card"

std::map<TString,TCanvas*> CName; 
std::map<TString,TH1F*> hName; 

double R=0; 

TFile *input_file;
TFile *gen_file; 
TFile *output_file;
TTree *tree; 
TH2F *dm_m; 
TH1F *dm; 
TH1F *m; 
TH1F *M_gen; 

TF1 *mass_function;

void Generate_LS(int job);
void open_file(TString name);
void initialize_tree();
double PDF_shape(double *x, double *par);
double get_dm(TH1F *dm_distribution);
void LS(int ups, int iy, int ipt);
double QED_Mass(double *x, double *par);
void fit_QED_mass(int ups);
void loop(int iy, int ipt);
void get_histogram(int ups, int iy, int ipt);
void CreateCanvas(TString Name,TString Title, int x, int y );
