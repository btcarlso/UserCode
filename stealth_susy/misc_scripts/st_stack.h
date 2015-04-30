/*
 *  st_stack.h
 *  
 *
 *  Created by Benjamin Carlson on 7/18/13.
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
std::map<TString, TH2F*> hName2D;
std::map<TString, TCanvas*> CName; 
std::map<TString,TFile*> fName;
std::map<TString,THStack*> stackName;


void CreateCanvas(TString Name,TString Title, int x, int y );
void CreateStack(TString Name,TString Title);

void GetHistogram(TString hist_name, TString file);
void plot_npv();
void plot_signal();
void st_stack();
void open_file(TString name);
void open_files();
void fill_stack(TString variable, bool blind);
void fill_stack_inclusive(TString variable, bool blind);
void GetHistogram(TString hist_name, TString file);