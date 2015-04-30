/*
 *  ratios.h
 *  
 *
 *  Created by Benjamin Carlson on 9/20/13.
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


void ratios();
void GetHistogram(TString hist_name, TString file);
void open_files();
void open_file(TString name);
void get_ratios();
void GetGraph(TString graph_name, TString file);
void CreateCanvas(TString Name,TString Title, int x, int y );
void scale_graph(TGraphAsymmErrors *gr, double SF);
double avg (TGraphAsymmErrors *gr);
void compare_ratios(int num, int den,int mu);
void compare_ratios_muon(int num, int den);
void setcolor(TGraphAsymmErrors *gr, int num);
void plot_ratios(TString file, int mu);