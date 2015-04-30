/*
 *  acceptance.h
 *  
 *
 *  Created by Benjamin Carlson on 10/29/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "/uscms/home/btcarlso/stealth_code/root_headers.h"

std::map<TString,TFile*> fName;
std::map<TString, TH1F*> hName;
std::map<TString, TCanvas*> CName; 

int Cx=700; 
int Cy=600; 


void acceptance();
void open_file(TString name);
void compute_acceptance(int mu);
void compute_sensitivity(int mu, int nJ, string signal);
void sensitivity_table(int mu);
void plot_sensitivity();
void plot_contamination();
void compute_contamination_Z();
void open_files();
void load_histograms();
void GetHistogram(TString hist_name, TString file);
void CreateCanvas(TString Name,TString Title, int x, int y );