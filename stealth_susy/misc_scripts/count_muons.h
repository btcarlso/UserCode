/*
 *  count_muons.h
 *  
 *
 *  Created by Benjamin Carlson on 9/24/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "/uscms/home/btcarlso/stealth_code/root_headers.h"
#include <iostream> 
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

Long64_t entries=-9; 

Float_t st=-9; 
Int_t nJets=-9; 
int nMuons;
TString FILENAME=""; 

std::map<TString,TFile*> FName; 
std::map<TString,double> Mu0; 
std::map<TString,double> Mu1; 
std::map<TString,double> Mu2; 
std::map<TString,double> Mug2;
std::map<TString,double> xs; 
std::map<TString,double> Ngen; 

void initialize_xs();
void initialize_tree(TTree *tree);
void open_file(TString name, TString full_name);
void event_loop(TTree *tree);
void analyze_file(TString file);
void print();
void count_muons();