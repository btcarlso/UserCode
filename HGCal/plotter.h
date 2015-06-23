//
//  plotter.h
//  
//
//  Created by Benjamin Carlson on 8/26/14.
//
//

#include "/uscms/home/btcarlso/stealth_code/root_headers.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "TObjArray.h"
#include "TDatime.h"
#include "TVirtualFitter.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector2.h"

int Nheat=16;
int Nsens=24;
std::vector<TGraph*> grList;
std::map<TString, TVector2*> pos;
std::map<TString, TGraph2D*> grName2D;
std::map<TString, TH2F*> hName2D;
std::map<TString, TCanvas*> CName;

void plotter();
int index(int imeas, int index);
int heaterIndex(int imeas, int iHeat);
int sensIndex(int imeas, int iSens);
TDatime StringsToTime(TString date, TString time);
void CreateGraph(TString name);
void CreateGraph2D(TString name);
void map_pos(TString name, double x, double y);
void map_all();
void CreateCanvas(TString name, int Cx,int Cy);