/*
 *  MC_smearing.h
 *  
 *
 *  Created by Benjamin Carlson on 9/26/12.
 *  Copyright 2012 Carnegie Mellon University. All rights reserved.
 *
 */



#include <TROOT.h>
#include <TEnv.h>
#include <TObject.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TMultiGraph.h> 
#include <TStopwatch.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TText.h>
#include <TPaveText.h>
#include <TAttMarker.h>
#include <TMinuit.h>
#include <TString.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TRandom3.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h> 
#include <TAxis.h>
#include <TPave.h>
#include <TPaveStats.h>

#include <iostream> 
#include <fstream>
#include <math.h>
#include <vector>

#include "/afs/cern.ch/work/b/bcarlson/public/common_code/bins_use.card" //bins_new.card is for the old binning, and switch to bins_use.card for new binning
#include "/afs/cern.ch/work/b/bcarlson/public/common_code/opt_params.C"
#include "/afs/cern.ch/work/b/bcarlson/public/common_code/names.C"


//string output_dir="/afs/cern.ch/work/b/bcarlson/private/smearing_out/new/";
string file_name=output_dir+FN;

TFile *MC= new TFile("/afs/cern.ch/work/b/bcarlson/public/data_files/rad_tail/mass_rad_tail.root","READ");
TF1 *mass_function;

TFile *outfile = new TFile(file_name.c_str(),"RECREATE");

TFile *data;//=new TFile("/afs/cern.ch/user/b/bcarlson/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Upsilon2MuMu/xsection/output_files/data_feb1.root","READ");///afs/cern.ch/work/b/bcarlson/public/data_files/cor15.root","READ");//data file 


void MC_smearing();
double get_dm(TH1D *h);
void peak(int ups);
void PDF(int peak, int iy, int ipt);
double PDF_shape(double *x, double *par);
double shape_m(double *x, double *par);
int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb);
