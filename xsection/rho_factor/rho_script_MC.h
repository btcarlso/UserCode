/*
 *  rho_script.h
 *  
 *
 *  Created by Benjamin Carlson on 5/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "TObject.h"
#include "TEnv.h"
#include <TROOT.h>
#include <TSystem.h>
#include "TClass.h"
#include "TStopwatch.h"

//formatting 
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include <TFitResult.h>
#include <THStack.h>
#include <TLegend.h>
#include <TAttMarker.h>
#include <TProfile.h>
#include <TChain.h>
#include <TGaxis.h> 
#include <TAxis.h>
#include <TEfficiency.h>
#include <TMarker.h>
#include <TPave.h>
#include <TPaveStats.h>

//random stuff
#include <TText.h>
#include <TError.h>

//io stuff
#include <TTree.h>
#include <TDirectory.h>
#include "TString.h"
#include "TFile.h"

//Math stuff
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TRandom3.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TMinuit.h>
#include "TVirtualFitter.h"

bool trigger_mumu; 
bool trigger_mu;
Double_t        fP_fX;
Double_t        fP_fY;
Double_t        fP_fZ;
Double_t        fE;

TLorentzVector *Ups_P; 
TLorentzVector *Mu1_P; 
TLorentzVector *Mu2_P;

TLorentzVector *Ups_P_reco; 
TLorentzVector *Mu1_P_reco; 
TLorentzVector *Mu2_P_reco;

int muPos_qual =-1; 
int muNeg_qual =-1; 
int NCand=-9;



Int_t D5_v1=0, D5_v2=0, D5_v3=0, D5_v5=0, D7_v1=0, D7_v4=0, D9_v1=0, D9_v4=0;
Int_t D0_v1=0, D0_v2=0, D0_v3=0, D0_v5=0, D0_v6=0, D0_v9=0;
Int_t DMu3_v1=0;
Int_t DUM0_v1=0, DUM0_v2=0, DUM0_v3=0, DUM0_v4=0, DUM0_v6=0;


Double_t vProb=0; 
Double_t distM1=0; 

Int_t HLT_Mu5_L2Mu2_v1; 
Int_t HLT_Dimuon10_Jpsi_Barrel_v6; 

Int_t runNb; 
float m=0, pt=0, ptMuP=0, ptMuM=0, etaMuP=0, etaMuM=0;
float y=0;
Int_t SeaGull=-1; 
Int_t Trig5=-1; 
Int_t Trig7=-1; 
Int_t Trig9=-1;


TFile *fEff = new TFile("MCEff2.root","READ");

TFile *fOut = new TFile("output_tree.root","RECREATE");

void set_output_tree();
void fill_tree();
void rho_script_MC();
void loop();
void set_branches(TTree *tree);
bool select_jpsi(); 
bool select_upsilons();
void chain();
int eta_to_bin(double eta);
double muonEff(double pt,double eta,string sys);
TTree *fTree; 
TChain *tree; 

