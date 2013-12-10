/*
 *  analyze_trees.h
 *  
 *
 *  Created by Benjamin Carlson on 6/28/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "/uscms/home/btcarlso/stealth_code/root_headers.h"


#include "Math/QuantFuncMathCore.h"
#include "TVirtualFitter.h"


Long64_t entries=-9; 

Float_t st=-9; 
Int_t nJets=-9; 
Int_t nBJets=-1; 
Int_t nPv=-9; 
int nMuons;
int nPhotons;
int nElectrons;
Float_t met=0;
Float_t met_phi=0;
Float_t mt=0; 
float muon_pT=-9; 

float weight=0; 

bool QCD=false; 

int nJetmin=1; 
int nJetmax=7;
double st_min=800;
double st_max=3000;

bool LO=false;

std::map<TString, TH1F*> hName;
std::map<TString, TH2F*> hName2D;
std::map<TString, TGraphAsymmErrors*> grName;
std::map<TString, float> xs; 
std::map<TString, float> Ngen; 
std::map<TString, TString> fileName; 
std::map<TString, TProfile*> prName; 

bool _MC=false; 
vector<TString> sample_list; 
TString output_file_name; 
TString selections="";

TFile *output_file;
//=new TFile("output_file_trigger_wJets_inclusive_1mu_xs.root","RECREATE"); 


TFile *fEff = new TFile("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root","READ"); 


TProfile *pr24 = new TProfile("efficiency_24","#epsilon_#mu;s_{T}",300,0,3000); 
TProfile *pr40 = new TProfile("efficiency_40","#epsilon_#mu;s_{T}",300,0,3000); 

bool fill_histos=true;
float calcHt();
float calcSt();
void setup_files(int jobnumber);
float event_rapidity();
void analyze_trees(int jobnumber);
void initialize_tree(TChain *tree);
void initialize();
void analyze_file(TString folder);
void event_loop(TChain *tree);

void print_jets();
void print_muons();
void remove_muon(int i);
void remove_jet(int ijet);


void fill_st();
void fill_jets(); 
void fill_muon(); 
void fill_rapidity();
void fill_deltaR();
void fill_bjet_mu_mass();
bool cuts();


//loops over objects, eg. Muons,electrons,photons jets

void divide_BW();
void quantile_ratio(TH1F *h_num, TH1F *h_den, TH1F *h_num_uW, TH1F *h_den_uW,TGraphAsymmErrors *g_rat);
void bookHisto();
void open_files();

void efficiency();
void nBtag();

double get_eff(double pt, double eta, TString mode);
void open_graph(TString name);
void get_muTrigEff();



void CreateHistogram(const char* name,   const char* title, const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, Double_t xLow, Double_t xUp,
					 Int_t nBinsY,Double_t yLow, Double_t yUp);
void CreateProfile(const char* name,   const char* title,
				   const char* xTitle, const char* yTitle,
				   Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void writeHisto();
//Muons


vector<float> *muon_px=0;
TBranch *b_muon_px;

vector<float> *muon_py=0;
TBranch *b_muon_py;

vector<float> *muon_pz=0;
TBranch *b_muon_pz;

vector<float> *muon_e=0;
TBranch *b_muon_e;

vector<float> *muon_charge=0;
TBranch *b_muon_charge;

//Electrons
vector<float> *electron_px=0;
TBranch *b_electron_px;

vector<float> *electron_py=0;
TBranch *b_electron_py;

vector<float> *electron_pz=0;
TBranch *b_electron_pz;

vector<float> *electron_e=0;
TBranch *b_electron_e;
//Photons
vector<float> *photon_px=0;
TBranch *b_photon_px;

vector<float> *photon_py=0;
TBranch *b_photon_py;

vector<float> *photon_pz=0;
TBranch *b_photon_pz;

vector<float> *photon_e=0;
TBranch *b_photon_e;
//Jets
vector<float> *jet_px=0;
TBranch *b_jet_px;

vector<float> *jet_py=0;
TBranch *b_jet_py;

vector<float> *jet_pz=0;
TBranch *b_jet_pz;

vector<float> *jet_e=0;
TBranch *b_jet_e;	

vector<bool> *jet_btag=0;
TBranch *b_jet_btag;

