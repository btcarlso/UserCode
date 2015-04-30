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

bool MakeEfficiency=true;
bool test=false; 

Long64_t entries=-9; 

Float_t st=-9; 
Float_t ht=0; 
Float_t deltaPt=0;
Float_t gen_st=0; 
Float_t gen_ht=0;

int eventNo;
int runNo;
int lumiNo; 
int NumInteractions;

Float_t puWt_nom;
Float_t puWt_up;
Float_t puWt_down;

double NG=0; //Ngen
double filter_eff=1;
double Nacc=0;
double weight_noXS=1;

float jec_sys=0;

Int_t nJets=-9;

int nLooseBJets=-1;
int nTightBJets=-1;
int nBJets=-1;

Int_t nb_truth=0;
Int_t nl_truth=0;
Int_t nc_truth=0;
Int_t n_split=0;

Int_t nPv=-9;
int nMuons;
int nPhotons;
int nElectrons;

int nLooseElectrons;
int nLooseMuons; 


Float_t met=0;
Float_t met_phi=0;
Float_t mt=0; 
float muon_pT=-9; 
TLorentzVector wP; 
TLorentzVector ttP;
TLorentzVector zP; 
TLorentzVector jetPsys;
TLorentzVector jetPttsys;

double st_top=0;


bool SS=false; 
bool dimuon_enriched=false;
bool Z_enriched_on=false;
bool Z_enriched_off=false;

bool tt_enriched=false;


float dimuon_mass=0; 
float dilepton_mass=0;
float dilepton_pt=0; 

float pTmuon1=0;
float pTmuon2=0;
float pTElectron1=0;

float pTbjet1=0;
float pTbjet2=0;

float weight=0; 
double B_weight;
double B_weightP;//sys+
double B_weightM;//sys-

double B_weightPbc;//sys+
double B_weightMbc;//sys-

double B_weightPl;//sys+
double B_weightMl;//sys-

double B_weight_tight;
double B_weight_tightP;//sys+
double B_weight_tightM;//sys-

double B_weight_tightPbc;//sys+
double B_weight_tightMbc;//sys-

double B_weight_tightPl;//sys+
double B_weight_tightMl;//sys-

double probMC0;
double probMC1;
double probMC2;
double probMC0_2;


double probMC0_SF;
double probMC1_SF;
double probMC2_SF;
double probMC0_2_SF;

double probMC0_SFP;
double probMC1_SFP;
double probMC2_SFP;
double probMC0_2_SFP;

double probMC0_SFM;
double probMC1_SFM;
double probMC2_SFM;
double probMC0_2_SFM;

bool QCD=false;


int nJetmin=1; 
int nJetmax=7;
int nBmax=2;
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
std::map<TString, TProfile2D*> prName2D;
std::map<TString, TEfficiency*> effName;
std::map<TString, TStopwatch> stopwatch;
std::map<TString,double> event_weight;

bool _MC=false;
bool _fastSim=false;
vector<TString> sample_list; 
TString output_file_name; 
TString selections="";

TFile *output_file;
//=new TFile("output_file_trigger_wJets_inclusive_1mu_xs.root","RECREATE"); 


TFile *fEff = new TFile("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root","READ");
TFile *fEff_elec=new TFile("CutBasedIdScaleFactors.root","READ");
TFile *fBtagEff;
TFile *fCorr = new TFile("corrections.root","READ"); 


TProfile *pr24 = new TProfile("efficiency_24","#epsilon_#mu;s_{T}",300,0,3000); 
TProfile *pr40 = new TProfile("efficiency_40","#epsilon_#mu;s_{T}",300,0,3000); 

TF1 *expo_corr;
double TTbar_corr; 
double TTbar_corrP;
double TTbar_corrM;

double sumW_EMu=0;
double sum_EMu=0;

bool fill_histos=true;
float calcSt();
string jet_label(int flv);
void fill_acceptance(int iBin);
void acceptance();
void setup_files(int jobnumber);
float event_rapidity();
int list_files(TString dirname, TString ext);
void analyze_trees(int jobnumber);
void initialize_tree(TChain *tree);
void initialize();
void analyze_file(TString folder);
void control_regions();
double getSF(int ijet, TString Atagger, TString sys);
double SF_light(TString Atagger, TString mode, TLorentzVector jetP);
void jetcomb_mass();
void kinematics();
bool basic_selections();
void event_loop(TChain *tree);

double read_btag_efficiency(TLorentzVector pJ, int flavor,string Tagger);
void fill_btag_efficiency(TLorentzVector pJ, int flavor, bool btagL, bool btagM, bool btagT);
void print_jets();
void print_muons();
void remove_muon(int i);
void remove_jet(int ijet);

void fill_inclusive(); 
void fill_discrete();
void pileup_reweight();
void gen_pt();

void fill_pileup();
void fill_categories();
void EventWeight();
void fill_ref_plots();
void fill_event_category();
//void check_jet_cleaning();
bool cuts();


//loops over objects, eg. Muons,electrons,photons jets

void quantile_ratio(TH1F *h_num, TH1F *h_den, TH1F *h_num_uW, TH1F *h_den_uW,TGraphAsymmErrors *g_rat);
void bookEff();
TString create_category(int nJ, bool inclusive, int Nmu, int Nel, int charge, int Nb, TString tag, TString sf, TString mass);
TString histName(TString variable, TString selection, TString btag_sel,TString scales);
void bookHistoRef();
void bookHisto();
void open_files();

void efficiency();
void nBtag();

double get_eff(double pt, double eta, TString mode);
void open_graph(TString name);
void get_muTrigEff();
void load_btag_sys();
void load_btagEff();
void load_FullFastSimSF();



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
void CreateProfile(const char* name,   const char* title,
                   const char* xTitle, const char* yTitle,
                   Int_t       nBinsX, Double_t    xLow, Double_t xUp,
                   Int_t   nBinsY, Double_t    yLow, Double_t yUp,Double_t zLow, Double_t zUp   );

void CreateProfile(const char* name,   const char* title,
				   const char* xTitle, const char* yTitle,
				   Int_t       nBinsX, Float_t *xBins);
void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, const Double_t* xBins,
					 Int_t nBinsY,Double_t yLow, Double_t yUp);
void CreateEfficiency(TString name, TString title);
void writeHisto();
//Muons

//tight muons
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

//loose muons

vector<float> *loose_muon_px=0;
TBranch *b_loose_muon_px;

vector<float> *loose_muon_py=0;
TBranch *b_loose_muon_py;

vector<float> *loose_muon_pz=0;
TBranch *b_loose_muon_pz;

vector<float> *loose_muon_e=0;
TBranch *b_loose_muon_e;

vector<float> *loose_muon_charge=0;
TBranch *b_loose_muon_charge;


//Tight Electrons
vector<float> *electron_px=0;
TBranch *b_electron_px;

vector<float> *electron_py=0;
TBranch *b_electron_py;

vector<float> *electron_pz=0;
TBranch *b_electron_pz;

vector<float> *electron_e=0;
TBranch *b_electron_e;

vector<float> *electron_charge=0;
TBranch *b_electron_charge;

//Loose Electrons

vector<float> *loose_electron_px=0;
TBranch *b_loose_electron_px;

vector<float> *loose_electron_py=0;
TBranch *b_loose_electron_py;

vector<float> *loose_electron_pz=0;
TBranch *b_loose_electron_pz;

vector<float> *loose_electron_e=0;
TBranch *b_loose_electron_e;

vector<float> *loose_electron_charge=0;
TBranch *b_loose_electron_charge;

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

vector<bool> *jet_unc=0;
TBranch *b_jet_unc;

vector<bool> *jet_bTagL=0;
TBranch *b_jet_bTagL;

vector<bool> *jet_bTagM=0;
TBranch *b_jet_bTagM;

vector<bool> *jet_bTagT=0;
TBranch *b_jet_bTagT;

vector<int> *jet_algFlavor=0;
TBranch *b_jet_algFlavor;

vector<int> *jet_phyFlavor=0;
TBranch *b_jet_phyFlavor;


//W
vector<float> *w_px=0;
TBranch *b_w_px;

vector<float> *w_py=0;
TBranch *b_w_py;

vector<float> *w_pz=0;
TBranch *b_w_pz;

vector<float> *w_e=0;
TBranch *b_w_e;	

//Z
vector<float> *z_px=0;
TBranch *b_z_px;

vector<float> *z_py=0;
TBranch *b_z_py;

vector<float> *z_pz=0;
TBranch *b_z_pz;

vector<float> *z_e=0;
TBranch *b_z_e;	

//top 
vector<float> *top_px=0;
TBranch *b_top_px;

vector<float> *top_py=0;
TBranch *b_top_py;

vector<float> *top_pz=0;
TBranch *b_top_pz;

vector<float> *top_e=0;
TBranch *b_top_e;	



