#include "/uscms/home/btcarlso/stealth_code/root_headers.h"
#include "Math/QuantFuncMathCore.h"
#include "TVirtualFitter.h"

void analyze_trees(int jobnumber);
void setup_files(int jobnumber);
void open_files();
void prepareNormalization();
void analyze_file(TString fName);
void initialize_tree(TChain *tree);
void event_loop(TChain *tree);
void bookHisto();
void writeHisto();

void fill_LQ();

std::vector<TString> sample_list;
bool isMC=false;
bool isMuEle=false;
//bool singleEle=false;
TString output_file_name;// sample_name;

bool _fastSim=false;
TString selections="";
bool QCD=false;
bool test=false;
double filter_eff=1;

TFile *output_file;
TFile *fBtagEff;

std::map<TString, TString> fileName;
std::map<TString, float> xs;

// btagging SF related stuff
std::map<TString, TEfficiency*> effName;
double read_btag_efficiency(TLorentzVector pJ, int flavor, string Tagger);
double getSF(TLorentzVector pJi, int algFlv, TString Atagger, TString sys);
double SF_light(TString Atagger, TString mode, TLorentzVector jetP);
void load_btag_sys();

// histogram booking

TString dir="/eos/uscms/store/user/btcarlso/trees/RPV/";

//global variables
Float_t mT(0.), hT(0.), sT(0.), dilepMass(0.);
float lepPt(0);

float dimuon_mass(0);
float dielectron_mass(0);
Int_t nBTags(-1);

int NG(0);

bool count_muons();
bool count_electrons();
void fill_jets();
void loop_leptons();

//cuts:
float weight(0.), lumiWeight(0.), lumi(19690.);
int nJetmin(1),nJetmax(8);
float m_muj(0);
std::vector<TLorentzVector> jetCollection;
std::vector<TLorentzVector> electronCollection;
std::vector<TLorentzVector> muonCollection;

float leptPt(0);

float B_weight(1);

float jetPt1_cut=125;
float jetPt2_cut=45;
float jetPt_cut=45;

float stMin(300);

float min_diLeptonMass=50;
float ZVeto(10);
float Zmass(91);

bool dimuon_enriched(0);
bool dielectron_enriched(0);
bool tt_enriched(0);
bool Z_enriched(0);
bool SS(0);

void print_event();

//------------
//branches
//-------------
Long64_t entries=-9;
Int_t nJets(-9), nPv(-9), nMuons, nElectrons;
int nLooseElectrons(0), nLooseMuons(0);
Float_t puWt_nom(1.), st(-9.), met(0.), met_phi(0.);

//Electrons
vector<float> *electron_px=0;
TBranch *b_electron_px;

vector<float> *electron_py=0;
TBranch *b_electron_py;

vector<float> *electron_pz=0;
TBranch *b_electron_pz;

vector<float> *electron_e=0;
TBranch *b_electron_e;

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

//vector<float> *electron_eta=0;
//TBranch *b_electron_eta;

vector<float> *electron_charge=0;
TBranch *b_electron_charge;
/*
vector<float> *electron_iso=0;
TBranch *b_electron_iso;

vector<float> *electron_mva=0;
TBranch *b_electron_mva;
*/
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

vector<bool> *jet_btagL=0;
TBranch *b_jet_btagL;

vector<bool> *jet_bTagM=0;
TBranch *b_jet_bTagM;

vector<bool> *jet_btagT=0;
TBranch *b_jet_btagT;

vector<int> *jet_algFlavor;
TBranch *b_jet_algFlavor;

vector<int> *jet_phyFlavor;
TBranch *b_jet_phyFlavor;

vector<float> *jet_unc=0;
TBranch *b_jet_unc;

// methods to find nBTags, Mt, Ht, St

void calcSt();
void calcMuJ();
void  calDilepMass();

Float_t  calcEleIdSF(Double_t eta, Double_t pt);
Float_t  calcEleTrigSF(Double_t eta, Double_t pt);
Float_t  electronSF(std::vector<TLorentzVector> eleP4);
Float_t  calcMuIdSF(Double_t eta, Double_t pt);
Float_t  muonSF(std::vector<TLorentzVector> muP4);
//TLorentzVector p4Lep(Double_t px,Double_t py, Double_t pz, Double_t e);

//histogram booking
std::map<TString, TH1F*> hName;
std::map<TString, TH2F*> hName2D;
std::map<TString, TProfile*> prName;

void load_btagEff();
bool preselection_selection();

//void BtagSF();

void createHistogram(const char* name,   const char* title, 
                     const char* xTitle, const char* yTitle,
                     Int_t       nBinsX, Double_t    xLow, Double_t xUp);
void createHistogram(const char* name,   const char* title,
                     const char* xTitle, const char* yTitle,
                     Int_t       nBinsX, const Float_t* xBins);
void createHistogram(const char* name,   const char* title,
                     const char* xTitle, const char* yTitle,
                     Int_t nBinsX, Double_t xLow, Double_t xUp,
                     Int_t nBinsY,Double_t yLow, Double_t yUp);
void createProfile(const char* name,     const char* title,
                   const char* xTitle,   const char* yTitle,
                   Int_t       nBinsX,   Double_t    xLow, Double_t xUp);


//std::map<TString, float> SF;
//TFile *output_file;
