// -*- C++ -*-
//
// Package:    JPsiAnalyzerPAT
// Class:      JPsiAnalyzerPAT
// 
/**\class JPsiAnalyzerPAT JPsiAnalyzerPAT.cc OctoberXTracking/JPsiAnalyzerPAT/src/JPsiAnalyzerPAT.cc
 
 Description: <one line class summary>
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author: Roberto Covarelli 
//         Created:  Fri Oct  9 04:59:40 PDT 2009
//$Id: JPsiAnalyzerPAT.cc,v 1.53.2.3 2012/05/02 20:19:51 eaguiloc Exp $
//
// based on: Onia2MuMu package V00-11-00
// changes done by: FT-HW

// system include files
#include <memory>
#include <fstream>
#include <ostream>
#include <iostream>
#include <math.h>

// ROOT/Roofit include files
#include <TStyle.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TVector3.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/BeamSpot/interface/BeamSpot.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h>

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//includes for MomentumScaleCalibration:
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/MomentumScaleCorrector.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/ResolutionFunction.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/BackgroundFunction.h"
#include "CondFormats/RecoMuonObjects/interface/MuScleFitDBobject.h"
#include "CondFormats/DataRecord/interface/MuScleFitDBobjectRcd.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/MuScleFitUtils.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace RooFit;

//
// class declaration
//
class JPsiAnalyzerPAT : public edm::EDAnalyzer {
public:
	explicit JPsiAnalyzerPAT(const edm::ParameterSet&);
	~JPsiAnalyzerPAT();
	
	
private:
	void bookHisto();
	double muonEff(double pt,double eta,string sys); 
	void CreateHistogram(const char* name,   const char* title,
						 const char* xTitle, const char* yTitle,
						 Int_t       nBinsX, Double_t    xLow, Double_t xUp);
	void CreateHistogram(const char* name,   const char* title,
						 const char* xTitle, const char* yTitle,
						 Int_t       nBinsX, const Float_t* xBins);
	
	void CreateHistogram2D(const char* name,   const char* title,
						   const char* xTitle, const char* yTitle,
						   Int_t nBinsX, Double_t xLow, Double_t xUp,
						   Int_t nBinsY,Double_t yLow, Double_t yUp);
	void CreateProfile(const char* name,   const char* title,
				  const char* xTitle, const char* yTitle,
				  Int_t       nBinsX, Double_t    xLow, Double_t xUp);
	
	int eta_to_bin(double eta);
	
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	void makeCuts() ;
	int theBestQQ();
	double M_uncertainty(TVector3 fMuon1Vect, TVector3 fMuon2Vect, double *fMuon1E, double *fMuon2E); 
	double correction_functionJPsi(const double & pt, const double & eta, const double & phi, const int chg, double *par);
	double correction_functionZ(const double & pt, const double & eta, const double & phi, const int chg,   double *parScale);
	double Pt_primeZ(TVector3 v1,int Q);
	double Pt_primeJPsi(TVector3 v1,int Q); 
	
	
	//Acceptance stuff
	
	void error_graph(TGraph *gr, TString name);
	void OpenPolarizationFiles(); 
	void HX_boost(); 
	void FillAcceptance(double w, bool pass, string method); 
	void acceptance(); 
	
	void fillHistograms(const pat::CompositeCandidate* aCand, const edm::Event& iEvent);
	void fillTreeAndDS(const pat::CompositeCandidate* aCand, const edm::Event&);
	bool isMuonInAccept(const pat::Muon* aMuon);
	bool selMuon(const pat::Muon* aMuon);
	bool selDimuon(const pat::CompositeCandidate* aCand);
	int getJpsiVarType(const double jpsivar, vector<double> vectbin);
	double CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode, TLorentzVector &DiMuon, TVector3 &Mu1, TVector3 &Mu2);
	
	// additional functions by f
	void resetDSVariables();
	void analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles);
	// void calcPol(TLorentzVector&, TLorentzVector&, std::vector< float >&, std::vector< float >& );
	void beginRun(const edm::Run &, const edm::EventSetup &);
	void hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup);
	void matchMuonToHlt(const pat::Muon*, const pat::Muon*);
	void muonStationDistance (const pat::CompositeCandidate* aCand);
	
	//variable for scale corrections:
	std::auto_ptr<MomentumScaleCorrector> corrector_;
	std::auto_ptr<ResolutionFunction> resolutionFunction_;
	
	
	// ROOT tree 
	TTree* tree_;//data; //*recoData;
	TFile* fOut_;
	TFile* fEff;
	TFile *fRho; 
	
	//histogram map
	std::map<TString, TH1F*> hName;
	std::map<TString, TH2F*> hName2D;
	std::map<TString, TGraph*> polName;
	std::map<TString, TProfile*> prName; 
	
	// SMALL dataset and RooRealVars
	TFile* fOut2_;
	RooDataSet* data;
	RooRealVar* Jpsi_MuScleMass;      
	RooRealVar* Jpsi_MuScleMassErr;
	RooRealVar* Jpsi_MassErr;
	RooRealVar* Jpsi_Pt;
	RooRealVar* Jpsi_Rap; 
	RooRealVar* Jpsi_ct;
	RooRealVar* Jpsi_ctErr;
	RooRealVar* Jpsi_ctTrue;      			
	RooCategory* Jpsi_PtType;
	RooCategory* Jpsi_RapType;
	RooCategory* Jpsi_MatchType;
	
	//1.) J/psi variables RECO
	// double JpsiMass, JpsiPt, JpsiRap;
	// double JpsiPx, JpsiPy, JpsiPz;
	TLorentzVector* JpsiP;
	double Jpsict, JpsictErr, JpsiVprob;
	int MCType; //GG, GT and TT
	//       int JpsiType,  JpsiCharge, MCType; //GG, GT and TT
	double JpsiMuScleMassCorr, JpsiMuScleMassErr, JpsiMassErr;
	double sigmaPtPos, sigmaPtNeg;
	double JpsiDistM1, JpsiDphiM1, JpsiDrM1;
	double JpsiDistM2, JpsiDphiM2, JpsiDrM2;
	
	
	double cosTheta, sinTheta,cosPhi;//HX frame variables
	
	double dca;
	//2.) muon variables RECO
	// double muPosPx, muPosPy, muPosPz;
	TLorentzVector* muPosP;
	// double muNegPx, muNegPy, muNegPz;
	TLorentzVector* muNegP;
	
	//3.) J/psi variables GEN
	// double JpsiMass_Gen, JpsiPt_Gen, JpsiRap_Gen;
	// double JpsiPx_Gen,   JpsiPy_Gen, JpsiPz_Gen;
	TLorentzVector* JpsiP_Gen;
	double Jpsict_Gen;
	
	//4.)muon variables GEN
	// double muPosPx_Gen, muPosPy_Gen, muPosPz_Gen;
	TLorentzVector* muPosP_Gen;
	// double muNegPx_Gen, muNegPy_Gen, muNegPz_Gen;
	TLorentzVector* muNegP_Gen;
	
	//5.) Event related variables
	unsigned int eventNb, runNb, lumiBlock, nPriVtx;
	int countTksOfPV;
	double vertexWeight, sumPTPV;
	int muPos_qual, muNeg_qual;
	int QQ_check;
	int NCand;
	int TriggerBit;
	
	//6.) POL variables
	// std::vector<std::string> polVarNames_;
	// std::vector<std::string> polVarNamesGen_;
	// std::map<std::string, double> mapPolVarsToValue_;
	// std::map<std::string, double> mapPolVarsToValueGen_;
	
	//7.) TriggerNames Map
	std::map<std::string, int> mapTriggerNameToIntFired_;
	// USAGE:
	// ---> For all single or symmetric double triggers:
	// 0 : event not firing the corresponding trigger
	// 1 : event firing the corresponding trigger and all RECO muons matched to the required number of HLT objects (1 or 2)
	// 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects
	// 
	// ---> For the asymmetric triggers:
	// 0 : event not firing the corresponding trigger
	// 1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, POSITIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
	// -1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, NEGATIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
	// 2 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, both matched to the tighter HLT object (usually a L3 muon)
	// 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects 
	std::map<std::string, int> mapTriggerNameToPrescaleFac_;
	
	Handle<pat::CompositeCandidateCollection > collAll;
	// Handle<TriggerResults> trigger;
	
	// data members
	InputTag       _patJpsi;
	bool           _writeTree;
	string         _treefilename; 
	bool           _writeDataSet; 
	string         _datasetname;
	vector<string> _triggerForDataset;
	double         _massMin;
	double         _massMax;
	vector<double> _ptbinranges;
	vector<double> _etabinranges;
	bool           _onlythebest;
	bool           _applycuts;
	bool           _applyExpHitcuts;
	bool           _applyDiMuoncuts;
	bool           _useBS;
	bool           _removeSignal;
	bool           _removeMuons;
	bool           _writeOutCands;
	int            _MassCorr;
	// bool           _JSON;
	int            _oniaPDG;
	InputTag       _genParticles;
	bool           _isMC;
	bool           _storeAllMCEvents;
	bool           _isPromptMC;
	
	InputTag      _triggerresults;  
	
	vector<const pat::CompositeCandidate*>   _thePassedCands;
	
	// number of events
	unsigned int nEvents;
	unsigned int passedTriggerResults_;
	unsigned int passedMuonSelectionCuts_;
	unsigned int passedTriggerMatch_;
	unsigned int passedTriggerResultsAnalyzer_;
	
	// limits 
	float JpsiMassMin;
	float JpsiMassMax;
	float JpsiCtMin;
	float JpsiCtMax;
	float JpsiPtMin;           // SET BY 
	float JpsiPtMax;           // DEFINITION
	float JpsiRapMin;          // OF BIN
	float JpsiRapMax;          // LIMITS
	
	math::XYZPoint RefVtx;
	ofstream* theTextFile;
	ofstream* JSON;
	
	int runtmp,lumitmp,count;
	int runmax,runmin;
	
	//muon distance
	PropagateToMuon prop1_, prop2_; 
	
	// Trigger Filter Studies
	edm::Handle< edm::TriggerResults> handleTriggerResults_;
	edm::InputTag tagTriggerResults_;
	bool          requireTriggerMatching_;
	HLTConfigProvider hltConfig_;
	bool hltConfigInit_;
	std::vector<std::string> HLTBitNames_;
	std::vector<std::string> HLTBitNames_SingleMu;
	std::vector<std::string> HLTLastFilterNames_SingleMu;
	std::vector<std::string> HLTBitNames_DoubleMu;
	std::vector<std::string> HLTLastFilterNames_DoubleMu;
	std::vector<std::string> HLTBitNames_MuL2Mu;
	std::vector<std::string> HLTLastFilterNames_MuL2Mu;
	std::vector<std::string> HLTBitNames_MuTrack;
	std::vector<std::string> HLTLastFilterNames_MuTrack;
	std::vector<std::string> HLTBitNames_MuTkMu;
	std::vector<std::string> HLTLastFilterNames_MuTkMu;
	std::map< std::string, std::string> mapTriggerToLastFilter_;
	
};    

// constants, enums and typedefs
enum {CS, HX, PHX, sGJ, GJ1, GJ2};
//

//
// constructors and destructor
//
JPsiAnalyzerPAT::JPsiAnalyzerPAT(const edm::ParameterSet& iConfig):
_patJpsi(iConfig.getParameter<InputTag>("src")),
_writeTree(iConfig.getParameter<bool>("writeTree")),
_treefilename(iConfig.getParameter<string>("treeFileName")),	
_writeDataSet(iConfig.getParameter<bool>("writeDataSet")),
_datasetname(iConfig.getParameter<string>("dataSetName")),
_triggerForDataset(iConfig.getParameter< vector<string> >("triggersForDataset")),
_massMin(iConfig.getParameter<double>("massMin")),
_massMax(iConfig.getParameter<double>("massMax")),
_ptbinranges(iConfig.getParameter< vector<double> >("pTBinRanges")),	
_etabinranges(iConfig.getParameter< vector<double> >("etaBinRanges")),
_onlythebest(iConfig.getParameter<bool>("onlyTheBest")),		
_applycuts(iConfig.getParameter<bool>("applyCuts")),
_applyExpHitcuts(iConfig.getUntrackedParameter<bool>("applyExpHitCuts",false)),
_applyDiMuoncuts(iConfig.getUntrackedParameter<bool>("applyDiMuonCuts",false)),
_useBS(iConfig.getParameter<bool>("useBeamSpot")),
_removeSignal(iConfig.getUntrackedParameter<bool>("removeSignalEvents",false)),
_removeMuons(iConfig.getUntrackedParameter<bool>("removeTrueMuons",false)),
_writeOutCands(iConfig.getUntrackedParameter<bool>("writeOutCandidates",false)),
_MassCorr(iConfig.getParameter<int>("massCorrectionMode")),
_oniaPDG(iConfig.getParameter<int>("oniaPDG")),
_genParticles(iConfig.getParameter<InputTag>("genParticles")),
_isMC(iConfig.getUntrackedParameter<bool>("isMC",false)),
_storeAllMCEvents(iConfig.getUntrackedParameter<bool>("storeAllMCEvents",false)),
_isPromptMC(iConfig.getUntrackedParameter<bool>("isPromptMC",false) ),
prop1_(iConfig.getParameter<edm::ParameterSet>("propagatorStation1")),
prop2_(iConfig.getParameter<edm::ParameterSet>("propagatorStation2")),
tagTriggerResults_(iConfig.getParameter<InputTag>("triggerResultsLabel")),
requireTriggerMatching_(iConfig.getUntrackedParameter<bool>("requireTriggerMatching",true)),
HLTBitNames_SingleMu(iConfig.getParameter< vector<string> >("HLTBitNames_SingleMu")),
HLTLastFilterNames_SingleMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_SingleMu")),
HLTBitNames_DoubleMu(iConfig.getParameter< vector<string> >("HLTBitNames_DoubleMu")),
HLTLastFilterNames_DoubleMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_DoubleMu")),
HLTBitNames_MuL2Mu(iConfig.getParameter< vector<string> >("HLTBitNames_MuL2Mu")),
HLTLastFilterNames_MuL2Mu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_MuL2Mu")),
HLTBitNames_MuTrack(iConfig.getParameter< vector<string> >("HLTBitNames_MuTrack")),
HLTLastFilterNames_MuTrack(iConfig.getParameter< vector<string> >("HLTLastFilterNames_MuTrack")),
HLTBitNames_MuTkMu(iConfig.getParameter< vector<string> >("HLTBitNames_MuTkMu")),
HLTLastFilterNames_MuTkMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_MuTkMu"))
{
	//now do what ever initialization is needed
	nEvents = 0; 
	// passedTriggerResults_=0;
	passedMuonSelectionCuts_=0;
	passedTriggerMatch_=0;
	// passedTriggerResultsAnalyzer_=0;
	
	JpsiMassMin = _massMin;
	JpsiMassMax = _massMax;
	JpsiCtMin = -2.0;
	JpsiCtMax = 3.5;
	
	if (_writeOutCands) theTextFile = new ofstream("passedCandidates.txt");
	
	if (HLTBitNames_SingleMu.size() != HLTLastFilterNames_SingleMu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
	std::vector<std::string>::iterator it2 = HLTLastFilterNames_SingleMu.begin();
	for(std::vector<std::string>::iterator it = HLTBitNames_SingleMu.begin(); it != HLTBitNames_SingleMu.end(); ++it){
		mapTriggerToLastFilter_[*it] = *it2;
		++it2;   HLTBitNames_.push_back(*it);
	}
	if (HLTBitNames_DoubleMu.size() != HLTLastFilterNames_DoubleMu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
	it2 = HLTLastFilterNames_DoubleMu.begin();
	for(std::vector<std::string>::iterator it = HLTBitNames_DoubleMu.begin(); it != HLTBitNames_DoubleMu.end(); ++it){
		mapTriggerToLastFilter_[*it] = *it2;
		++it2;   HLTBitNames_.push_back(*it);
	}
	if (2*HLTBitNames_MuL2Mu.size() != HLTLastFilterNames_MuL2Mu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
	it2 = HLTLastFilterNames_MuL2Mu.begin();
	for(std::vector<std::string>::iterator it = HLTBitNames_MuL2Mu.begin(); it != HLTBitNames_MuL2Mu.end(); ++it){
		std::string theName = *it;
		mapTriggerToLastFilter_[theName] = *it2;
		++it2;   HLTBitNames_.push_back(theName);  
		theName += "_special";
		mapTriggerToLastFilter_[theName] = *it2;
		++it2;    
	}
	if (HLTBitNames_MuTrack.size() != HLTLastFilterNames_MuTrack.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
	it2 = HLTLastFilterNames_MuTrack.begin();
	for(std::vector<std::string>::iterator it = HLTBitNames_MuTrack.begin(); it != HLTBitNames_MuTrack.end(); ++it){
		mapTriggerToLastFilter_[*it] = *it2;
		++it2;   HLTBitNames_.push_back(*it);
	}
	if (HLTBitNames_MuTkMu.size() != HLTLastFilterNames_MuTkMu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
	it2 = HLTLastFilterNames_MuTkMu.begin();
	for(std::vector<std::string>::iterator it = HLTBitNames_MuTkMu.begin(); it != HLTBitNames_MuTkMu.end(); ++it){
		mapTriggerToLastFilter_[*it] = *it2;
		++it2;  HLTBitNames_.push_back(*it);
	}
	
	for(std::vector<std::string>::iterator it = HLTBitNames_.begin(); it != HLTBitNames_.end(); ++it){
		mapTriggerNameToIntFired_[*it] = -9999;
	}
}


JPsiAnalyzerPAT::~JPsiAnalyzerPAT()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	if (_writeOutCands) theTextFile->close();
}

// ------------ method called once each job just before starting event loop  ------------

void JPsiAnalyzerPAT::bookHisto(){
	
	TH1F *h=(TH1F*)fRho->Get("rho_dRPtE_TH1"); 
	hName["rho_dRPtE_TH1"]=h; 
	
	CreateHistogram("rho_pt_num","#rho ","p_{T} [GeV]","#rho",50,10,100);
	CreateHistogram("rho_pt_num_uW","#rho ","p_{T} [GeV]","#rho",50,10,100);

	CreateHistogram("rho_pt_den","#rho Denominator","p_{T} [GeV]","#rho",50,10,100);
	
	CreateHistogram("rho_dR_num","#rho ","#DeltaR","#rho",100,0,3);
	CreateHistogram("rho_dR_den","#rho Denominator","#DeltaR","#rho",100,0,3);
	
	CreateHistogram("rho_DM1_num","#rho ","dist M1 [cm]","#rho",150,0,1500); 
	CreateHistogram("rho_DM1_den","#rho Denominator","dist M1 [cm]","#rho",150,0,1500); 

	CreateHistogram("rho_dRPtE_num","#rho ","#DeltaR_{#Deltap_{T}}^{elliptic}","#rho",100,0,3);
	CreateHistogram("rho_dRPtE_num_uW","#rho unweighted","#DeltaR_{#Deltap_{T}}^{elliptic}","#rho",100,0,3); //unweighted
	CreateHistogram("rho_dRPtE_den","#rho Denominator","#DeltaR_{#Deltap_{T}}^{elliptic}","#rho",100,0,3); 
	
	CreateHistogram("rho_cosTheta_num","#rho","cos(#theta)","#rho",100,-1,1); 
	CreateHistogram("rho_cosTheta_num_uW","#rho","cos(#theta)","#rho",100,-1,1); 
	CreateHistogram("rho_cosTheta_den","#rho","cos(#theta)","#rho",100,-1,1); 

	
	
	CreateHistogram2D("rho_dR_pt_num","#rho ","p_{T} [GeV]","#DeltaR",50,10,100,50,0,3);
	CreateHistogram2D("rho_dR_pt_den","#rho Denominator","p_{T} [GeV]","#DeltaR",50,10,100,50,0,3);
	
	CreateHistogram2D("rho_dRPtE_pt_num","#rho ","p_{T} [GeV]","#DeltaR_{#Deltap_{T}}^{elliptic}",50,10,100,50,0,3);
	CreateHistogram2D("rho_dRPtE_pt_den","#rho Denominator","p_{T} [GeV]","#DeltaR_{#Deltap_{T}}^{elliptic}",50,10,100,50,0,3);
	
	CreateHistogram2D("rho_y_pt_num","#rho ","p_{T} [GeV]","|y|",50,10,100,50,0,1.25);
	CreateHistogram2D("rho_y_pt_den","#rho Denominator","p_{T} [GeV]","|y|",50,10,100,50,0,1.25);
	
	CreateHistogram2D("rho_Dphi_Deta_num","#rho #Delta#phi vs. #Delta|#eta|","#Delta|#eta|","#Delta#Phi",50,0,3.15,50,0,3); 
	CreateHistogram2D("rho_Dphi_Deta_den","#rho #Delta#phi vs. #Delta|#eta|","#Delta|#eta|","#Delta#Phi",50,0,3.15,50,0,3); 
	
	CreateHistogram2D("rho_pt_distM1_num","#rho #Delta#phi vs. p_{T}","p_{T}(#mu#mu) [GeV]","dist M1 [cm]",50,10,100,100,0,1000); 
	CreateHistogram2D("rho_pt_distM1_den","#rho #Delta#phi vs. p_{T}","p_{T}(#mu#mu) [GeV]","dist M1 [cm]",50,10,100,100,0,1000); 
	
	CreateHistogram("vertex_num","#epsilon_{vp} ","p_{T} [GeV]","#epsilon_{vp}",50,10,100);
	CreateHistogram("vertex_den","#epsilon_{vp} Denominator","p_{T} [GeV]","#epsilon_{vp}",50,10,100);
	
	CreateHistogram("sg_num","#epsilon_{sg} ","p_{T} [GeV]","#epsilon_{sg}",50,10,100);
	CreateHistogram("sg_den","#epsilon_{sg} Denominator","p_{T} [GeV]","#epsilon_{sg}",50,10,100);
	
	CreateHistogram("sg_num_dRPtE","#epsilon_{sg} ","#DeltaR_{#Deltap_{T}}^{elliptic}","#epsilon_{sg}",100,0,3);
	CreateHistogram("sg_den_dRPtE","#epsilon_{sg} Denominator","#DeltaR_{#Deltap_{T}}^{elliptic}","#epsilon_{sg}",100,0,3);
	
	CreateHistogram2D("DeltaPhi_DeltaEta_gen","Gen #Delta#phi vs. #Delta#|#eta|","#Delta|#eta|","#Delta#Phi",50,0,0.7,50,0,0.7); 
	
	CreateHistogram2D("eta1_eta2_y0","|#eta|_{1} vs |#eta|_{2}, |y|<0.6", "|#eta|_{2}","|#eta|_{1}",50,0,2.5,50,0,2.5); 
	CreateHistogram2D("eta1_eta2_y1","|#eta|_{1} vs |#eta|_{2}, 0.6<|y|<1.2", "|#eta|_{2}","|#eta|_{1}",50,0,2.5,50,0,2.5); 
	
	CreateHistogram("muon_dz","Muon |dz|","|dz|","",100,0,50);
	
	CreateHistogram("eff_mass_num","#epsilon_{M_{#mu#mu}}","M_{#mu#mu}","#epsilon",500,8.45,11.55); 
	CreateHistogram("eff_mass_den","#epsilon_{M_{#mu#mu}}","M_{#mu#mu}","#epsilon",500,8.45,11.55); 

	CreateHistogram("M_mumu","M^{gen}_{#mu#mu}","M^{gen}_{#mu#mu}","Events",2500,8.5,11.5);
	
	CreateHistogram("DeltaZeta","Delta #zeta","#Delta#Zeta [MeV]","Events",200,-100,100); 
	
	float fPTbin[]={10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,43,46,50,55,60,70,100};
	 
	const int fNpt=sizeof(fPTbin)/sizeof(float)-1;
	CreateHistogram("bin_width","bin_width","bin","width",fNpt,fPTbin); 
	for(int ipt=0; ipt<fNpt; ipt++){
		CreateHistogram2D(Form("dm_m_y0_pt%d",ipt),Form("#zeta vs M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "#zeta [MeV]",250, 8.5,11.5,4000,0,1000);
		CreateHistogram2D(Form("dm_m_y1_pt%d",ipt),Form("#zeta vs M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "#zeta [MeV]",250, 8.5,11.5,4000,0,1000); 

		CreateHistogram(Form("m_y0_pt%d",ipt),Form("M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "Events/10MeV",250, 8.5,11.5); 
		CreateHistogram(Form("m_y1_pt%d",ipt),Form("M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "Events/10MeV",250, 8.5,11.5); 
		
		//weighted 1/eff
		CreateHistogram2D(Form("dm_mw_y0_pt%d",ipt),Form("#zeta vs M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "#zeta [MeV]",250, 8.5,11.5,4000,0,1000);
		CreateHistogram2D(Form("dm_mw_y1_pt%d",ipt),Form("#zeta vs M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "#zeta [MeV]",250, 8.5,11.5,4000,0,1000); 
		
		CreateHistogram(Form("mw_y0_pt%d",ipt),Form("M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "Events/10MeV",250, 8.5,11.5); 
		CreateHistogram(Form("mw_y1_pt%d",ipt),Form("M_{#mu#mu}, %.0f<p_{T}<%.0f",fPTbin[ipt],fPTbin[ipt+1]), "M_{#mu#mu} [GeV]", "Events/10MeV",250, 8.5,11.5); 
			
		
	}
	for(int iy=0; iy<2; iy++){
		
		CreateHistogram(Form("acceptance_passed_unPol_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		CreateHistogram(Form("acceptance_total_unPol_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		
		CreateHistogram(Form("acceptance_passed_transverse_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		CreateHistogram(Form("acceptance_total_transverse_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		
		CreateHistogram(Form("acceptance_passed_longitudinal_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		CreateHistogram(Form("acceptance_total_longitudinal_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		
		CreateHistogram(Form("acceptance_passed_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		CreateHistogram(Form("acceptance_total_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		
		CreateHistogram(Form("acceptance_passedEp_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		CreateHistogram(Form("acceptance_totalEp_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		
		CreateHistogram(Form("acceptance_passedEm_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
		CreateHistogram(Form("acceptance_totalEm_y%d",iy),"Acceptance Passed","p_{T} [GeV]","A",90,10,100);
	}
	
	CreateHistogram2D("acceptance_passed","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	CreateHistogram2D("acceptance_total","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	
	CreateHistogram2D("acceptance_passed_unPol","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	CreateHistogram2D("acceptance_total_unPol","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	
	CreateHistogram2D("acceptance_passed_transverse","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	CreateHistogram2D("acceptance_total_transverse","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	
	CreateHistogram2D("acceptance_passed_longitudinal","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	CreateHistogram2D("acceptance_total_longitudinal","Acceptance Passed","p_{T} [GeV]","A",90,10,100,50,0,1.25);
	
	for(int iy=0; iy<2; iy++){
		CreateProfile(Form("lambda_theta_y%d",iy),"#lambda_{#theta}","p_{T} [GeV]","<#lambda_{#theta}>",180,10,100); 
		CreateProfile(Form("lambda_phi_y%d",iy),"#lambda_{#phi}","p_{T} [GeV]","<#lambda_{#theta#phi}>",180,10,100); 
		CreateProfile(Form("lambda_theta_phi_y%d",iy),"#lambda_{#theta#phi}","p_{T} [GeV]","<#lambda_{#theta#phi}>",180,10,100); 

		
		CreateProfile(Form("lambda_thetaEp_y%d",iy),"#lambda_{#theta}","p_{T} [GeV]","<#lambda_{#theta}>",180,10,100); 
		CreateProfile(Form("lambda_phiEp_y%d",iy),"#lambda_{#phi}","p_{T} [GeV]","<#lambda_{#theta#phi}>",180,10,100); 
		CreateProfile(Form("lambda_theta_phiEp_y%d",iy),"#lambda_{#theta#phi}","p_{T} [GeV]","<#lambda_{#theta#phi}>",180,10,100); 
		
		CreateProfile(Form("lambda_thetaEm_y%d",iy),"#lambda_{#theta}","p_{T} [GeV]","<#lambda_{#theta}>",180,10,100); 
		CreateProfile(Form("lambda_phiEm_y%d",iy),"#lambda_{#phi}","p_{T} [GeV]","<#lambda_{#theta#phi}>",180,10,100); 
		CreateProfile(Form("lambda_theta_phiEm_y%d",iy),"#lambda_{#theta#phi}","p_{T} [GeV]","<#lambda_{#theta#phi}>",180,10,100); 
	}
	
}


void JPsiAnalyzerPAT::beginJob()
{
	
	//read scale-correction parameters from the database:
	// edm::ESHandle<MuScleFitDBobject> dbObject;
	// iSetup.get<MuScleFitDBobjectRcd>().get(dbObject);
	// corrector_.reset(new MomentumScaleCorrector( dbObject.product() ) );
	
	
    //std::cout << "[JPsiAnalyzerPAT] --- beginJob " << std::endl;
	
	fEff = new TFile("src/data/MCEff2.root","READ"); 
	fRho = new TFile("src/data/rho_factor.root","READ"); 
	if(!fEff->IsOpen()){
		cout <<"Efficiency file failed to open." << endl;
		return;
	}
	bookHisto();
	OpenPolarizationFiles(); 
	if (_writeTree) {
		fOut_ = new TFile(_treefilename.c_str(), "RECREATE");
		fOut_->cd();

		// TTree
		//load Branches
		tree_ = new TTree ("data", "CMSSW Quarkonia J/psi Polarization+Trigger Tree");
		
		JpsiP = new TLorentzVector();
		muPosP = new TLorentzVector();
		muNegP = new TLorentzVector();
		JpsiP_Gen = new TLorentzVector();
		muPosP_Gen = new TLorentzVector();
		muNegP_Gen = new TLorentzVector();
		
		// Event variables
		tree_->Branch("eventNb",             &eventNb,             "eventNb/I");
		tree_->Branch("runNb",               &runNb,               "runNb/I");
		tree_->Branch("lumiBlock",           &lumiBlock,           "lumiBlock/I");
		tree_->Branch("nPriVtx",             &nPriVtx,             "nPriVtx/I");
		tree_->Branch("vertexWeight",         &vertexWeight,         "vertexWeight/D");
		tree_->Branch("sumPTPV",             &sumPTPV,             "sumPTPV/D");
		tree_->Branch("countTksOfPV",        &countTksOfPV,        "countTksOfPV/I");
		
		// Jpsi Variables
		//     tree_->Branch("JpsiType",   &JpsiType,  "JpsiType/I");
		tree_->Branch("JpsiP",  "TLorentzVector", &JpsiP);
		tree_->Branch("JpsiMuScleMassCorr",   &JpsiMuScleMassCorr,  "JpsiMuScleMassCorr/D");
		tree_->Branch("JpsiMuScleMassErr",   &JpsiMuScleMassErr,  "JpsiMuScleMassErr/D");
		tree_->Branch("JpsiMassErr",   &JpsiMassErr,  "JpsiMassErr/D");
		tree_->Branch("sigmaPtPos",   &sigmaPtPos,  "sigmaPtPos/D");
		tree_->Branch("sigmaPtNeg",   &sigmaPtNeg,  "sigmaPtNeg/D");
		// tree_->Branch("JpsiPt",     &JpsiPt,    "JpsiPt/D");
		// tree_->Branch("JpsiRap",    &JpsiRap,   "JpsiRap/D");
		//     tree_->Branch("JpsiCharge", &JpsiCharge,"JpsiCharge/I");
		// tree_->Branch("JpsiPx",     &JpsiPx,    "JpsiPx/D");
		// tree_->Branch("JpsiPy",     &JpsiPy,    "JpsiPy/D");
		// tree_->Branch("JpsiPz",     &JpsiPz,    "JpsiPz/D");
		tree_->Branch("Jpsict",     &Jpsict,    "Jpsict/D");
		tree_->Branch("JpsictErr",  &JpsictErr, "JpsictErr/D");
		tree_->Branch("JpsiVprob",  &JpsiVprob, "JpsiVprob/D");
		//muon distance 
		tree_->Branch("JpsiDistM1",   &JpsiDistM1,    "JpsiDistM1/D");
		tree_->Branch("JpsiDphiM1",   &JpsiDphiM1,    "JpsiDphiM1/D");
		tree_->Branch("JpsiDrM1",     &JpsiDrM1,      "JpsiDrM1/D");
		tree_->Branch("JpsiDistM2",   &JpsiDistM2,    "JpsiDistM2/D");
		tree_->Branch("JpsiDphiM2",   &JpsiDphiM2,    "JpsiDphiM2/D");
		tree_->Branch("JpsiDrM2",     &JpsiDrM2,      "JpsiDrM2/D");
		
		tree_->Branch("muPosP", "TLorentzVector", &muPosP);
		tree_->Branch("muNegP", "TLorentzVector", &muNegP);
		// tree_->Branch("muPosPx",    &muPosPx,   "muPosPx/D");
		// tree_->Branch("muPosPy",    &muPosPy,   "muPosPy/D");
		// tree_->Branch("muPosPz",    &muPosPz,   "muPosPz/D");
		// tree_->Branch("muNegPx",    &muNegPx,   "muNegPx/D");
		// tree_->Branch("muNegPy",    &muNegPy,   "muNegPy/D");
		// tree_->Branch("muNegPz",    &muNegPz,   "muNegPz/D");
		
		tree_->Branch("muPos_qual",&muPos_qual,"muPos_qual/I");
		tree_->Branch("muNeg_qual",&muNeg_qual,"muNeg_qual/I");
		tree_->Branch("QQ_check",&QQ_check,"QQ_check/I");
		tree_->Branch("TroggerBit",&TriggerBit,"TriggerBit/I");
		tree_->Branch("NCand",&NCand,"NCand/I");

		
		tree_->Branch("DCA",&dca,"DCA/D");
		
		//add HLT Variables to TTree
		for(std::vector< std::string >:: iterator it = HLTBitNames_.begin(); it != HLTBitNames_.end(); ++it){
			std::string hlt_name= *it;
			tree_->Branch(hlt_name.c_str(), &(mapTriggerNameToIntFired_[*it]), (hlt_name + "/I").c_str());
			tree_->Branch((hlt_name+"_PreScale").c_str(), &(mapTriggerNameToPrescaleFac_[*it]), (hlt_name+"_PreScale" + "/I").c_str());
		}
		
		//add Generator Information
		if(_isMC){
			tree_->Branch("MCType",         &MCType,        "MCType/I");
			tree_->Branch("JpsiP_Gen",  "TLorentzVector", &JpsiP_Gen);
			tree_->Branch("Jpsict_Gen",     &Jpsict_Gen,    "Jpsict_Gen/D");
			tree_->Branch("muPosP_Gen",  "TLorentzVector", &muPosP_Gen);
			tree_->Branch("muNegP_Gen",  "TLorentzVector", &muNegP_Gen);
		}
	}
	if (_writeDataSet) {
		
		fOut2_ = new TFile(_datasetname.c_str(), "RECREATE");
		fOut2_->cd();
		
		Jpsi_PtType = new RooCategory("Jpsi_PtType","Category of Pt");
		Jpsi_RapType = new RooCategory("Jpsi_RapType","Category of Rap");
		
		JpsiPtMin = _ptbinranges[0];  cout << "Pt min = " << JpsiPtMin << endl;
		JpsiPtMax = _ptbinranges[_ptbinranges.size()-1];  cout << "Pt max = " << JpsiPtMax << endl;
		
		for(unsigned int i=0;i<_ptbinranges.size()-1;i++){
			char catname[100];
			sprintf(catname,"P%d",i+1);
			Jpsi_PtType->defineType(catname,i+1); 
			cout << "Pt bin " << i+1 << ": Min = " << _ptbinranges[i] << " Max = " << _ptbinranges[i+1] << endl;   
		}
		
		JpsiRapMin = _etabinranges[0];  cout << "Rap min = " << JpsiRapMin << endl;
		JpsiRapMax = _etabinranges[_etabinranges.size()-1];  cout << "Rap max = " << JpsiRapMax << endl;
		
		for(unsigned int i=0;i<_etabinranges.size()-1;i++){
			char catname[100];
			sprintf(catname,"E%d",i+1);
			Jpsi_RapType->defineType(catname,i+1); 
			cout << "Rap bin " << i+1 << ": Min = " << _etabinranges[i] << " Max = " << _etabinranges[i+1] << endl;   
		}
		
		Jpsi_MatchType = new RooCategory("Jpsi_MatchType","Category of matching");
		
		Jpsi_MatchType->defineType("unmatched",0);
		Jpsi_MatchType->defineType("matched",1);
		
		Jpsi_MuScleMass = new RooRealVar("Jpsi_MuScleMass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
		Jpsi_MuScleMassErr = new RooRealVar("Jpsi_MuScleMassErr","J/psi mass error",0,1.,"GeV/c^{2}");
		Jpsi_MassErr = new RooRealVar("Jpsi_MassErr","J/psi vtx mass error",0,1.,"GeV/c^{2}");
		Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/psi pt",JpsiPtMin,JpsiPtMax,"GeV/c");
		Jpsi_Rap = new RooRealVar("Jpsi_Rap","J/psi eta",-JpsiRapMax,JpsiRapMax);
		Jpsi_ct = new RooRealVar("Jpsi_ct","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
		Jpsi_ctErr = new RooRealVar("Jpsi_ctErr","J/psi ctau error",-1.,1.,"mm");
		Jpsi_ctTrue = new RooRealVar("Jpsi_ctTrue","J/psi ctau true",-100.,JpsiCtMax,"mm"); 		
		
		RooArgList varlist(*Jpsi_MuScleMass,*Jpsi_MuScleMassErr, *Jpsi_MassErr,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap,*Jpsi_MatchType);
		varlist.add(*Jpsi_ctTrue);   varlist.add(*Jpsi_PtType);
		varlist.add(*Jpsi_RapType);  varlist.add(*Jpsi_ctErr);
		
		data = new RooDataSet("data","A sample",varlist);
	}
}

double JPsiAnalyzerPAT::correction_functionJPsi(const double & pt, const double & eta, const double & phi, const int chg, double *par) {
	// correction function JPsi 
	double etaCorr = 0.;
	if( (chg < 0 && eta > par[4]) || (chg > 0 && eta < -par[4]) ) {
		etaCorr = par[1]+par[2]*fabs(fabs(eta)-par[4])+par[3]*(fabs(eta)-par[4])*(fabs(eta)-par[4]);
	}
	
	//don't forget par[3] = 1.6
	
	double ptCorr = 0.;
	if( pt < par[7] ) {
		ptCorr = par[5]*(pt - par[7]) + par[6]*(pt - par[7])*(pt - par[7]);
	}
	
	//don't forget par[6] = 6
	
	return par[0]*pt*(1 + etaCorr + ptCorr);
} 

double JPsiAnalyzerPAT::correction_functionZ(const double & pt, const double & eta, const double & phi, const int chg,   double *parScale) {    
    double ampl(0), phase(0), twist(0), ampl2(0), freq2(0), phase2(0);
	// very bwd bin
    if ( eta  < parScale[4] ) {
		ampl = parScale[1]; phase = parScale[2]; ampl2 = parScale[21]; freq2 = parScale[22]; phase2 = parScale[23];
		twist = parScale[3]*(eta-parScale[4])+parScale[7]*(parScale[4]-parScale[8])+parScale[11]*parScale[8]; 
		// bwd bin
    } else if ( parScale[4] <= eta && eta < parScale[8] ) {
		ampl = parScale[5]; phase = parScale[6];
		twist = parScale[7]*(eta-parScale[8])+parScale[11]*parScale[8] ; 
		// barrel bin
    } else if ( parScale[8] <= eta && eta < parScale[12] ) {
		ampl = parScale[9]; phase = parScale[10];
		twist = parScale[11]*eta; 
		// fwd bin
    } else if ( parScale[12] <= eta && eta < parScale[16] ) {
		ampl = parScale[13]; phase = parScale[14];
		twist = parScale[15]*(eta-parScale[12])+parScale[11]*parScale[12]; 
		// very fwd bin
    } else if ( parScale[16] < eta ) {
		ampl = parScale[17]; phase = parScale[18]; ampl2 = parScale[24]; freq2 = parScale[25]; phase2 = parScale[26];
		twist = parScale[19]*(eta-parScale[16])+parScale[15]*(parScale[16]-parScale[12])+parScale[11]*parScale[12]; 
    }
    
    // apply the correction
    double curv = (1.+parScale[0])*((double)chg/pt
									-twist
									-ampl*sin(phi+phase)
									-ampl2*sin(freq2*phi+phase2)
									-0.5*parScale[20]);
    return 1./((double)chg*curv);
}

double JPsiAnalyzerPAT::Pt_primeZ(TVector3 v1,int Q){
	//Z corrections on 2011 Data 
	double parScale[]={-0.000322,0.001403,0.98565,-0.000792,-2.1,0.000535,0.825151,-0.000234,-1.5,0.000148,-1.47645,-2.8e-05,1.5,0.000182,0.84186,-0.000425,2.1,0.000733,1.93751,-0.001493,4.8e-05,0.00048,2,1.79589,0.000586,2,-0.518118};
	double parScale_MC[]={0.00123,0.001113,0.168873,0.000124,-2.1,0.000336,-0.999195,0.000174,-1.5,0.000177,-1.76348,2.6e-05,1.5,0.000377,-1.08773,-0.000101,2.1,0.000879,-1.43165,6.4e-05,4.4e-05,0.00087,2,-1.94623,0.000651,2,2.54445};
	
	double muon_pt=v1.Perp();
	double muon_eta=v1.Eta();
	double muon_phi=v1.Phi();
	double corrected_pt=0; 
	if(!_isMC){
		//DATA
		corrected_pt=correction_functionZ(muon_pt,muon_eta,muon_phi,Q,parScale);
	}
	if(_isMC){
		corrected_pt=correction_functionZ(muon_pt,muon_eta,muon_phi,Q,parScale_MC);
	}
	return corrected_pt; 
	
}

double JPsiAnalyzerPAT::Pt_primeJPsi(TVector3 v1,int Q){
	//Momentum correction using JPsi, Pt<10GeV, maybe extrapolate to 20GeV? 
	
	//Mode==4 (data)	
	double parScaleJPsi_it1[]={1.00074,0.000103333, -0.000947069, -0.00723064, 1.65256, 0., 0.000181535, 6}; // JPsi full 2011 data set, corrections for first iteration
	double parScaleJPsi_it2[]={1.00014, 2.2114e-05, -6.29854e-05,-0.000515906, 1.65164, 0., -7.40772e-06, 6}; // JPsi corrections for 2nd iteration
	
	double parScaleJPsi_MC[]={1.00006,-0.00042742,0.00011569,0.00393103,1.60458,0,5.87978e-05,6};
	
	double muon_pt=v1.Perp();
	double muon_eta=v1.Eta();
	double muon_phi=v1.Phi();
	double corrected_pt=0; 
	if(!_isMC){
		//Data
		corrected_pt=correction_functionJPsi(muon_pt,muon_eta,muon_phi,Q,parScaleJPsi_it1);
		corrected_pt=correction_functionJPsi(corrected_pt,muon_eta,muon_phi,Q,parScaleJPsi_it2); 
	}
	if(_isMC){
		corrected_pt=correction_functionJPsi(muon_pt,muon_eta,muon_phi,Q,parScaleJPsi_MC);
	}
	
	return corrected_pt; 
}


double JPsiAnalyzerPAT::CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode, TLorentzVector &DiMuon, TVector3 &Mu1, TVector3 &Mu2){  
	
	//MuScle Fit corrections
	//1) correct the momentum scale:
	
	bool print=0; 
	
	double corrPt1 = (*corrector_)(mu1);
	double corrPt2 = (*corrector_)(mu2);
	
	
	
	
	//   cout << "original pT1 " << (mu1.innerTrack()->momentum()).Rho() << " corrected pT1 " << corrPt1 << endl;
	//   cout << "original pT2 " << (mu2.innerTrack()->momentum()).Rho() << " corrected pT2 " << corrPt2 << endl;
	
	const double mumass = 0.105658;
	TLorentzVector mu1Corr, mu2Corr; 
	mu1Corr.SetPtEtaPhiM(corrPt1, mu1.innerTrack()->eta(), mu1.innerTrack()->phi(), mumass);
	mu2Corr.SetPtEtaPhiM(corrPt2, mu2.innerTrack()->eta(), mu2.innerTrack()->phi(), mumass);
	TLorentzVector onia = mu1Corr+mu2Corr;
	JpsiMuScleMassCorr = onia.M();
	
	//JPsiZ method
	TLorentzVector mu1Corr_M2, mu2Corr_M2; //method 2
	TVector3 p1; 
	TVector3 p2; 
	
	p1.SetPtEtaPhi(mu1.innerTrack()->pt(),mu1.innerTrack()->eta(),mu1.innerTrack()->phi()); 
	p2.SetPtEtaPhi(mu2.innerTrack()->pt(),mu2.innerTrack()->eta(),mu2.innerTrack()->phi()); 
	double pt_prime=0; 
	
	if(p1.Perp()<=15) pt_prime=Pt_primeJPsi(p1,mu1.charge()); 
	else pt_prime=Pt_primeZ(p1,mu1.charge()); 
	mu1Corr_M2.SetPtEtaPhiM(pt_prime,mu1.innerTrack()->eta(),mu1.innerTrack()->phi(),mumass); 
	if(p2.Perp()<=15) pt_prime=Pt_primeJPsi(p2,mu2.charge()); 
	else pt_prime=Pt_primeZ(p2,mu2.charge()); 
	mu2Corr_M2.SetPtEtaPhiM(pt_prime,mu2.innerTrack()->eta(),mu2.innerTrack()->phi(),mumass); 
	
	p1=mu1Corr_M2.Vect();
	p2=mu2Corr_M2.Vect();
	
	TLorentzVector Upsilon = mu1Corr_M2+mu2Corr_M2; 
	if(print){
		cout << "Delta Onia-Upsilon: " << Upsilon.M()-onia.M() << endl;
		cout << "ptCorr: " << onia.Perp() << endl; 
		cout << "ptCorr JPsi/Z: " << Upsilon.Perp() << endl; 
	}
	DiMuon=Upsilon; 
	Mu1=p1; 
	Mu2=p2; 
	
	//2) calculate the error on the mass
	double ptEtaPhiE_1[4] = {corrPt1, mu1.innerTrack()->eta(), mu1.innerTrack()->phi(), 0.};//E will be calculated automatically
	double ptEtaPhiE_2[4] = {corrPt2, mu2.innerTrack()->eta(), mu2.innerTrack()->phi(), 0.};//E will be calculated automatically
	
	JpsiMuScleMassErr = MuScleFitUtils::massResolution(MuScleFitUtils::fromPtEtaPhiToPxPyPz(ptEtaPhiE_1), 
													   MuScleFitUtils::fromPtEtaPhiToPxPyPz(ptEtaPhiE_2), 
													   *resolutionFunction_);
	
	//3) save also the resolution on the single muon pT
	double sigmaPT1 = resolutionFunction_->sigmaPt(mu1, 0);
	double sigmaPT2 = resolutionFunction_->sigmaPt(mu2, 0);
	if(mu1.charge() > 0){sigmaPtPos = sigmaPT1; sigmaPtNeg = sigmaPT2;}
	else if(mu1.charge() < 0){sigmaPtPos = sigmaPT1; sigmaPtNeg = sigmaPT1;}
	
	//cout << "mass " << JpsiMassCorr << " errMass " << JpsiMassErr << endl;
	return JpsiMuScleMassCorr;
	
}

// ------------ method called to for each event  ------------
void
JPsiAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	bool print=0; 
	nEvents++;
	if(print) std::cout << "Event: " << nEvents << std::endl; 
	// reset TTree Variables
	if (_writeTree) this->resetDSVariables();
	
	// check HLT TriggerReuslts
	this->hltReport(iEvent, iSetup);
	
	bool trigOK = false;
	for (unsigned int iTrig = 0 ; iTrig < HLTBitNames_.size() ; iTrig++) {
		if (mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == 3) trigOK = true;
	}
	if(trigOK)TriggerBit=1; 
	else TriggerBit=0;
	
	if (requireTriggerMatching_ && !trigOK && !_storeAllMCEvents) return;
	
	// Event related infos
	eventNb= iEvent.id().event() ;
	runNb=iEvent.id().run() ;
	lumiBlock= iEvent.luminosityBlock() ;
	
	Handle<reco::VertexCollection> privtxs;
	iEvent.getByLabel("offlinePrimaryVertices", privtxs);
	nPriVtx = privtxs->size();
	VertexCollection::const_iterator privtx;
	
	if ( privtxs->begin() != privtxs->end() ) {
		privtx=privtxs->begin();
		RefVtx = privtx->position();
	} else {
		RefVtx.SetXYZ(0.,0.,0.);
	}
	// }
	
	try {iEvent.getByLabel(_patJpsi,collAll);} 
	catch (...) {cout << "J/psi not present in event!" << endl;}
	
	_thePassedCands.clear();
	
	// APPLY CUTS
	this->makeCuts();
	
	bool storeEvent = false;
	if(print) std::cout << "cand dize: " << _thePassedCands.size() << std::endl; 
	// BEST J/PSI? 
	if (_onlythebest && _thePassedCands.size()>0) {  // yes, fill simply the best
		int iBest = theBestQQ();
		if (iBest > -1){
			fillTreeAndDS(_thePassedCands.at(iBest), iEvent);
			//fillHistograms(_thePassedCands.at(count),iEvent);
			passedMuonSelectionCuts_++;
			storeEvent=true;
		}
	} else {   // no, fill all candidates passing cuts (possibly wrong-sign)
		for( unsigned int count = 0; count < _thePassedCands.size(); count++) { 
			fillTreeAndDS(_thePassedCands.at(count),iEvent);
			//fillHistograms(_thePassedCands.at(count),iEvent);
			passedMuonSelectionCuts_++;
		}
	}
	
	//! FILL GENERATOR COLLECTION and store the event
	if ( _storeAllMCEvents || storeEvent ) {
		if(print) cout << "store Event: " << endl; 
		Handle<reco::GenParticleCollection> genParticles;
		iEvent.getByLabel( _genParticles, genParticles );
		if ( genParticles.isValid() )
		{
			//std::cout << "------ analyze GENERATED JPsis:" << std::endl;
			this->analyzeGenerator( genParticles );
		}
				
		
		// Write all Branches to the Tree ONLY 
		// - for the best candidate
		// - for the opposite sign
		
	
		if (_writeTree) tree_->Fill();
	}
	
	if(_isMC){
		//MC related histograms 
		if(print)cout << "Fill Gen level Mass: " << endl; 
		TLorentzVector pMP;
		pMP.SetPxPyPzE(muPosP_Gen->Px(),muPosP_Gen->Py(),muPosP_Gen->Pz(),muPosP_Gen->E());
		TLorentzVector pMN;
		pMN.SetPxPyPzE(muNegP_Gen->Px(),muNegP_Gen->Py(),muNegP_Gen->Pz(),muNegP_Gen->E());
		TLorentzVector Mmumu=pMP+pMN;
		
		hName["M_mumu"]->Fill(Mmumu.M()); 
		if(print) cout << " size of cand collection: " << _thePassedCands.size() << endl; 
		for( unsigned int count = 0; count < _thePassedCands.size(); count++) { 
			if(print)cout <<  "fillHistograms: " << endl; 
			fillHistograms(_thePassedCands.at(count),iEvent);
		}
		
	}//_isMC

	
	
}

double JPsiAnalyzerPAT::M_uncertainty(TVector3 fMuon1Vect, TVector3 fMuon2Vect, double *fMuon1E, double *fMuon2E){
	
	double MMUON=0.105658;
	
	double pt1 = fMuon1Vect.Perp();
	double pt2 = fMuon2Vect.Perp();
	
	double fMuon1Eta=fMuon1Vect.Eta();
	double fMuon2Eta=fMuon2Vect.Eta();
	double fMuonPhi12=fMuon1Vect.DeltaPhi(fMuon2Vect); 
	
	double Pl1 = fMuon1Vect.Pz();
	double Pl2 = fMuon2Vect.Pz();
	
	double fMuon1PtE=fMuon1E[0];
	double fMuon2PtE=fMuon2E[0]; 
	
	double fMuon1EtaE=fMuon1E[1];
	double fMuon2EtaE=fMuon2E[1];
	
	double fMuon1PhiE=fMuon1E[2];
	double fMuon2PhiE=fMuon2E[2]; 
	
	double dphi=fMuonPhi12; 
	double dphiE = sqrt(fMuon1PhiE*fMuon1PhiE+fMuon2PhiE*fMuon2PhiE);
	
	double sigma_Pl1 = sqrt(pow(TMath::SinH(fMuon1Eta)*fMuon1PtE,2)+pow(TMath::CosH(fMuon1Eta)*pt1*fMuon1EtaE,2));
	double sigma_Pl2= sqrt(pow(TMath::SinH(fMuon2Eta)*fMuon2PtE,2)+pow(TMath::CosH(fMuon2Eta)*pt2*fMuon2EtaE,2));
	
	double sigmaE1=sqrt(pow(pt1*fMuon1PtE,2)+pow(Pl1*sigma_Pl1,2));
	
	double sigmaE2=sqrt(pow(pt2*fMuon2PtE,2)+pow(Pl2*sigma_Pl2,2));
	
	double E1=sqrt(fMuon1Vect.Mag()*fMuon1Vect.Mag()+MMUON*MMUON);
	
	double E2=sqrt(fMuon2Vect.Mag()*fMuon2Vect.Mag()+MMUON*MMUON);
	
	double dm_pt1 = pt1*pow(TMath::CosH(fMuon1Eta),2)*(E2/E1)-pt2*TMath::Cos(dphi)-TMath::SinH(fMuon1Eta)*Pl2;
	double dm_pt2 = pt2*pow(TMath::CosH(fMuon2Eta),2)*(E1/E2)-pt1*TMath::Cos(dphi)-TMath::SinH(fMuon2Eta)*Pl1;
	
	double dm_eta1 = fMuon1Vect.Mag()*(Pl1*(E2/E1)-Pl2); 
	double dm_eta2 = fMuon2Vect.Mag()*(Pl2*(E1/E2)-Pl1);
	
	double dm_phi = pt1*pt2*TMath::Sin(dphi);
	
	double M = sqrt(2)*sqrt(MMUON*MMUON+E1*E2-pt1*pt2*TMath::Cos(dphi)-Pl1*Pl2);
	double sigmaM = (1/M)*sqrt(pow(dm_pt1*fMuon1PtE,2)+pow(dm_pt2*fMuon2PtE,2)+pow(dm_eta1*fMuon1EtaE,2)+pow(dm_eta2*fMuon2EtaE,2)+pow(dm_phi*dphiE,2));
	return sigmaM;
}


void JPsiAnalyzerPAT::error_graph(TGraph *gr, TString name){
	//cout << "Create Error graph: " << name << endl; 
	
	TGraph *grEp = new TGraph(gr->GetN()); 
	grEp->SetName(name+"P"); 

	TGraph *grEm = new TGraph(gr->GetN()); 
	grEm->SetName(name+"M"); 
	double *x_ = gr->GetX();
	double *y_ = gr->GetY();
	
	for(int i=0; i<gr->GetN(); i++){
		double EY=gr->GetErrorY(i);
	//	cout << "x: " << x_[i] << "yE+: " << y_[i]+EY << endl;
		grEp->SetPoint(i,x_[i],y_[i]+EY);
		grEm->SetPoint(i,x_[i],y_[i]-EY);
	}
	
	polName[grEp->GetName()]=grEp; 
	polName[grEm->GetName()]=grEm; 
	
}

void JPsiAnalyzerPAT::OpenPolarizationFiles(){

	TFile *F; 
	TGraph *gr; 
	std::map<TString, TFile*> FName; 
	
	for(int ups=1; ups<=3; ups++){
		for(int sigma=1; sigma<=3; sigma++){
			TString file_name=Form("src/data/TGraphResults_%dSUps_%dsigma.root",ups,sigma);
			F = new TFile(file_name,"READ"); 
			if(F->IsOpen()!=1) {
				cout << "File: " << file_name << " Failed to open." << endl; 
				return; 
			}
			FName[file_name]=F;
		}
		TString file_name=Form("src/data/TGraphResults_%dSUps_SU_1sigma.root",ups);
		F = new TFile(file_name,"READ"); 
	
		if(F->IsOpen()!=1) {
			cout << "File: " << file_name << " Failed to open." << endl; 
			return; 
		}
		
		FName[file_name]=F; 
	}
	
	
	string frame_list[]={"HX","CS","PX"};
	
	for(int ups=1; ups<=3; ups++){
		for(int sigma=1; sigma<=3; sigma++){
			TString name=Form("src/data/TGraphResults_%dSUps_%dsigma.root",ups,sigma); 
		//	cout << "File: " << name << endl; 
			for(int iy=1; iy<=2; iy++){
				for(int frame=0; frame<3; frame++){
					
					//rapidity |y|<0.6 and 0.6<|y|<1.2
					//lambda theta, lambda phi, lambda theta phi 
					// computed in HX,CS and PX frames
				//	cout << Form("lth_%s_rap%d",frame_list[frame].c_str(),iy) << endl; 
					gr = (TGraph*)FName[name]->Get(Form("lth_%s_rap%d",frame_list[frame].c_str(),iy));
					polName[Form("lth_%s_%dSUps_%dSigma_y%d",frame_list[frame].c_str(),ups,sigma,iy-1)]=gr; 
					error_graph(gr,Form("lth_%s_%dSUps_%dSigma_y%d",frame_list[frame].c_str(),ups,sigma,iy-1)); 
					
				//	cout << Form("lph_%s_rap%d",frame_list[frame].c_str(),iy) << endl; 
					gr = (TGraph*)FName[name]->Get(Form("lph_%s_rap%d",frame_list[frame].c_str(),iy));
					polName[Form("lph_%s_%dSUps_%dSigma_y%d",frame_list[frame].c_str(),ups,sigma,iy-1)]=gr; 
					error_graph(gr,Form("lph_%s_%dSUps_%dSigma_y%d",frame_list[frame].c_str(),ups,sigma,iy-1)); 
					
				//	cout << Form("ltp_%s_rap%d",frame_list[frame].c_str(),iy) << endl; 
					gr = (TGraph*)F->Get(Form("ltp_%s_rap%d",frame_list[frame].c_str(),iy));
					polName[Form("ltp_%s_%dSUps_%dSigma_y%d",frame_list[frame].c_str(),ups,sigma,iy-1)]=gr; 
					error_graph(gr,Form("ltp_%s_%dSUps_%dSigma_y%d",frame_list[frame].c_str(),ups,sigma,iy-1)); 
				}
				
			}//rapidity bin
		}//sigma loop
	}//upsilon loon 
	
	for (std::map<TString,TGraph*>::iterator it=polName.begin(); it!=polName.end(); it++) {
		//cout << "TGraphs in memory: " << it->first << endl; 
	}
	
}

void JPsiAnalyzerPAT::HX_boost(){
	bool print=0; 
	if(print) cout << "HX_boost: " << endl; 
	TLorentzVector gCand; 
	TLorentzVector HXMuPlus; 
	TLorentzVector h1; h1.SetPxPyPzE(0,1,0,0); 
	
	gCand.SetPtEtaPhiE(JpsiP_Gen->Perp(),JpsiP_Gen->Eta(),JpsiP_Gen->Phi(),JpsiP_Gen->E()); 
	HXMuPlus.SetPtEtaPhiE(muPosP_Gen->Perp(),muPosP_Gen->Eta(),muPosP_Gen->Phi(), muPosP_Gen->E()); 
	
	TLorentzRotation boost(-gCand.BoostVector()); 
	HXMuPlus *=boost; 
	
	 cosTheta = HXMuPlus.Vect().Dot(gCand.Vect())/(HXMuPlus.Vect().Mag()*gCand.Vect().Mag());
	 sinTheta = HXMuPlus.Vect().Cross(gCand.Vect()).Mag()/(HXMuPlus.Vect().Mag()*gCand.Vect().Mag());
//	double sin2Theta = 2*sinTheta*cosTheta;
	
	h1 *= boost; //h2 *= boost; h3 *= boost;
	 cosPhi = HXMuPlus.Vect().Dot(h1.Vect().Unit())/HXMuPlus.Vect().Mag();
//	double cos2Phi = 2*cosPhi*cosPhi - 1;
	
	if(print) cout << "cosTheta: " << cosTheta << " sinTheta: " << sinTheta << " cosPhi: " << cosPhi << endl; 
	
	
}

void JPsiAnalyzerPAT::FillAcceptance(double w, bool pass, string method){
	int iY=-1; 
	if(fabs(JpsiP_Gen->Rapidity())<0.6)iY=0; 
	else iY=1; 
	
	bool print=0; 
	
	if(print) cout << "FillAcceptance: " << endl; 
	
	TString pass_name = Form("acceptance_passed%s_y%d",method.c_str(),iY); 
	TString total_name = Form("acceptance_total%s_y%d",method.c_str(),iY); 

	if(print)cout << "Fill 1D acceptance: " << endl; 
	if(pass)hName[pass_name]->Fill(JpsiP_Gen->Perp(),w); 
	hName[total_name]->Fill(JpsiP_Gen->Perp(),w); 
	
	if(method=="Ep" || method=="Em") method=""; 
	
	pass_name=Form("acceptance_passed%s",method.c_str());
	total_name=Form("acceptance_total%s",method.c_str()); 

	if(print) cout << "Fill 2D acceptance: " << endl; 
	if(pass)hName2D[pass_name]->Fill(JpsiP_Gen->Perp(),fabs(JpsiP_Gen->Rapidity()),w);
	hName2D[total_name]->Fill(JpsiP_Gen->Perp(),fabs(JpsiP_Gen->Rapidity()),w); 
	
}

void JPsiAnalyzerPAT::acceptance(){
	bool print=0; 
	if(print)cout << "acceptance: " << endl; 
	double cos2Phi=2*cosPhi*cosPhi-1; 
	double sin2Theta=2*sinTheta*cosTheta; 
	
	string frame="HX"; 
	int sigma=1; 
	int iY=-1; 
	if(fabs(JpsiP_Gen->Rapidity())<0.6)iY=0; 
	else iY=1; 
	int ups=-1; 
	   
	if(_oniaPDG==553) ups=1; 
	if(_oniaPDG==10553) ups=2; 
	if(_oniaPDG==20553) ups=3; 
	
	TString nameTheta=Form("lth_%s_%dSUps_%dSigma_y%d",frame.c_str(),ups,sigma,iY);
	TString nameThetaEp=Form("lth_%s_%dSUps_%dSigma_y%dP",frame.c_str(),ups,sigma,iY);
	TString nameThetaEm=Form("lth_%s_%dSUps_%dSigma_y%dM",frame.c_str(),ups,sigma,iY);

	TString namePhi=Form("lph_%s_%dSUps_%dSigma_y%d",frame.c_str(),ups,sigma,iY);
	TString namePhiEp=Form("lph_%s_%dSUps_%dSigma_y%dP",frame.c_str(),ups,sigma,iY);
	TString namePhiEm=Form("lph_%s_%dSUps_%dSigma_y%dM",frame.c_str(),ups,sigma,iY);

	TString nameThetaPhi=Form("ltp_%s_%dSUps_%dSigma_y%d",frame.c_str(),ups,sigma,iY);
	TString nameThetaPhiEp=Form("ltp_%s_%dSUps_%dSigma_y%dP",frame.c_str(),ups,sigma,iY);
	TString nameThetaPhiEm=Form("ltp_%s_%dSUps_%dSigma_y%dM",frame.c_str(),ups,sigma,iY);

	if(print) cout << "Name theta: " << nameTheta << endl; 
	
	double pt=JpsiP_Gen->Perp(); 
	double bin_center = polName[nameTheta]->GetX()[polName[nameTheta]->GetN()-1]; 
	if(pt>bin_center)pt=bin_center; // only use value up to bin center of last bin
	
	if(print) cout << "Lambda-Theta: " << endl; 
	
	double lambdaTheta = polName[nameTheta]->Eval(pt);
	double lambdaTheta_Ep = polName[nameThetaEp]->Eval(pt);
	double lambdaTheta_Em = polName[nameThetaEm]->Eval(pt);

	if(print) cout << "Lambda-Phi: " << endl; 
	
	double lambdaPhi = polName[namePhi]->Eval(pt);
	double lambdaPhi_Ep = polName[namePhiEp]->Eval(pt);
	double lambdaPhi_Em = polName[namePhiEm]->Eval(pt);
	
	if(print) cout << "Lambda theta phi: " << endl; 
	
	double lambdaThetaPhi = polName[nameThetaPhi]->Eval(pt);
	double lambdaThetaPhi_Ep = polName[nameThetaPhiEp]->Eval(pt);
	double lambdaThetaPhi_Em = polName[nameThetaPhiEm]->Eval(pt);

	prName[Form("lambda_theta_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaTheta,1); 
	prName[Form("lambda_phi_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaPhi,1); 
	prName[Form("lambda_theta_phi_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaThetaPhi,1); 
	
	prName[Form("lambda_thetaEp_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaTheta_Ep,1); 
	prName[Form("lambda_phiEp_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaPhi_Ep,1); 
	prName[Form("lambda_theta_phiEp_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaThetaPhi_Ep,1); 
	
	prName[Form("lambda_thetaEm_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaTheta_Em,1); 
	prName[Form("lambda_phiEm_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaPhi_Em,1); 
	prName[Form("lambda_theta_phiEm_y%d",iY)]->Fill(JpsiP_Gen->Perp(),lambdaThetaPhi_Em,1); 
	
	if(print) cout << "pT: " << pt << " lambdaTheta: " << lambdaTheta << " lambdaPhi: " << lambdaPhi << " lambdaThetaPhi: " << lambdaThetaPhi << endl; 
	
	//compute weights from polarization parameters 
	double w1 = ((3/(4*TMath::Pi()))/(3+lambdaTheta))*(1+(lambdaTheta*cosTheta*cosTheta)+(lambdaPhi*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi*sin2Theta*cosPhi));
	double w1Ep = ((3/(4*TMath::Pi()))/(3+lambdaTheta_Ep))*(1+(lambdaTheta_Ep*cosTheta*cosTheta)+(lambdaPhi_Ep*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi_Ep*sin2Theta*cosPhi));
	double w1Em = ((3/(4*TMath::Pi()))/(3+lambdaTheta_Em))*(1+(lambdaTheta_Em*cosTheta*cosTheta)+(lambdaPhi_Em*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi_Em*sin2Theta*cosPhi));
	
	if(print) cout << "w: " << w1 << " wEp: " << w1Ep << " wEm: " << w1Em << endl; 
	
	bool muon1_good = false; 
	bool muon2_good = false; 
	
	bool pass = false; 
	
	double etaP=TMath::Abs(muPosP_Gen->Eta());
	double etaM=TMath::Abs(muNegP_Gen->Eta());
	
	double PtP=muPosP_Gen->Perp();
	double PtM=muNegP_Gen->Perp();
	
	if( ((etaP<1.2 && PtP>4.5) || (etaP>1.2 && etaP<1.4 && PtP>3.5) || (etaP>1.4 && etaP<1.6 && PtP>3.0)) && etaP<1.6 )
		muon1_good=true;
	if( ((etaM<1.2 && PtM>4.5) || (etaM>1.2 && etaM<1.4 && PtM>3.5) || (etaM>1.4 && etaM<1.6 && PtM>3.0)) && etaM<1.6 )
		muon2_good=true;
	
	if(muon1_good && muon2_good) pass = true; 
	
	if(print) cout << "Fill Acceptance measured polarization: " << endl; 
	FillAcceptance(w1,pass,""); //fill with measured weights
	FillAcceptance(w1Ep,pass,"Ep"); // fill with measured weights + 1 sigma
	FillAcceptance(w1Em,pass,"Em"); // fill with measured weights - 1 sigma
	
	if(print) cout << "Fill Acceptance unpolarized: " << endl; 
	FillAcceptance(1,pass,"_unPol"); //fill unpolarized
	
	lambdaTheta = 1; 
	lambdaPhi=0; 
	lambdaThetaPhi=0; 
	w1 = ((3/(4*TMath::Pi()))/(3+lambdaTheta))*(1+(lambdaTheta*cosTheta*cosTheta)+(lambdaPhi*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi*sin2Theta*cosPhi));

	if(print) cout << "Fill Acceptance transverse polarization: " << endl; 
	FillAcceptance(w1,pass, "_transverse"); //transverse polarization
	lambdaTheta=-1;
	w1 = ((3/(4*TMath::Pi()))/(3+lambdaTheta))*(1+(lambdaTheta*cosTheta*cosTheta)+(lambdaPhi*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi*sin2Theta*cosPhi));
	
	
	if(print) cout << "Fill Acceptance longitudinal polarization: " << endl; 
	FillAcceptance(w1,pass,"_longitudinal"); //longitudinal polarization
	
	if(print) cout << "end acceptance. " << endl; 
		
}

void JPsiAnalyzerPAT::fillHistograms(const pat::CompositeCandidate* aCand, const edm::Event& iEvent){
	//cout << "JpsiP_Gen: " << JpsiP_Gen->Perp()<< endl; 
	//cout << "MuP_Gen: " << muPosP_Gen->Perp() << endl; 
	//cout << "MuN_Gen: " << muNegP_Gen->Perp() << endl; 
	
	int iY=-1; 
	int iPT=-1; 
	
	bool print=0; 
	bool opp_charge=false; 
	bool ups_kin=false; 
	bool muon_qual=false; 
	bool seagull=false; 
	bool dimuon_qual = false; 
	bool rapidity_check=false;
	
	bool muon1_good=false;
	bool muon2_good=false;
	
	double vP=aCand->userFloat("vProb");
	bool vP_check=false;
	if(vP>0.01) vP_check=true;
	
	bool trigOK = false;
	if(TriggerBit==1)trigOK=true;
	
	
	if(muPos_qual==1 && muNeg_qual==1) muon_qual=true; 
	
	if(QQ_check==-1) opp_charge=true;
	
	double etaP=TMath::Abs(muPosP_Gen->Eta());
	double etaM=TMath::Abs(muNegP_Gen->Eta());
	
	double PtP=muPosP_Gen->Perp();
	double PtM=muNegP_Gen->Perp();
	
	//kinematic terms 
	double DeltaPhi=TMath::Abs(muPosP_Gen->DeltaPhi(*muNegP_Gen)); 
	double DeltaEta=TMath::Abs(muPosP_Gen->Eta()-muNegP_Gen->Eta());
	double DeltaR=muPosP_Gen->DeltaR(*muNegP_Gen);
	double DeltaPt=TMath::Abs(PtP-PtM);
	double DeltaRPtE=TMath::Sqrt(TMath::Power(0.00157*DeltaPt,2)+TMath::Power(1.2*DeltaPhi,2)+TMath::Power(DeltaEta,2));//deltaR elliptic_deltaPt
	
	hName2D["DeltaPhi_DeltaEta_gen"]->Fill(DeltaEta,DeltaPhi);
	if(fabs(JpsiP_Gen->Rapidity())<=0.6) hName2D["eta1_eta2_y0"]->Fill(etaP,etaM); 
	if(fabs(JpsiP_Gen->Rapidity())>0.6 && fabs(JpsiP_Gen->Rapidity())<1.2) hName2D["eta1_eta2_y1"]->Fill(etaP,etaM); 

	//compare dimuon mass with gen mass
	//fix kinematic cuts
	if(JpsiP_Gen->Perp()>10 && JpsiP_Gen->Perp()<100 && JpsiP_Gen->M()>8.5 && JpsiP_Gen->M() <11.5)
		ups_kin=true; 
	
	if(fabs(JpsiP_Gen->Rapidity())<=0.6)
		rapidity_check=true;
	
	if( ((etaP<1.2 && PtP>4.5) || (etaP>1.2 && etaP<1.4 && PtP>3.5) || (etaP>1.4 && etaP<1.6 && PtP>3.0)) && etaP<1.6 )
		muon1_good=true;
	if( ((etaM<1.2 && PtM>4.5) || (etaM>1.2 && etaM<1.4 && PtM>3.5) || (etaM>1.4 && etaM<1.6 && PtM>3.0)) && etaM<1.6 )
		muon2_good=true;
	
	if(muPosP_Gen->DeltaPhi(*muNegP_Gen)<0)
		seagull=true; 
	
	if(etaP>1.6 && muon1_good) cout << "muon1 "  << etaP << endl;  
	if(etaM>1.6 && muon2_good) cout << "muon2 "  << etaM << endl; 

	
	HX_boost(); //compute boosted parameters
	
	if(seagull && fabs(JpsiP_Gen->Rapidity())<1.2 && ups_kin) acceptance(); 
	
	
	if(muon1_good==0 || muon2_good==0) return;
	if(muPos_qual==-9 || muNeg_qual==-9) return;
	
	double weight= muonEff(muPosP_Gen->Perp(),etaP,"nom")*muonEff(muNegP_Gen->Perp(),etaM,"nom"); //eff1*eff2
	
	if(print) cout << "weight: " << weight << endl; 
	
	double rhoMeasured = hName["rho_dRPtE_TH1"]->Interpolate(DeltaRPtE);
	if(DeltaRPtE>1.7) rhoMeasured=1; 
	
	double dimuon_efficiency = 1./(muonEff(muPosP_Gen->Perp(),etaP,"nom")*muonEff(muNegP_Gen->Perp(),etaM,"nom")*rhoMeasured); 

	
	if(!opp_charge) return; 
	if(!ups_kin) return; 
	
	if(vP_check && seagull && muon_qual && rapidity_check) hName["eff_mass_den"]->Fill(JpsiP->M());
	if(vP_check && seagull && muon_qual && rapidity_check && trigOK) hName["eff_mass_num"]->Fill(JpsiP->M());


	if(seagull && vP_check && fabs(JpsiP->Rapidity())<1.2){
			// fill mass plots 
		double Y=fabs(JpsiP->Rapidity());
		double pt=JpsiP->Perp(); 
		double M=JpsiP->M(); 
		double zeta=1000*JpsiMassErr; 
		
		if(Y<0.6) iY=0; 
		else iY=1; 
		
		iPT=hName["bin_width"]->FindBin(pt)-1; 
		
		TLorentzVector pMP;
		pMP.SetPxPyPzE(muPosP_Gen->Px(),muPosP_Gen->Py(),muPosP_Gen->Pz(),muPosP_Gen->E());
		TLorentzVector pMN;
		pMN.SetPxPyPzE(muNegP_Gen->Px(),muNegP_Gen->Py(),muNegP_Gen->Pz(),muNegP_Gen->E());
		TLorentzVector Mmumu=pMP+pMN;
		
		double mumu_gen=Mmumu.M(); 
		if ((mumu_gen<9.4 || mumu_gen>9.459) && pt>10 && pt<100 && Y<1.2){
		//	cout << "pT: " << pt << " " << Y <<  " " << Form("dm_m_y%d_pt%d",iY,iPT) << " " << Form("m_y%d_pt%d",iY,iPT) << endl; 
			hName2D[Form("dm_m_y%d_pt%d",iY,iPT)]->Fill(M,zeta); 
			hName[Form("m_y%d_pt%d",iY,iPT)]->Fill(M); 
			
			hName2D[Form("dm_mw_y%d_pt%d",iY,iPT)]->Fill(M,zeta,dimuon_efficiency); 
			hName[Form("mw_y%d_pt%d",iY,iPT)]->Fill(M,dimuon_efficiency);
			
		}
		
	}
	
	
	
	if(JpsiP->M()<8.7 || JpsiP->M()>11.3) return;
	
	
	if(vP_check && seagull) hName2D["rho_y_pt_den"]->Fill(JpsiP_Gen->Perp(),fabs(JpsiP_Gen->Rapidity()));
	if(vP_check && trigOK && muon_qual && seagull) hName2D["rho_y_pt_num"]->Fill(JpsiP_Gen->Perp(),fabs(JpsiP_Gen->Rapidity()),1./weight);

	if(!rapidity_check) return; 
	
	if(seagull && trigOK) 
		hName["sg_num"]->Fill(JpsiP_Gen->Perp()); 
	if(trigOK)
		hName["sg_den"]->Fill(JpsiP_Gen->Perp()); 
	
	if(!seagull) return; 
	
	
	if(rapidity_check){
		hName["vertex_den"]->Fill(JpsiP_Gen->Perp());
		if(vP>0.01)hName["vertex_num"]->Fill(JpsiP_Gen->Perp());
	}
	
	if(vP_check) {
		hName["rho_cosTheta_den"]->Fill(cosTheta);
		hName["rho_DM1_den"]->Fill(JpsiDistM1); 
		hName["rho_dRPtE_den"]->Fill(DeltaRPtE);
		hName["rho_pt_den"]->Fill(JpsiP_Gen->Perp()); 
		hName["rho_dR_den"]->Fill(DeltaR); 
		hName2D["rho_Dphi_Deta_den"]->Fill(DeltaEta,DeltaPhi);
		hName2D["rho_dR_pt_den"]->Fill(JpsiP_Gen->Perp(),DeltaR); 
		hName2D["rho_dRPtE_pt_den"]->Fill(JpsiP_Gen->Perp(),DeltaRPtE); 
		hName2D["rho_pt_distM1_den"]->Fill(JpsiP_Gen->Perp(),JpsiDistM1); 
	}

	
	//&& muon_qual
	if(vP_check && trigOK && muon_qual){
		hName["rho_cosTheta_num"]->Fill(cosTheta,1./weight);
		hName["rho_cosTheta_num_uW"]->Fill(cosTheta);
	  
		hName["rho_DM1_num"]->Fill(JpsiDistM1,1./weight);
		
		hName["rho_dRPtE_num"]->Fill(DeltaRPtE,1./weight);
		hName["rho_dRPtE_num_uW"]->Fill(DeltaRPtE);//unweighted for computing uncertainties
		
		hName["rho_pt_num"]->Fill(JpsiP_Gen->Perp(),1./weight); 	
		hName["rho_pt_num_uW"]->Fill(JpsiP_Gen->Perp()); 	

		hName["rho_dR_num"]->Fill(DeltaR,1./weight); 	
		hName2D["rho_Dphi_Deta_num"]->Fill(DeltaEta,DeltaPhi,1./weight);
		hName2D["rho_dR_pt_num"]->Fill(JpsiP_Gen->Perp(),DeltaR,1./weight); 
		hName2D["rho_dRPtE_pt_num"]->Fill(JpsiP_Gen->Perp(),DeltaRPtE,1./weight); 

		hName2D["rho_pt_distM1_num"]->Fill(JpsiP_Gen->Perp(),JpsiDistM1,1./weight); 
		hName["eff_mass_num"]->Fill(JpsiP->M()); 
	}
	
	
}

void
JPsiAnalyzerPAT::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){
	
    //init HLTConfigProvider
    const std::string pro = tagTriggerResults_.process();
    bool changed = true;
	
    //bool init(const edm::Run& iRun, const edm::EventSetup& iSetup, const std::string& processName, bool& changed);
    hltConfigInit_ = false;
    if( hltConfig_.init(iRun, iSetup, pro, changed) ) hltConfigInit_ = true;
	
    prop1_.init(iSetup);
    prop2_.init(iSetup);
	
	//read scale-correction parameters from the database:
	edm::ESHandle<MuScleFitDBobject> dbObject;
	iSetup.get<MuScleFitDBobjectRcd>().get(dbObject);
	corrector_.reset(new MomentumScaleCorrector( dbObject.product() ) );
	//resolutionFunction_.reset(new ResolutionFunction( iConfig.getUntrackedParameter<std::string>("ResolutionsIdentifier") ) );
	resolutionFunction_.reset(new ResolutionFunction( "Resol_JPsi_19pb" ) ); //H: which identifiers are available?
	std::cout << "resolutionFunction_ = " << &*resolutionFunction_ << std::endl;
	
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalyzerPAT::endJob() {
	
	cout << "Total number of events = " << nEvents << endl;
	// cout << "Analyzed runs from  " << runmin << "  to  " << runmax << endl; 
	cout << "============================================================" << endl;
	// cout << "Total number of passed candidates TRIGGER RESULTS ANALYZER   = " << passedTriggerResultsAnalyzer_ << endl;
	cout << "Total number of passed candidates MUON SELECTION CUTS        = " << passedMuonSelectionCuts_ << endl;
	cout << "Total number of passed candidates TRIGGER MATCH              = " << passedTriggerMatch_ << endl;
	
	/* if (_JSON){
	 cout << "JSON file produced" << endl;
	 *JSON << lumitmp <<"]]}";
	 JSON->close();
	 }*/
	
	// Write TTree to File
	if (_writeTree) {
		fOut_->cd();
		tree_->Write();
		
		for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
			hName[it->first]->Write();
		}
		for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++) {
			hName2D[it->first]->Write();
		}
		
		for (std::map<TString,TProfile*>::iterator it=prName.begin(); it!=prName.end(); it++) {
			prName[it->first]->Write();
		}
		fOut_->Close();
	}
	if (_writeDataSet) {
		fOut2_->cd();
		data->Write();
		fOut2_->Close();
	}
}

												  
												  
//! Fill the TTree with all RECO variables
void 
JPsiAnalyzerPAT::fillTreeAndDS(const pat::CompositeCandidate* aCand, const edm::Event& iEvent){
	bool print=0; 
	if(print) std::cout << "FillTreeAndDS" << std::endl;
	const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
	const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));
	
	
	//1.) continue only if we have opposite sign muons <- REMOVED
	//
	if (muon1->charge()*muon2->charge() >= 0) {
		if(muon1->charge() == 0) {
			printf("pat::Muon with zero charge!\n");   return;
		}
		if(muon2->charge() == 0) {
			printf("pat::Muon with zero charge!\n");   return;
		}
		// return;
	}
	
	const pat::Muon *muonPos = 0, *muonNeg = 0;
	if(muon1->charge() > 0){ muonPos = muon1; muonNeg = muon2;}
	else if(muon1->charge() < 0){ muonPos = muon2; muonNeg = muon1;}
	
	muPos_qual=selMuon(muonPos); 
	muNeg_qual=selMuon(muonNeg);
	QQ_check=muonPos->charge()*muonNeg->charge(); 
	
	TrackRef iTrackP = muonPos->innerTrack();
	if(print){
		cout << "Positive Track: " << endl; 
		cout << "pT Error: " << iTrackP->ptError() << endl;
		cout << "eta Error: " << iTrackP->etaError() << endl;
		cout << "phi Error: " << iTrackP->phiError() << endl;
	}	
		
	double muPosErr[]={iTrackP->ptError(),iTrackP->etaError(),iTrackP->phiError()};
	TrackRef iTrackM = muonNeg->innerTrack();
	if(print){
		cout << "Negative track: " << endl; 
		cout << "pT Error: " << iTrackM->ptError() << endl;
		cout << "eta Error: " << iTrackM->etaError() << endl;
		cout << "phi Error: " << iTrackM->phiError() << endl;
	}
	double muNegErr[]={iTrackM->ptError(),iTrackM->etaError(),iTrackM->phiError()};
	
	TVector3 fMuonPosVect, fMuonNegVect;
	
	float theMass = aCand->mass();
	JpsiMassErr = aCand->userFloat("MassErr");
	
	TLorentzVector DiMuon; 

	float theMassErr = 0.;
	if (_MassCorr!=0){
		float CMass = CorrectMass(*muonPos,*muonNeg,_MassCorr, DiMuon, fMuonPosVect, fMuonNegVect);
		if (CMass!=0.0){
			if(print)   cout << "uncorrected mass " << theMass << " corrected mass " << CMass << endl; 
			theMass = CMass;
			theMassErr = JpsiMuScleMassErr;
		}
	}
	if(print) cout << "fill delta Zeta plot: " << endl; 
	double M_u=M_uncertainty(fMuonPosVect,fMuonNegVect,muPosErr,muNegErr); 
	hName["DeltaZeta"]->Fill((JpsiMassErr-M_u)*1000);

	float theRapidity =DiMuon.Rapidity();
	if(print) cout << "Rapidity: " << theRapidity << endl; 
	
	// if (!_useRapidity) theRapidity = theRapidity;
	
	float theCtau; 
	if (_useBS) {theCtau = 10.*aCand->userFloat("ppdlBS");}
	else {theCtau = 10.*aCand->userFloat("ppdlPV");}
	
	float theCtauErr; 
	if (_useBS) {theCtauErr = 10.*aCand->userFloat("ppdlErrBS");}
	else {theCtauErr = 10.*aCand->userFloat("ppdlErrPV");}
	
	dca = aCand->userFloat("DCA");
	
	// MC matching
	reco::GenParticleRef genJpsi = aCand->genParticleRef();
	bool isMatched = (genJpsi.isAvailable() && genJpsi->pdgId() == _oniaPDG);
	
	// Input DataSet Type: P, NP, BG J/psi
	if (isMatched && _isPromptMC) MCType= 0;
	if (isMatched && _isPromptMC == false) MCType=1;
	if (!isMatched && _isMC) MCType=2;
	
	if (isMatched && _removeSignal) return;
	
	reco::GenParticleRef genMu1 = muon1->genParticleRef();
	reco::GenParticleRef genMu2 = muon2->genParticleRef();
	bool isMuMatched = (genMu1.isAvailable() && genMu2.isAvailable() && 
						genMu1->pdgId()*genMu2->pdgId() == -169 && 
						genMu1->momentum().rho() > 2.5 && genMu2->momentum().rho() > 2.5);
	if (isMuMatched && _removeMuons) return;
	if(print) std::cout << "vProb: " << aCand->userFloat("vProb") << std::endl; 
	//store the number of tracks attached to the primary vertex selected by the dimuon:
	if(aCand->hasUserFloat("vertexWeight"))
		vertexWeight = aCand->userFloat("vertexWeight");
	if(aCand->hasUserFloat("sumPTPV"))
		sumPTPV = aCand->userFloat("sumPTPV");
	if(aCand->hasUserInt("countTksOfPV"))
		countTksOfPV = aCand->userInt("countTksOfPV");
	
	if (_writeOutCands) *theTextFile << iEvent.id().run() << "\t" << iEvent.luminosityBlock() << "\t" << iEvent.id().event() << "\t" << theMass << "\n";
	
	if(!_isMC) JpsiP->SetPxPyPzE(DiMuon.Px(),DiMuon.Py(),DiMuon.Pz(),DiMuon.E());
	else JpsiP->SetPxPyPzE(aCand->px(),aCand->py(),aCand->pz(),aCand->energy()); 
	
	//aCand->px(),aCand->py(),aCand->pz(),aCand->energy());
	Jpsict=theCtau;
	JpsictErr=theCtauErr;
	Jpsict_Gen=10.*aCand->userFloat("ppdlTrue");
	JpsiVprob=aCand->userFloat("vProb");
	this->muonStationDistance(aCand);
	Double_t muMass = 0.105658;
	
	if(!_isMC){
		muPosP->SetVectM(fMuonPosVect,muMass);
		muNegP->SetVectM(fMuonNegVect,muMass); 
	}
	else {
		float f_muPosPx, f_muPosPy, f_muPosPz;
		float f_muNegPx, f_muNegPy, f_muNegPz;
		f_muPosPx = muonPos->px();
		f_muPosPy = muonPos->py();
		f_muPosPz = muonPos->pz();
		f_muNegPx = muonNeg->px();
		f_muNegPy = muonNeg->py();
		f_muNegPz = muonNeg->pz();
		Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
		Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
		muPosP->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos); 
		muNegP->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);
	}
	
	if (_writeDataSet) {
		
		bool trigOK = false;
		for (unsigned int iTrig = 0 ; iTrig < _triggerForDataset.size() ; iTrig++) {
			if (mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == 1 ||
				mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == -1 ||
				mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == 2 ) trigOK = true;
		}
		//	isMuonInAccept(muon1) && isMuonInAccept(muon2) &&
		if (theMass > JpsiMassMin && theMass < JpsiMassMax && 
			theCtau > JpsiCtMin && theCtau < JpsiCtMax && 
			aCand->pt() > JpsiPtMin && aCand->pt() < JpsiPtMax && 
			fabs(theRapidity) > JpsiRapMin && fabs(theRapidity) < JpsiRapMax &&
			trigOK) {
			
			//       int ss=999;
			//       if (muon1->charge() + muon2->charge() == 0) ss=0;
			//       if (muon1->charge() + muon2->charge() == 2) ss=1;
			//       if (muon1->charge() + muon2->charge() == -2) ss=2;
			
			Jpsi_Pt->setVal(aCand->pt()); 
			Jpsi_Rap->setVal(theRapidity); 
			Jpsi_MuScleMass->setVal(theMass);
			Jpsi_MuScleMassErr->setVal(theMassErr);
			Jpsi_MassErr->setVal(aCand->userFloat("MassErr"));
			Jpsi_ct->setVal(theCtau);
			Jpsi_ctErr->setVal(theCtauErr);
			// cout << "Type = " << theCat << " pt = " << aCand->pt() << " eta = " << theRapidity << endl;
			// cout << " PPDL = " << theCtau << " Mother = " << aCand->userInt("momPDGId") << " PPDL true = " << 10.*aCand->userFloat("ppdlTrue") << endl;
			Jpsi_MatchType->setIndex((int)isMatched,kTRUE);
			Jpsi_ctTrue->setVal(10.*aCand->userFloat("ppdlTrue"));
			
			Jpsi_PtType->setIndex(getJpsiVarType(aCand->pt(),_ptbinranges),kTRUE);
			Jpsi_RapType->setIndex(getJpsiVarType(fabs(theRapidity),_etabinranges),kTRUE);
			// Fill RooDataSet
			RooArgSet varlist_tmp(*Jpsi_MuScleMass,*Jpsi_MuScleMassErr,*Jpsi_MassErr,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap,*Jpsi_MatchType);   // temporarily remove tag-and-probe weights
			varlist_tmp.add(*Jpsi_ctTrue);   varlist_tmp.add(*Jpsi_PtType);
			varlist_tmp.add(*Jpsi_RapType);  varlist_tmp.add(*Jpsi_ctErr);
			data->add(varlist_tmp);
		}
	}
	
	if(print) cout << "end of fillTreeAndDS: " << endl;
	
}

void JPsiAnalyzerPAT::makeCuts() {
	bool print=0; 
	if(print) std::cout << "makeCuts: "<< std::endl; 
	if(print) std::cout << "collAll.isValid(): " << collAll.isValid() << std::endl; 

	if (collAll.isValid()) {
		
		for(vector<pat::CompositeCandidate>::const_iterator it=collAll->begin();
			it!=collAll->end();++it) {
			if(print) std::cout << "iteration: " << std::endl; 
			NCand++;
			const pat::CompositeCandidate* cand = &(*it);
			// cout << "Now checking candidate of type " << theJpsiCat << " with pt = " << cand->pt() << endl;
			const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
			const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));
			
			// PAT trigger match, 2 muons to the last filter used in the HLT path (new way)
			this->matchMuonToHlt(muon1, muon2);
			muonStationDistance(cand);

			bool trigOK = false;
			for (unsigned int iTrig = 0 ; iTrig < HLTBitNames_.size() ; iTrig++) {
				if (mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == 1 ||
					mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == -1 ||
					mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == 2 ) trigOK = true;
			}
			
			if(trigOK)TriggerBit=1;
			else TriggerBit=0;
			if (requireTriggerMatching_ && !_isMC && !trigOK) continue;

			if(print)std::cout << "PassedCand.push_back: " << std::endl; 
			
			// some counter
			passedTriggerMatch_++;
			
			if (muon1->charge() + muon2->charge() != 0)  continue;
			bool Mu=false;  
			bool DiMu=false; 
			bool Acc=false; 
			
			if (_applycuts && (selMuon(muon1) && selMuon(muon2))) Mu=true; 
				 
			if(selDimuon(cand) ) DiMu=true; 
			if(isMuonInAccept(muon1) && isMuonInAccept(muon2)) Acc=true;
			
			if(!_isMC && (Mu==false || DiMu==false || Acc==false)) continue; 
			
			_thePassedCands.push_back(cand);

			
		}
	}
	
	return;
}

int JPsiAnalyzerPAT::theBestQQ() {
	
	int iBest = -1;
	float maxVProb = -1;
	
	for( unsigned int i = 0; i < _thePassedCands.size(); i++) { 
		if (_thePassedCands.at(i)->userFloat("vProb") > maxVProb) {
			maxVProb = _thePassedCands.at(i)->userFloat("vProb");
			iBest = (int)i;
		}
	}
	
	return iBest;
	
}

bool
JPsiAnalyzerPAT::isMuonInAccept(const pat::Muon* aMuon) {
	// *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
	float eta=fabs(aMuon->eta());
	float pt=fabs(aMuon->pt()); 
	
	bool r1=false;
	bool r2=false;
	bool r3=false; 
	
	if(eta<1.6){
		if(eta>1.4 && eta<1.6 && pt>3.0) r3=true; 
		if(eta>1.2 && eta<1.4 && pt>3.5) r2=true;
		if(eta<1.2 && pt>4.5) r1=true; 
	}
	bool qual=false; 
	if(r1 || r2 || r3) qual=true; 
	
	return qual; 
	
}

bool
JPsiAnalyzerPAT::selMuon(const pat::Muon* aMuon) {
	
	TrackRef iTrack = aMuon->innerTrack();
	const reco::HitPattern& p = iTrack->hitPattern();
	const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
	const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();
	//cout << "pTError: " << iTrack->ptError() << endl;
	TrackRef gTrack = aMuon->globalTrack();
	//   bool globalOK = true;
	//   if (gTrack.isNonnull()) {
	//     const reco::HitPattern& q = gTrack->hitPattern();
	//     globalOK = gTrack->chi2()/gTrack->ndof() < 20.0 && q.numberOfValidMuonHits() > 0 ;
	//   }
	
	bool trackOK = false;
	int NTracks=static_cast<int>(iTrack->found()); 
	
	// cooler way of cutting on tracks
	float fHits = static_cast<double>(iTrack->found()) / static_cast<double> (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
	if (_applyExpHitcuts) {
		trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
		// old way of cutting on tracks  
	} else if(NTracks>10) trackOK=true; 
	
	if(		trackOK &&
	   // 	  globalOK &&
	   iTrack->chi2()/iTrack->ndof() < 1.8 &&
	   aMuon->muonID("TrackerMuonArbitrated") &&
	   aMuon->muonID("TMOneStationTight") &&
	   p.pixelLayersWithMeasurement() > 1 &&
	   fabs(iTrack->dxy(RefVtx)) < 3.0 ) hName["muon_dz"]->Fill(TMath::Abs(iTrack->dz(RefVtx)));
	
	return (// isMuonInAccept(aMuon) &&
			trackOK &&
			// 	  globalOK &&
			iTrack->chi2()/iTrack->ndof() < 1.8 &&
			aMuon->muonID("TrackerMuonArbitrated") &&
			aMuon->muonID("TMOneStationTight") &&
			p.pixelLayersWithMeasurement() > 1 &&
			fabs(iTrack->dxy(RefVtx)) < 3.0 &&
			fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool 
JPsiAnalyzerPAT::selDimuon(const pat::CompositeCandidate* aCand) {
	
	if (!_applyDiMuoncuts) return true;
	const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
	const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));
	
	TVector3 p1,p2;
	
	p1.SetXYZ(muon1->px(), muon1->py(), muon1->pz()); 
	p2.SetXYZ(muon2->px(), muon2->py(), muon2->pz()); 
	
	bool SG=false; 
	
	if(muon1->charge()*p1.DeltaPhi(p2)<0) SG=true; 
	
	return ( aCand->userFloat("vProb") > 0.01 && aCand->mass()>8.4 && aCand->mass()<11.6 && aCand->pt()>8 && aCand->pt()<102 && SG);
}

int 
JPsiAnalyzerPAT::getJpsiVarType(const double jpsivar, vector<double> vectbin) {
	
	for(unsigned int i=0;i<vectbin.size()-1;i++) {
		if(jpsivar > vectbin[i] && jpsivar < vectbin[i+1]) return i+1;
	}
	
	return -999;
}

// reset the global DataSet variables
void
JPsiAnalyzerPAT::resetDSVariables(){
	
    //reset J/psi RECO variables
    JpsiMuScleMassCorr=-9999.;
    JpsiMuScleMassErr=-9999.;
    JpsiMassErr=-9999.;
    sigmaPtPos=-9999.;
    sigmaPtNeg=-9999.;
     Jpsict=-9999.;
    JpsictErr=-9999.;
    Jpsict_Gen=-9999.;
    JpsiVprob=-9999.;
    JpsiDistM1=-9999.;
    JpsiDphiM1=-9999.; 
    JpsiDrM1=-9999.;
    JpsiDistM2=-9999.;
    JpsiDphiM2=-9999.;
    JpsiDrM2=-9999.;
	
	    if(_isMC){
        MCType=-1;
		
        }
	
    //reset EVENT information
	muPos_qual=-9; 
	muNeg_qual=-9; 
	QQ_check=-9;
	TriggerBit=-9;
	NCand=-9;
    eventNb= 0 ;
    runNb= 0 ;
    nPriVtx= 0 ;
    lumiBlock= 0 ;
    vertexWeight = -999.;
    sumPTPV = -999.;
    countTksOfPV = -999;
	
    //reset Trigger Variables
    for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToIntFired_.begin(); clearIt != mapTriggerNameToIntFired_.end(); clearIt++){
        clearIt->second=0;
    }
    for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToPrescaleFac_.begin(); clearIt != mapTriggerNameToPrescaleFac_.end(); clearIt++){
        clearIt->second=-1;
    }
	
    JpsiP->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muPosP->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muNegP->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    JpsiP_Gen->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muPosP_Gen->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muNegP_Gen->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    
}
void JPsiAnalyzerPAT::CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins)
{
	TH1F* h = new TH1F(name, title, nBinsX, xBins);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}

void JPsiAnalyzerPAT::CreateProfile(const char* name,   const char* title,
									const char* xTitle, const char* yTitle,
									Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TProfile* pr = new TProfile(name, title, nBinsX, xLow, xUp);
	
	pr->GetXaxis()->SetTitle(xTitle);
	pr->GetYaxis()->SetTitle(yTitle);
	
	prName[name] = pr;
}

void JPsiAnalyzerPAT::CreateHistogram(const char* name,   const char* title,
									  const char* xTitle, const char* yTitle,
									  Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}
											

void JPsiAnalyzerPAT::CreateHistogram2D(const char* name,   const char* title,
										const char* xTitle, const char* yTitle,
										Int_t nBinsX, Double_t xLow, Double_t xUp,
										Int_t nBinsY,Double_t yLow, Double_t yUp)
{
	TH2F* h = new TH2F(name, title, nBinsX, xLow,xUp,nBinsY, yLow,yUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName2D[name] = h;
}
												  


int JPsiAnalyzerPAT::eta_to_bin(double eta){
	
	//binning for single muon efficiencies
	int bin=-9; 
	if(TMath::Abs(eta)<0.2) bin=0; 
	if(TMath::Abs(eta)>0.2 && TMath::Abs(eta)<0.3) bin=1; 
	if(TMath::Abs(eta)>0.3 && TMath::Abs(eta)<0.6) bin=2;
	if(TMath::Abs(eta)>0.6 && TMath::Abs(eta)<0.8) bin=3;
	if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.0) bin=4; 
	if(TMath::Abs(eta)>1.0 && TMath::Abs(eta)<1.2) bin=5;
	if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<1.4) bin=6; 
	if(TMath::Abs(eta)>1.4 && TMath::Abs(eta)<1.6) bin=7; 
	
	return bin; 
}

double JPsiAnalyzerPAT::muonEff(double pt,double eta,string sys="nom"){
	bool print=0;
	int bin=eta_to_bin(eta); 
	if(bin<0) cout << "eta: " << eta << endl;
	if(bin<0) return 1; 
	
	string name =Form("gEff_MC_PT_AETA%d",bin); 
	//string name= Form("gEff_nominal_AETA%d",bin);//nominal
	string nameP=Form("gEff_totsys_p_AETA%d",bin);
	string nameM=Form("gEff_totsys_m_AETA%d",bin);
	if(print) cout << name << endl; 
	fEff->cd();
	if(fEff->IsOpen()!=1) {
		cout << "Efficiency File not open. " << endl; 
		return 0; 
	}
	
	//TEfficiency *eff = (TEfficiency*)fEff->FindObjectAny(Form("totEff_MCTRUTH_PT_AETA%d",bin));
	
	TGraph *gr=(TGraph*)fEff->Get(name.c_str()); 
	//TGraph *grP=(TGraph*)fEff->Get(nameP.c_str()); 
	//TGraph *grM=(TGraph*)fEff->Get(nameM.c_str()); 

	double emu=gr->Eval(pt)*0.99; 
	//double emu=eff->GetEfficiency(eff->FindFixBin(pt)); 
	
	//double emuP=grP->Eval(pt)/SFp; // positive systematic 
	//double emuM=grM->Eval(pt)/SFm; // negative systematic 
	
	double RV=-9; 
	
	if(sys=="nom") RV=emu; 
	//if(sys=="minus") RV=emuM; 
	//if(sys=="plus") RV=emuP; 
	
	if(print) cout << "Pt: " << pt << "effmu: " << emu << endl; 
	
	fOut_->cd();
	//delete eff; 
	delete gr;
	//delete grP; 
	//delete grM; 
	if(RV<0) cout << "return value <0" << endl; 
	return RV; 
}






//! fill Generator Information
void
JPsiAnalyzerPAT::analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
	bool print=0; 
	if(print) std::cout << "analyzeGenerator: " << std::endl; 
	
    using namespace trigger;
	
    std::vector < const reco::Candidate* > genMuons;
    //bool genjpsi= false;
    reco::Candidate::size_type nrD;
	
    //int count= 0;
    for( size_t i = 0; i < genParticles->size(); ++ i )
    {
        // std::cout << "analyzeGenerator: " << i << std::endl;
        const reco::Candidate & cand = (*genParticles)[ i ];
        int Mc_particleID = cand.pdgId();
        if (abs(Mc_particleID) == _oniaPDG && cand.status()==2 )//&& cand.pt() >= 1)
        {
			//          std::cout << "------::analyzeGenerator:: gen JPsi's: pt=" << cand.pt() << "; eta=" << cand.eta() << "; phi=" << cand.phi() << std::endl;
			
            //Fill global TTree variables
            // JpsiMass_Gen=cand.mass();
            // JpsiPt_Gen=cand.pt();
            // JpsiRap_Gen=cand.rapidity();
            // JpsiPx_Gen=cand.px();
            // JpsiPy_Gen=cand.py();
            // JpsiPz_Gen=cand.pz();
			Double_t enGen = sqrt(cand.px()*cand.px() + cand.py()*cand.py() + cand.pz()*cand.pz() + cand.mass()*cand.mass());
			JpsiP_Gen->SetPxPyPzE(cand.px(),cand.py(),cand.pz(),enGen);
			
            //Jpsict_Gen=10.*cand.userFloat("ppdlTrue"));
			
            nrD= cand.numberOfDaughters();
            int count_muon=0;
            for(reco::Candidate::size_type t=0; t < nrD; t++){
                const reco::Candidate* muon= cand.daughter(t);
                int pID = muon->pdgId();
				//              std::cout << "------::analyzeGenerator:: gen JPsi's daughter pdgId: " << pID << std::endl;
				
                if (abs(pID) == 13 && cand.daughter(t)->status()==1)
                {
                    genMuons.push_back(muon);
					//                  std::cout << "------::analyzeGenerator:: gen JPsi's daughter #: " << count_muon << std::endl;
					//                  std::cout << " muon" << count_muon << " pt=     " << muon->pt() << std::endl;
					//                  std::cout << " muon" << count_muon << " eta=     " << muon->eta() << std::endl;
					//                  std::cout << " moun" << count_muon << " phi=     " << muon->phi() << std::endl;
                    count_muon++;
                }
            }
			
			
            if ( genMuons.empty() ) break;
			
            const reco::Candidate* muon1= genMuons.front();
            const reco::Candidate* muon2= genMuons.back();
			
            // look for opposite charge gen muon pair
            if (muon1->charge()*muon2->charge() <= 0){
                const reco::Candidate *muonPos = 0, *muonNeg = 0;
				
                if(muon1->charge() > 0){ muonPos = muon1; muonNeg = muon2;}
                else if(muon1->charge() < 0){ muonPos = muon2; muonNeg = muon1;}
				
                float f_muPosPx, f_muPosPy, f_muPosPz;
                float f_muNegPx, f_muNegPy, f_muNegPz;
				
                f_muPosPx = muonPos->px();
                f_muPosPy = muonPos->py();
                f_muPosPz = muonPos->pz();
				
                f_muNegPx = muonNeg->px();
                f_muNegPy = muonNeg->py();
                f_muNegPz = muonNeg->pz();
				
                // fill global TTree variables - gen muon
                // muPosPx_Gen=muonPos->px();
                // muPosPy_Gen=muonPos->py();
                // muPosPz_Gen=muonPos->pz();
				
                // muNegPx_Gen=muonNeg->px();
                // muNegPy_Gen=muonNeg->py();
                // muNegPz_Gen=muonNeg->pz();
				
                // fill Polarization variables - gen muons
                Double_t muMass = 0.105658;
				
                Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
                // TLorentzVector *muPos = new TLorentzVector();
                muPosP_Gen->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);
				if(print) std::cout << "mu_px: " << f_muPosPx << " mu_py: " << f_muPosPy << " mu_pz: " << f_muPosPz << std::endl; 
                Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
                // TLorentzVector *muNeg = new TLorentzVector();
                muNegP_Gen->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);
				
                //! Fill Polarization Variables;
                // std::vector< float > thisCosTh, thisPhi;
                // thisCosTh.resize(6); thisPhi.resize(6);
                // this->calcPol(*muPosP_Gen, *muNegP_Gen, thisCosTh, thisPhi);
				
            }
        } // end loop over genParticles
    }
}

void
JPsiAnalyzerPAT::hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup)
{
	
    std::map<std::string, bool> mapTriggernameToTriggerFired;
    std::map<std::string, unsigned int> mapTriggernameToHLTbit;
    // std::map<std::string, unsigned int> mapTriggerNameToPrescaleFac;
	
    for(std::vector<std::string>::const_iterator it= HLTBitNames_.begin(); it !=HLTBitNames_.end(); ++it){
        mapTriggernameToTriggerFired[*it]=false;
        mapTriggernameToHLTbit[*it]=1000;
        // mapTriggerNameToPrescaleFac[*it]=0;
    }
	
    // HLTConfigProvider
    if ( hltConfigInit_ ) {
        
        //! Use HLTConfigProvider
		const unsigned int n= hltConfig_.size();
		for (std::map<std::string, unsigned int>::iterator it = mapTriggernameToHLTbit.begin(); it != mapTriggernameToHLTbit.end(); it++) {
			unsigned int triggerIndex= hltConfig_.triggerIndex( it->first );
			if (triggerIndex >= n) {
				//std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " not available in config!" << std::endl;
			}
			else {
				it->second= triggerIndex;
				//std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " available in config!" << std::endl;
			}
		}
    }
    
    // Get Trigger Results
    try {
		iEvent.getByLabel( tagTriggerResults_, handleTriggerResults_ );
		//cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResult is present in current event" << endl;
    }
    catch(...) {
		//cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults NOT present in current event" << endl;
    }
    if ( handleTriggerResults_.isValid() ){
		//cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults IS valid in current event" << endl;
		
		// loop over Trigger Results to check if paths was fired
		for(std::vector< std::string >::iterator itHLTNames= HLTBitNames_.begin(); itHLTNames != HLTBitNames_.end(); itHLTNames++){
			const std::string triggerPathName =  *itHLTNames;
			//std::cout << "[FloJPsiAnalyzer::hltReport] --- TriggerName --- TriggerName LOOP" << std::endl;
			
			if ( mapTriggernameToHLTbit[triggerPathName] < 1000 ) {
				if (handleTriggerResults_->accept( mapTriggernameToHLTbit[triggerPathName] ) ){
					mapTriggerNameToIntFired_[triggerPathName] = 3;
				}
				
				//-------prescale factor------------
				if (!_isMC) {
					const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerPathName));
					//std::cout << "[FloJPsiAnalyzer::prescalvalues] --- TriggerName"<<triggerPathName<<" prescales first "<< prescales.first <<" prescales second "<< prescales.second <<std::endl;
					mapTriggerNameToPrescaleFac_[triggerPathName] = prescales.first * prescales.second;
				}
			}
		}
    } else cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerResults NOT valid in current event" << endl;
}

void
JPsiAnalyzerPAT::matchMuonToHlt(const pat::Muon* muon1, const pat::Muon* muon2)
{
	
    std::string HLTL3MuCollName = "hltL3MuonCandidates::" + tagTriggerResults_.process();
    std::string HLTL2MuCollName = "hltL2MuonCandidates::" + tagTriggerResults_.process();
    std::string HLTTrackCollName = "hltMuTrackJpsiCtfTrackCands::" + tagTriggerResults_.process();
    std::string HLTTkMuCollName = "hltMuTkMuJpsiTrackerMuonCands::" + tagTriggerResults_.process();
    
    //! Loop over Trigger Paths and match muons to last Filter/collection
    for ( std::map<std::string, int>::iterator it = mapTriggerNameToIntFired_.begin(); it != mapTriggerNameToIntFired_.end(); it ++ ) {
		
        std::string triggerName = it->first;
		
        //! just use Triggers which are in TriggerResults; value == 3
        if ( it->second != 3 ) continue;
		
        std::string hltLastFilterName = mapTriggerToLastFilter_[triggerName];
		
        const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter( hltLastFilterName );
        const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter( hltLastFilterName );
        bool pass1 = mu1HLTMatches.size() > 0;
        bool pass2 = mu2HLTMatches.size() > 0; 
		
        // treat "MuX_TrackX" Trigger separately: Match by Tracker collection: hltMuTrackJpsiCtfTrackCands
		std::vector<std::string>::iterator theCheck = std::find(HLTBitNames_MuTrack.begin(),HLTBitNames_MuTrack.end(),triggerName);
		if (theCheck != HLTBitNames_MuTrack.end()) {
			bool matchedMu3[2] = {false, false}, matchedTrack[2] = {false, false};
			for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
				if (mu1HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[0] = true;	     
				if (mu1HLTMatches[k].collection() == HLTTrackCollName ) matchedTrack[0] = true;
			}
			for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
				if (mu2HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[1] = true;
				if (mu2HLTMatches[k].collection() == HLTTrackCollName ) matchedTrack[1] = true;
			}
			if( matchedMu3[0] && matchedTrack[1] )
				mapTriggerNameToIntFired_[triggerName] = 1;
			else if( matchedMu3[1] && matchedTrack[0] )
				mapTriggerNameToIntFired_[triggerName] = -1;
			if( matchedMu3[0] && matchedTrack[1] && matchedMu3[1] && matchedTrack[0] )
				mapTriggerNameToIntFired_[triggerName] = 2;
        }
		
		//         // treat "MuX_TkMuX" Trigger separately: Match by Tracker collection:hltMuTkMuJpsiTrackerMuonCands
		theCheck = std::find(HLTBitNames_MuTkMu.begin(),HLTBitNames_MuTkMu.end(),triggerName);
		if (theCheck != HLTBitNames_MuTkMu.end()) {
			
			bool matchedMu3[2] = {false, false}, matchedTrack[2] = {false, false};
			for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
				if (mu1HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[0] = true;
				if (mu1HLTMatches[k].collection() == HLTTkMuCollName ) matchedTrack[0] = true;
			}
			for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
				if (mu2HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[1] = true;
				if (mu2HLTMatches[k].collection() == HLTTkMuCollName ) matchedTrack[1] = true;
			}
			if( matchedMu3[0] && matchedTrack[1] )
				mapTriggerNameToIntFired_[triggerName] = 1;
			else if( matchedMu3[1] && matchedTrack[0] )
				mapTriggerNameToIntFired_[triggerName] = -1;
			if( matchedMu3[0] && matchedTrack[1] && matchedMu3[1] && matchedTrack[0] )
				mapTriggerNameToIntFired_[triggerName] = 2;
        }
		
		// treat "MuX_L2MuX" Trigger separately: Match by L2 collection: hltL2MuonCandidates and on a different SaveTag'ed filter
        theCheck = std::find(HLTBitNames_MuL2Mu.begin(),HLTBitNames_MuL2Mu.end(),triggerName);
		if (theCheck != HLTBitNames_MuL2Mu.end()) {
			
	        std::string triggerNameSpecial = triggerName + "_special";
			std::string hltLastFilterName2 = mapTriggerToLastFilter_[triggerNameSpecial];
			
			const pat::TriggerObjectStandAloneCollection mu1Level2Matches = muon1->triggerObjectMatchesByFilter( hltLastFilterName2 );
			const pat::TriggerObjectStandAloneCollection mu2Level2Matches = muon2->triggerObjectMatchesByFilter( hltLastFilterName2 );
			
			bool matchedMu3[2] = {false, false}, matchedMu2[2] = {false, false};
			for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
				if (mu1HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[0] = true;
			}
			for (unsigned k = 0; k < mu1Level2Matches.size(); ++k) {
				if (mu1Level2Matches[k].collection() == HLTL2MuCollName ) matchedMu2[0] = true;                }
			for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
				if (mu2HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[1] = true;
			}
			for (unsigned k = 0; k < mu2Level2Matches.size(); ++k) {
				if (mu2Level2Matches[k].collection() == HLTL2MuCollName ) matchedMu2[1] = true;
			}
			if( matchedMu3[0] && matchedMu2[1] )
				mapTriggerNameToIntFired_[triggerName] = 1;
			else if( matchedMu3[1] && matchedMu2[0] )
				mapTriggerNameToIntFired_[triggerName] = -1;
			if( matchedMu3[0] && matchedMu2[1] && matchedMu3[1] && matchedMu2[0] )
				mapTriggerNameToIntFired_[triggerName] = 2;
        }
		
        // All the other Paths match by last filter:
		theCheck = std::find(HLTBitNames_DoubleMu.begin(),HLTBitNames_DoubleMu.end(),triggerName);
		if (theCheck != HLTBitNames_DoubleMu.end() && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
		theCheck = std::find(HLTBitNames_SingleMu.begin(),HLTBitNames_SingleMu.end(),triggerName);
		if (theCheck != HLTBitNames_SingleMu.end() && (pass1 == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
    }
}

void 
JPsiAnalyzerPAT::muonStationDistance (const pat::CompositeCandidate* aCand) 
{
	
	const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
	const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));
	const reco::Candidate &d1 = *muon1; 
	const reco::Candidate &d2 = *muon2;
	const reco::RecoCandidate *mu1 = dynamic_cast<const reco::RecoCandidate *>(&d1);
	const reco::RecoCandidate *mu2 = dynamic_cast<const reco::RecoCandidate *>(&d2);
    
	// Propagate to station 1
	TrajectoryStateOnSurface prop1_Stat1 = prop1_.extrapolate(*mu1);
	TrajectoryStateOnSurface prop2_Stat1 = prop1_.extrapolate(*mu2);
	if (prop1_Stat1.isValid() && prop2_Stat1.isValid()) {
		JpsiDphiM1 = deltaPhi<float>(prop1_Stat1.globalPosition().phi(), prop2_Stat1.globalPosition().phi());
		JpsiDrM1   = hypot(JpsiDphiM1, std::abs<float>(prop1_Stat1.globalPosition().eta() - prop2_Stat1.globalPosition().eta()));
		JpsiDistM1 = (prop1_Stat1.globalPosition()-prop2_Stat1.globalPosition()).mag();
	}
	// Propagate to station 2
	TrajectoryStateOnSurface prop1_Stat2 = prop2_.extrapolate(*mu1);
	TrajectoryStateOnSurface prop2_Stat2 = prop2_.extrapolate(*mu2);
	if (prop1_Stat2.isValid() && prop2_Stat2.isValid()) {
		JpsiDphiM2 = deltaPhi<float>(prop1_Stat2.globalPosition().phi(), prop2_Stat2.globalPosition().phi());
		JpsiDrM2   = hypot(JpsiDphiM2, std::abs<float>(prop1_Stat2.globalPosition().eta() - prop2_Stat2.globalPosition().eta()));
		JpsiDistM2 = (prop1_Stat2.globalPosition()-prop2_Stat2.globalPosition()).mag();
	}
	
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalyzerPAT);
