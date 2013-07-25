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
//$Id: JPsiAnalyzerPAT.cc,v 1.1 2013/04/17 03:38:09 bcarlson Exp $
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
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>
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
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void makeCuts() ;
      int theBestQQ();
      void fillTreeAndDS(const pat::CompositeCandidate* aCand, const edm::Event&);
      bool isMuonInAccept(const pat::Muon* aMuon);
      bool selMuon(const pat::Muon* aMuon);
      bool selDimuon(const pat::CompositeCandidate* aCand);
      int getJpsiVarType(const double jpsivar, vector<double> vectbin);
      double CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode);

	  int eta_to_bin(double eta);
	  double muonEff(double pt,double eta,string sys); // efficiency stuff 
	
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
	  TFile* fEff; // file with single muon efficiencies

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
	
	  unsigned int Ngen; 
	  unsigned int N_triggered; 
	
	  unsigned int HLT_trig; 
	  unsigned int HLT_tot; 
	
	  unsigned int Ndenom; 
	  unsigned int Nnum;
	
	  bool pass_cuts;

	  //Histograms
	  TH1D *trigger_pass; 
	  TH1D *trigger_passP;
 	  TH1D *trigger_passM;
	  TH1D *kin_pass;
	  TH1D *gen_mass; 
	  TH1D *trigg_bits; 
	  TH1D *diff_p; 
	  TH1D *diff_n; 
	
	  TH2D *difftheta_diffphi; 	 
	
      TH1D *DR_fail; 
	  TH1D *DR_pass; 
	
	  TH1D *deltapt_fail; 
	  TH1D *deltapt_match; 
	
	  TH1D *iTrack_hist;
	TH1F *fHits_hist;
	
	TH1D *chi2_hist;
	TH1D *tma_hist;
	TH1D *tmt_hist;
	TH1D *tpm_hist;
	
	TH1D *SG_num;
	TH1D *SG_den; 
	
	
	TH1D *dm;
	TH1D *Mass; 
	
	TH2D *dm_m; 
	
      // limits 
      float JpsiMassMin;
      float JpsiMassMax;
      float JpsiCtMin;
      float JpsiCtMax;
      float JpsiPtMin;           // SET BY 
      float JpsiPtMax;           // DEFINITION
      float JpsiRapMin;          // OF BIN
      float JpsiRapMax;          // LIMITS

	  //event kinematics
	  double candpt;
	  double candy;
	
	
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

   Ngen=0; 
   N_triggered=0; 	
	
	Ndenom=0; 
	Nnum=0; 
	
	HLT_trig=0; 
	HLT_tot=0; 
	
	pass_cuts=false; 
	
  JpsiMassMin = _massMin;
  JpsiMassMax = _massMax;
  JpsiCtMin = -2.0;
  JpsiCtMax = 3.5;
	
	candpt=-9;
	candy=-99;

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
void
JPsiAnalyzerPAT::beginJob()
{

  //read scale-correction parameters from the database:
  // edm::ESHandle<MuScleFitDBobject> dbObject;
  // iSetup.get<MuScleFitDBobjectRcd>().get(dbObject);
  // corrector_.reset(new MomentumScaleCorrector( dbObject.product() ) );

	
	//fEff = new TFile("/afs/hephy.at/user/k/knuenz/public/forBenjamin/singleMuTruthEff_18Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins200MeV.root","READ"); //eff_data.root MC Efficiency eff_all uses data scale
	fEff=new TFile("MCEff2.root","READ"); 
    std::cout << "[JPsiAnalyzerPAT] --- beginJob " << std::endl;
  if (_writeTree) {
    fOut_ = new TFile(_treefilename.c_str(), "RECREATE");
    fOut_->cd();

	//Create Histos
	  
	trigger_pass = new TH1D("trigger_pass", "Passed Triggers: passedTriggerMatch_", 9,10,100);
	kin_pass = new TH1D("kin_pass","Passed Muon: passedMuonSelectionCuts_", 9,10,100); 
	trigger_passP = new TH1D("trigger_passP", "Passed Triggers: passedTriggerMatch_, sys+", 9,10,100);
	trigger_passM = new TH1D("trigger_passM", "Passed Triggers: passedTriggerMatch_, sys-", 9,10,100);

	gen_mass = new TH1D("gen_mass", "Generator Level Mass Distribution", 100, 8.7,11.2);   
	trigg_bits = new TH1D("trigg_bits", "Trigger Bits", 100, 0, 5); 
	diff_p = new TH1D("diff_p", "|Muon-Gen| P_{T}, MeV",100, 0, 100);   
	diff_n = new TH1D("diff_n", "Muon-Gen| P_{T}, MeV", 100, 0, 100); 
	deltapt_fail = new TH1D("deltapt_fail", "#Delta P_{T} Fail", 100, -100, 100);   
	deltapt_match = new TH1D("deltapt_match", "#Delta P_{T} Match",100, -100, 100);   
	DR_fail = new TH1D("DR_fail", "DR Fail", 100, 0, 5);   
	DR_pass = new TH1D("DR_pass", "DR Pass", 100, 0, 5);   
	  iTrack_hist = new TH1D("iTrack_hist", "iTrack Found", 20,0,50);   
	  fHits_hist = new TH1F("fHits_hits","fHits Hist", 100,0,1); 
	  
	  difftheta_diffphi = new TH2D("difftheta_diffphi", "#Delta#theta vs. #Delta#phi; #Delta#phi; #Delta#theta",100, -3.14, 3.13, 100, -3.14, 3.14);   
	  
	  chi2_hist = new TH1D("chi2_hist","#Chi^{2}/NDOF",100,0,10); 
	  tma_hist = new TH1D("tma_hist", "TrackerMuonArbitrated",3,-1.5,1.5); 
	  tmt_hist = new TH1D("tmt_hist","TMOneStationTight",3,-1.5,1.5); 
	  tpm_hist = new TH1D("tpm_hist","pixel Layers With Measurement >1", 10,0,10); 
	  
	  SG_num = new TH1D("SG_num","Seagull Efficiency Num", 10,10,100); 
	  SG_den = new TH1D("SG_den","Seagull Efficiency Den", 10,10,100); 
	  
	  dm = new TH1D("dm","#zeta MeV, 10<P_{T}<50; #zeta MeV; Events",800,0,200);
	  Mass = new TH1D("Mass","M_{#mu#mu} GeV",250 ,8.7,11.2);
	  dm_m = new TH2D("dm_m", "#zeta vs M_{#mu#mu};M_{#mu#mu} GeV; #zeta MeV", 250,8.7,11.2,800,0,200); 
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
        // tree_->Branch("JpsiMass_Gen",   &JpsiMass_Gen,  "JpsiMass_Gen/D");
        // tree_->Branch("JpsiPt_Gen",     &JpsiPt_Gen,    "JpsiPt_Gen/D");
        // tree_->Branch("JpsiRap_Gen",    &JpsiRap_Gen,   "JpsiRap_Gen/D");
        // tree_->Branch("JpsiPx_Gen",     &JpsiPx_Gen,    "JpsiPx_Gen/D");
        // tree_->Branch("JpsiPy_Gen",     &JpsiPy_Gen,    "JpsiPy_Gen/D");
        // tree_->Branch("JpsiPz_Gen",     &JpsiPz_Gen,    "JpsiPz_Gen/D");
        tree_->Branch("Jpsict_Gen",     &Jpsict_Gen,    "Jpsict_Gen/D");
        tree_->Branch("muPosP_Gen",  "TLorentzVector", &muPosP_Gen);
        tree_->Branch("muNegP_Gen",  "TLorentzVector", &muNegP_Gen);
        // tree_->Branch("muPosPx_Gen",    &muPosPx_Gen,   "muPosPx_Gen/D");
        // tree_->Branch("muPosPy_Gen",    &muPosPy_Gen,   "muPosPy_Gen/D");
        // tree_->Branch("muPosPz_Gen",    &muPosPz_Gen,   "muPosPz_Gen/D");
        // tree_->Branch("muNegPx_Gen",    &muNegPx_Gen,   "muNegPx_Gen/D");
        // tree_->Branch("muNegPy_Gen",    &muNegPy_Gen,   "muNegPy_Gen/D");
        // tree_->Branch("muNegPz_Gen",    &muNegPz_Gen,   "muNegPz_Gen/D");
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


double JPsiAnalyzerPAT::CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode){  

	std::cout<<"JPsiAnlyazerPAT::CorrectMass. " << endl; 
  //MuScle Fit corrections
  //1) correct the momentum scale:
  double corrPt1 = (*corrector_)(mu1);
  double corrPt2 = (*corrector_)(mu2);
   cout << "original pT1 " << (mu1.innerTrack()->momentum()).Rho() << " corrected pT1 " << corrPt1 << endl;
   cout << "original pT2 " << (mu2.innerTrack()->momentum()).Rho() << " corrected pT2 " << corrPt2 << endl;

  const double mumass = 0.105658;
  TLorentzVector mu1Corr, mu2Corr; 
  mu1Corr.SetPtEtaPhiM(corrPt1, mu1.innerTrack()->eta(), mu1.innerTrack()->phi(), mumass);
  mu2Corr.SetPtEtaPhiM(corrPt2, mu2.innerTrack()->eta(), mu2.innerTrack()->phi(), mumass);
  TLorentzVector onia = mu1Corr+mu2Corr;
  JpsiMuScleMassCorr = onia.M();

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
	if(print) std::cout << "JPsiAnalyzerPAT::analyze" << endl; 
   nEvents++;
	if(print) cout << "event: " << nEvents << endl; 	
	candpt=-9; 
   // reset TTree Variables
   if (_writeTree) this->resetDSVariables();

   // check HLT TriggerReuslts
   this->hltReport(iEvent, iSetup);

   bool trigOK = false;
   for (unsigned int iTrig = 0 ; iTrig < HLTBitNames_.size() ; iTrig++) {
	   if(print) cout << "HLTBinNames: " <<  HLTBitNames_.at(iTrig) << endl; 
	   if(print) cout << "Trig Bit: " << mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] << endl;
     if (mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == 3) trigOK = true;
   }
   if(print) cout << "requireTriggerMatching_: " << requireTriggerMatching_ << endl; 
   if(print) cout << "trigOK : " << trigOK << endl; 
   if(print) cout << "storeAllMCEvents: " << _storeAllMCEvents << endl; 
   
	//if _storeAllMCEvents is set to true in the config file, then the function will never end here 	
   if (requireTriggerMatching_ && !trigOK && !_storeAllMCEvents) return;

   // Event related infos
   eventNb= iEvent.id().event() ;
   runNb=iEvent.id().run() ;
   lumiBlock= iEvent.luminosityBlock() ;

   Handle<reco::VertexCollection> privtxs;
   iEvent.getByLabel("offlinePrimaryVertices", privtxs);
   nPriVtx = privtxs->size();
   VertexCollection::const_iterator privtx;

	//Get gen particle info 
	pass_cuts=false; 
	Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByLabel( _genParticles, genParticles );
	if ( genParticles.isValid() )
	{
		if(print) std::cout << "------ analyze GENERATED Upsilons:" << std::endl;
		this->analyzeGenerator( genParticles );
		if(pass_cuts==0) return; // if the events fail the cuts, don't do anything... 
	}
	//if the events do pass the cuts, continue on... and fill histograms etc 
	Ngen++;
	
	
   if ( privtxs->begin() != privtxs->end() ) {
     privtx=privtxs->begin();
     RefVtx = privtx->position();
   } else {
     RefVtx.SetXYZ(0.,0.,0.);
   }
   // }

   try {iEvent.getByLabel(_patJpsi,collAll);} 
   catch (...) {cout << "Upsilon not present in event!" << endl;}

   _thePassedCands.clear();

	//makeCuts fills the candidates vector with any good candidates 
   // APPLY CUTS
   this->makeCuts();
	bool one_fail=false; 
	for (unsigned int iTrig = 0 ; iTrig < HLTBitNames_.size() ; iTrig++) {
		trigg_bits->Fill(mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)]); 
		if(mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)]!=1) one_fail=true; 
		if(one_fail){
			//cout << "HLTBinNames: " <<  HLTBitNames_.at(iTrig) << endl; 
			//cout << "Trig Bit: " << mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] << endl;
			
		}
		
	}
	
	if(print){
		//one_fail && _thePassedCands.size()==0
		cout << "Upsilon: " << endl;
		JpsiP_Gen->Print(); 
		cout << endl; 
		cout << "Muon+: " << endl;
		muPosP_Gen->Print();
		cout << endl; 
		cout << "Muon-: " << endl;
		muNegP_Gen->Print();
		cout << endl;
	}
	
	if(print && one_fail && _thePassedCands.size()>0) cout << "Passed cuts. " << endl; 
	
   bool storeEvent = false;
   if(print) cout << "_onlythebest: " << _onlythebest << endl; 
   if(print) cout << "_thePassedCands.size(): " << _thePassedCands.size() << endl; 	
	// BEST J/PSI? 
	if (_onlythebest && _thePassedCands.size()>0) {  // yes, fill simply the best
		if(print) cout << "Fill the best, with at least 1 candidate. " << endl; 
		int iBest = theBestQQ();
		if(print) cout << "iBest: " << iBest <<endl; 
		if (iBest > -1){
			if(print) cout << "iBest >-1 loop. " << endl; 
			fillTreeAndDS(_thePassedCands.at(iBest), iEvent);
			passedMuonSelectionCuts_++;
			N_triggered++;
			// std::cout << "Before filling muon: " << candpt << endl; 
			if(candpt<0) cout << "Error in pt. " << endl; 
			storeEvent=true;
		}
	} else {   // no, fill all candidates passing cuts (possibly wrong-sign)
		//loop over all candidates, instead of just filling one. 
		if(print) cout << "Filling all " << _thePassedCands.size() << " candidates. " << endl; 
		for( unsigned int count = 0; count < _thePassedCands.size(); count++) { 
			if(print ) cout << "Filled candidate that was not the best. " << endl; 
			fillTreeAndDS(_thePassedCands.at(count),iEvent);
			passedMuonSelectionCuts_++;
			// std::cout << "before filling muon things. " << candpt << endl; 
			if(candpt<0) cout << "Error in pt. " << endl; 
			
		}
   }
  if(print) cout << "_storeAllMCEvents: " << _storeAllMCEvents << endl; 
  if(print) cout << "storeEvent: " << storeEvent << endl; 
  
	//! FILL GENERATOR COLLECTION and store the event
	if ( _storeAllMCEvents || storeEvent ) {
		//either we store all the events, or we just store the events "best" candidates, but either way we store events 
		if (_writeTree) tree_->Fill();
		
		// Write all Branches to the Tree ONLY 
		// - for the best candidate
		// - for the opposite sign
	}
	if(print) cout << endl << endl; 
}

void
JPsiAnalyzerPAT::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

	std::cout << "beginRun: " << endl; 
    //init HLTConfigProvider
    const std::string pro = tagTriggerResults_.process();
    bool changed = true;
	
	cout << "Process: " << pro << endl;
	
	
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
  std::cout << "endJob." << endl; 
  cout << "Total number of events = " << nEvents << endl;
  // cout << "Analyzed runs from  " << runmin << "  to  " << runmax << endl; 
  cout << "============================================================" << endl;
  // cout << "Total number of passed candidates TRIGGER RESULTS ANALYZER   = " << passedTriggerResultsAnalyzer_ << endl;
  cout << "Total number of passed candidates MUON SELECTION CUTS        = " << passedMuonSelectionCuts_ << endl;
  cout << "Total number of passed candidates TRIGGER MATCH              = " << passedTriggerMatch_ << endl;
  
	cout << "My N_Triggered: " << N_triggered << endl; 
	cout << "My N_gen: " << Ngen << endl; 
	
	cout << "Nnum: " << Nnum << endl; 
	cout << "Ndenom: " << Ndenom << endl; 
	
	cout << "rho(Nnum/Ndenom): " << static_cast<double>(Nnum)/static_cast<double> (Ndenom) << endl; 
	
  cout << "Rho (?): " << static_cast<double> (N_triggered)/static_cast<double>(Ngen) << endl; 
	
	cout << "HLT_triggered: " << HLT_trig << endl; 
	cout << "HLT_tot: " << HLT_tot << endl; 
	
	trigger_pass->Sumw2();
	kin_pass->Sumw2();
	
	for(int i=1; i<kin_pass->GetNbinsX(); i++){
		kin_pass->SetBinError(i,0); 
		trigger_passP->SetBinError(i,0);
		trigger_passM->SetBinError(i,0);
	}
	
	TH1D *rho_pt=(TH1D*)trigger_pass->Clone();
	TH1D *rho_ptP=(TH1D*)trigger_passP->Clone();
	TH1D *rho_ptM=(TH1D*)trigger_passM->Clone();

	rho_pt->SetName("rho_pt"); 
	rho_ptP->SetName("rho_ptP"); 
	rho_ptM->SetName("rho_ptM"); 
	rho_pt->SetTitle("#rho=#epsilon_{#mu#mu}/#epsilon_{#mu}*#epsilon_{#mu}"); 
	
	rho_pt->Divide(kin_pass); 
	rho_ptP->Divide(kin_pass); 
	rho_ptM->Divide(kin_pass); 
		
	for (int i=1; i<=rho_pt->GetNbinsX(); i++) {
		double rho_i=rho_pt->GetBinContent(i);
		double trigE_i=trigger_pass->GetBinError(i);
		double trig_i=trigger_pass->GetBinContent(i); 
		double passE_i=kin_pass->GetBinError(i);
		double pass_i=kin_pass->GetBinContent(i); 
		double err_i=rho_i*TMath::Sqrt(TMath::Power((trigE_i/trig_i),2)+TMath::Power((passE_i/pass_i),2));
		//rho_pt->SetBinError(i,err_i);
	}

	
  /* if (_JSON){
    cout << "JSON file produced" << endl;
    *JSON << lumitmp <<"]]}";
    JSON->close();
    }*/

  // Write TTree to File
  if (_writeTree) {
    fOut_->cd();
	
	  chi2_hist->Write();
	  tma_hist->Write();
	  tmt_hist->Write();
	  tpm_hist->Write(); 
	  dm->Write();
	  Mass->Write();
	  dm_m->Write();
	  
	  iTrack_hist->Write();
	  fHits_hist->Write();
	kin_pass->Write();
	rho_pt->Write();   
	rho_ptM->Write();
    rho_ptP->Write();
	trigger_pass->Write(); 
	gen_mass->Write();
	trigg_bits->Write();
	  diff_p->Write(); 
	  diff_n->Write();
	  deltapt_fail->Write();
	  difftheta_diffphi->Write();
	  deltapt_match->Write();
	  DR_fail->Write();
	  DR_pass->Write();
	  
	  SG_den->Write(); 
	  SG_num->Write(); 
	  
    tree_->Write();
    fOut_->Close();
  }
  if (_writeDataSet) {
    fOut2_->cd();
    data->Write();
    fOut2_->Close();
  }
}

//! Fill the TTree with all RECO variables

/*
double JPsiAnalyzerPAT::compute_dm(TVector3 fMuon1Vect, TVector3 fMuon2Vect){
	bool print =1; 
	if(fMuon1Vect.X()==fMuon1Vect.Y()==fMuon1Vect.Z())
		cout << "No data for muon1." << endl; 
	if(fMuon2Vect.X()==fMuon2Vect.Y()==fMuon2Vect.Z())
		cout <<"No data for muon2."<<endl; 
	
	double pt1 = fMuon1Vect.Perp();
	double pt2 = fMuon2Vect.Perp();
	
	double Pl1 = fMuon1Vect.Pz();
	double Pl2 = fMuon2Vect.Pz();
	
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
	if(print) cout <<"M: " << M << " sigmaM: " << sigmaM << endl; 
	
	return sigmaM;
	
}
*/

void 
JPsiAnalyzerPAT::fillTreeAndDS(const pat::CompositeCandidate* aCand, const edm::Event& iEvent){
 //	std::cout << "fillTreeAndDS: " << endl;
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

  //
  float theMass = aCand->mass();
  JpsiMassErr = aCand->userFloat("MassErr");

  float theMassErr = 0.;
  if (_MassCorr!=0){
    float CMass = CorrectMass(*muon1,*muon2,_MassCorr);
    if (CMass!=0.0){
      //      cout << "uncorrected mass " << theMass << " corrected mass " << CMass << endl; 
      theMass = CMass;
      theMassErr = JpsiMuScleMassErr;
    }
  }

  float theRapidity = aCand->rapidity();
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

	if(isMuMatched==0) cout << "isMuMatched: " << isMuMatched << endl;
	if(isMatched==0) cout << "isMatched: " << isMatched << endl;
	
	
  //store the number of tracks attached to the primary vertex selected by the dimuon:
  if(aCand->hasUserFloat("vertexWeight"))
    vertexWeight = aCand->userFloat("vertexWeight");
  if(aCand->hasUserFloat("sumPTPV"))
    sumPTPV = aCand->userFloat("sumPTPV");
  if(aCand->hasUserInt("countTksOfPV"))
     countTksOfPV = aCand->userInt("countTksOfPV");

  if (_writeOutCands) *theTextFile << iEvent.id().run() << "\t" << iEvent.luminosityBlock() << "\t" << iEvent.id().event() << "\t" << theMass << "\n";

  // write out JPsi RECO information
  // JpsiPt=aCand->pt();
  // JpsiRap=theRapidity;
//   JpsiCharge=theCharge;
  // std::cout << "[JPsiAnalyzerPAT::fillTreeAndDS] ----- JpsiCharge: " << theCharge << std::endl;
  // JpsiPx=aCand->px();
  // JpsiPy=aCand->py();
  // JpsiPz=aCand->pz();
  JpsiP->SetPxPyPzE(aCand->px(),aCand->py(),aCand->pz(),aCand->energy());
  Jpsict=theCtau;
  JpsictErr=theCtauErr;
  Jpsict_Gen=10.*aCand->userFloat("ppdlTrue");
  JpsiVprob=aCand->userFloat("vProb");
  this->muonStationDistance(aCand);
  
  // write out Muon RECO information
  float f_muPosPx, f_muPosPy, f_muPosPz;
  float f_muNegPx, f_muNegPy, f_muNegPz;
  f_muPosPx = muonPos->px();
  f_muPosPy = muonPos->py();
  f_muPosPz = muonPos->pz();
  f_muNegPx = muonNeg->px();
  f_muNegPy = muonNeg->py();
  f_muNegPz = muonNeg->pz();
  // muPosPx= f_muPosPx ;
  // muPosPy= f_muPosPy ;
  // muPosPz= f_muPosPz ;
  // muNegPx= f_muNegPx ;
  // muNegPy= f_muNegPy ;
  // muNegPz= f_muNegPz ;
  
  //write out Calculated Polarization variables
  Double_t muMass = 0.105658;
  
  Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
  // TLorentzVector *muPosP = new TLorentzVector();
  muPosP->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);
  
  Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
  // TLorentzVector *muNegP = new TLorentzVector();
  muNegP->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);
  
  //! Fill Polarization Variables;
  // std::vector< float > thisCosTh, thisPhi;
  // thisCosTh.resize(6); thisPhi.resize(6);
  // this->calcPol(*muPosP, *muNegP, thisCosTh, thisPhi);

  if (_writeDataSet) {
    
    bool trigOK = false;
    for (unsigned int iTrig = 0 ; iTrig < _triggerForDataset.size() ; iTrig++) {
      if (mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == 1 ||
	  mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == -1 ||
	  mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == 2 ) trigOK = true;
    }
	  
    if (theMass > JpsiMassMin && theMass < JpsiMassMax && 
	theCtau > JpsiCtMin && theCtau < JpsiCtMax && 
	aCand->pt() > JpsiPtMin && aCand->pt() < JpsiPtMax && 
	fabs(theRapidity) > JpsiRapMin && fabs(theRapidity) < JpsiRapMax &&
	isMuonInAccept(muon1) && isMuonInAccept(muon2) &&
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
	
	double SF=-9; //data to MC scale factor
	double SFp=-9; //SF uncertainty 
	double SFm=-9; 
	
	if(bin!=1){
		SF=1.019;

	}

	if(bin==1){
		SF=1.032; 

	}
	
	//SF=0.99;
	//SFp=0.99;
	//SFm=0.99;
	
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

void JPsiAnalyzerPAT::makeCuts() {
	bool print=0; 
	if(print) std::cout << "Begin JPsiAnalyzerPAT::makeCuts() " << endl; 
	if(print) cout << "collAll.isValid: " << collAll.isValid() << endl; 
	int size=0; 
	if (collAll.isValid()) {
		
		for(vector<pat::CompositeCandidate>::const_iterator it=collAll->begin();
			it!=collAll->end();++it) {
			size++; 
			const pat::CompositeCandidate* cand = &(*it);
			// cout << "Now checking candidate of type " << theJpsiCat << " with pt = " << cand->pt() << endl;
			const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
			const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));
			
			double mE = cand->userFloat("MassErr");
			double M = cand->mass();
			
			//	cout << "mass error: " << mE << endl;
			//cout << "Mass: " << M << endl;
			
			TVector3 Mu1(muon1->px(),muon1->py(), muon1->pz()); 
			TVector3 Mu2(muon2->px(),muon2->py(),muon2->pz()); 

			if(muon1->charge()>0){
				diff_p->Fill(TMath::Abs(Mu1.Perp()-muPosP_Gen->Perp())*1000 ); 	
				diff_n->Fill(TMath::Abs(Mu2.Perp()-muNegP_Gen->Perp())*1000 ); 
			}
			else {
				diff_p->Fill(TMath::Abs(Mu2.Perp()-muPosP_Gen->Perp())*1000 ); 
				diff_n->Fill(TMath::Abs(Mu1.Perp()-muNegP_Gen->Perp())*1000 ); 
			}

			// PAT trigger match, 2 muons to the last filter used in the HLT path (new way)
			this->matchMuonToHlt(muon1, muon2);
			
			bool trigOK = false;
			for (unsigned int iTrig = 0 ; iTrig < HLTBitNames_.size() ; iTrig++) {
				if(print) cout << "HLTBinNames: " <<  HLTBitNames_.at(iTrig) << endl; 
				if(print) cout << "Trig Bit: " << mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] << endl; 
				if (mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == 1 ||
					mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == -1 ||
					mapTriggerNameToIntFired_[HLTBitNames_.at(iTrig)] == 2 ) trigOK = true;
			}
			if(print) cout << "trigOK: " << trigOK << endl;
			
			
			reco::GenParticleRef genMu1 = muon1->genParticleRef();
			reco::GenParticleRef genMu2 = muon2->genParticleRef();
			if(print) cout << "Mu1 Available: " << genMu1.isAvailable() << endl; 
			if(print) cout << "Mu2 Available: "  <<genMu2.isAvailable() << endl;
			
			bool isMuMatched = (genMu1.isAvailable() && genMu2.isAvailable() && 
								genMu1->pdgId()*genMu2->pdgId() == -169 && 
								genMu1->momentum().rho() > 2.5 && genMu2->momentum().rho() > 2.5);
			
			double dpt=muPosP_Gen->Perp()-muNegP_Gen->Perp(); 
			
			double theta1=muPosP_Gen->Theta();
			double theta2=muNegP_Gen->Theta();
			
			double dtheta=theta1-theta2; 
			double dphi=muPosP_Gen->DeltaPhi(*muNegP_Gen); 
			
			
			double DR=muPosP_Gen->DeltaR(*muNegP_Gen); 
			
					/*
			if(_isMC && isMuMatched && trigOK) {
				deltapt_match->Fill(dpt); 
				trigger_pass->Fill(JpsiP_Gen->Perp()); 
				Nnum++;}
			*/
						
			if(_isMC && isMuMatched && !trigOK) {
				deltapt_fail->Fill(dpt); 
				difftheta_diffphi->Fill(dphi,dtheta);
				DR_fail->Fill(DR); 
			}
			
			//if (requireTriggerMatching_ && !trigOK) continue;
			if(print) cout << "Candidate OK. " << endl; 
			//	cout << "cand pt: " << cand->pt() << endl; 
			candpt=cand->pt();
			
			// some counter
			

			
			passedTriggerMatch_++;
			//trigger_pass->Fill(cand->pt()); 
			
			bool opp_charge=false; 
			bool ups_kin=false; 
			bool muon1_good=false; 
			bool muon2_good=false; 
			bool seagull=false; 
			bool muon_qual = false; 
			bool dimuon_qual = false; 
			bool rapidity_check=false; 
			
			if(selMuon(muon1) && selMuon(muon2)) muon_qual=true; 
			
			if (muon1->charge() + muon2->charge() == 0) {	  
				opp_charge=true; 
			}//muon charge
			
			double etaP=TMath::Abs(muPosP_Gen->Eta());
			double etaM=TMath::Abs(muNegP_Gen->Eta());
			
			double PtP=muPosP_Gen->Perp();
			double PtM=muNegP_Gen->Perp();
			
			//compare dimuon mass with gen mass
			//fix kinematic cuts
			if(JpsiP_Gen->Perp()>10 && JpsiP_Gen->Perp()<100 && JpsiP_Gen->M()>8.7 && JpsiP_Gen->M() <11.2)
				ups_kin=true; 
			
			if(JpsiP_Gen->Rapidity()<=0.6)
				rapidity_check=true;
			
			if( ((etaP<1.2 && PtP>4.5) || (etaP>1.2 && etaP<1.4 && PtP>3.5) || (etaP>1.4 && etaP<1.6 && PtP>3.0)) && etaP<1.6)
				muon1_good=true;
			if( ((etaM<1.2 && PtM>4.5) || (etaM>1.2 && etaM<1.4 && PtM>3.5) || (etaM>1.4 && etaM<1.6 && PtM>3.0)) && etaM<1.6)
				muon2_good=true;
			
			if(muPosP_Gen->DeltaPhi(*muNegP_Gen)<0)
				seagull=true; 
			if(selDimuon(cand)) dimuon_qual=true; 
			//I think I need to add a line. selMuon(muon1), selMuon(muon2)
			
			//If an event passes all the dimuon kinematic cuts (we don't want efficiencies for dimuons outside of the kinematic regions of interest) 
			//then compute the single muon efficiency, 
			
			//emu=N(PS,offline,quality,trigger)/N(PS)
			//then compute emumu=N(PS,offline,quality,trigger)/N(PS)
			//but here quality requires BOTH muons to pass the selMuon function.
			
			if(TMath::Abs(etaP)>1.6 ) cout << "muon1 check: " << muon1_good << " " << etaP << endl;  
			if(TMath::Abs(etaM)>1.6 ) cout << "muon2 check: " << muon2_good << " " << etaM << endl; 
			
			//Normally need to include the following statemets: 
				/*		
			if(_isMC && isMuMatched && trigOK && opp_charge && ups_kin && muon1_good && muon2_good && seagull && muon_qual && dimuon_qual && rapidity_check){
				deltapt_match->Fill(dpt);
				double weight=muonEff(muPosP_Gen->Perp(),etaP,"nom")*muonEff(muNegP_Gen->Perp(),etaM,"nom"); //eff1*eff2
				//double weightP=muonEff(muPosP_Gen->Perp(),etaP,"plus")*muonEff(muNegP_Gen->Perp(),etaM,"plus");
				//double weightM=muonEff(muPosP_Gen->Perp(),etaP,"minus")*muonEff(muNegP_Gen->Perp(),etaM,"minus");
				trigger_pass->Fill(JpsiP_Gen->Perp(),1./weight); 	
				trigger_passP->Fill(JpsiP_Gen->Perp(),1./weight);
				trigger_passM->Fill(JpsiP_Gen->Perp(),1./weight);
				DR_pass->Fill(DR); 
				_thePassedCands.push_back(cand); 
			}
				 */
			
			if(_isMC && isMuMatched && seagull) {
				kin_pass->Fill(JpsiP_Gen->Perp()); 
				Ndenom++;
			}
			
			if(_isMC && isMuMatched && seagull && trigOK) 
				SG_num->Fill(JpsiP_Gen->Perp()); 
			if(_isMC && isMuMatched && trigOK)
				SG_den->Fill(JpsiP_Gen->Perp()); 
			
			//&& muon_qual
			if(_isMC && isMuMatched && trigOK && muon_qual && seagull){
				double weight=muonEff(muPosP_Gen->Perp(),etaP,"nom")*muonEff(muNegP_Gen->Perp(),etaM,"nom"); //eff1*eff2
				trigger_pass->Fill(JpsiP_Gen->Perp(),1./weight); 	
				if(JpsiP_Gen->Perp()>10 && JpsiP_Gen->Perp()<50){
					Mass->Fill(M,1./weight);
					if(M>9.26 && M<=9.66)dm->Fill(mE*1000,1./weight); 
					dm_m->Fill(M,mE*1000,1./weight);
				}
				
			}
			
						
		} // it thingy
		
	
	} // for loop 
    if(print) cout << "size of collection: " << size << endl; 
  if(print) std::cout << "End JPsiAnalyzerPAT::makeCuts() " << endl;
	
  return;
}

int JPsiAnalyzerPAT::theBestQQ() {
	//gets candidate with the maximum vertex Probability 
	//std::cout << "theBestQQ. " << endl; 
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
	//std::cout << "isMuoninAccept." << endl; 
   // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )

	bool r1=false; 
	bool r2=false; 
	bool r3=false; 
	
	double eta=TMath::Abs(aMuon->eta());
	double pt=aMuon->pt();
	
	
	r1=eta<1.6 && eta>1.4 && pt>3.0; 
	r2=eta<1.4 && eta>1.2 && pt>3.5; 
	r3=eta<1.2 && pt>4.5;
	
	bool RR=false;
	
	if (r1 || r2 || r3) {
		RR=true;
	}
	
	return RR;
	
   // *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
   // by just returning TRUE
   //  return true;
}

bool
JPsiAnalyzerPAT::selMuon(const pat::Muon* aMuon) {
	
	//std::cout << "selMuon." << endl; 
  TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
  const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

  TrackRef gTrack = aMuon->globalTrack();

  bool trackOK = false;
  int NTracks=static_cast<int>(iTrack->found()); 

  // cooler way of cutting on tracks
	float fHits = static_cast<double>(iTrack->found()) / static_cast<double> (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
  if (_applyExpHitcuts) {
    trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
  // old way of cutting on tracks  
  } else if(NTracks>10) trackOK=true; 
	  
	//cout << "fHits: " << fHits << endl; 
	fHits_hist->Fill(fHits);
	double chi2_dof=static_cast<double> (iTrack->chi2())/static_cast<double> (iTrack->ndof());
	bool TCHI;
	if(chi2_dof<1.8) TCHI=true;
	bool TMA=aMuon->muonID("TrackerMuonArbitrated");
	bool TMT=aMuon->muonID("TMOneStationTight");
	bool TPM=p.pixelLayersWithMeasurement() > 1;
	bool TDXY= fabs(iTrack->dxy(RefVtx)) < 3.0; 
	bool TDZ= fabs(iTrack->dz(RefVtx)) < 15.0 ;
	
	if(TCHI && TMA && TMT && TPM && TDXY && TDZ) iTrack_hist->Fill(static_cast<double>(NTracks)); 
	if(trackOK && TMA && TMT && TPM && TDXY && TDZ) chi2_hist->Fill(chi2_dof); 
	if(trackOK && TCHI && TMT && TPM && TDXY && TDZ) tma_hist->Fill(TMA);
	if(trackOK && TCHI && TMA && TPM && TDXY && TDZ) tmt_hist->Fill(TMT); 
	if(trackOK && TCHI && TMT && TMA  && TDXY && TDZ) tpm_hist->Fill(p.pixelLayersWithMeasurement()); 
		
	return (trackOK && TCHI && aMuon->muonID("TrackerMuonArbitrated") && aMuon->muonID("TMOneStationTight") && (p.pixelLayersWithMeasurement() > 1) && (fabs(iTrack->dxy(RefVtx)) < 3.0) && (fabs(iTrack->dz(RefVtx)) < 15.0 ));
	
}

bool 
JPsiAnalyzerPAT::selDimuon(const pat::CompositeCandidate* aCand) {
	//std::cout<< "selDimuon." << endl; 
	/*
	if (!_applyDiMuoncuts){ 
		if(aCand->pt()>10 && aCand->pt()<100 && aCand->mass()>8.7 && aCand->mass()<11.2) return true;
	}
	 */
	bool check_dca=false; 
	bool check_vProb=false; 
	bool good_dimuon=false; 
	double dca = aCand->userFloat("DCA");
	double vProb=aCand->userFloat("vProb");
	if(dca<0.5) check_dca=true; 
	//normally set vProb>0.01
	if(vProb>0.01) check_vProb=true;
	
	if(_applyDiMuoncuts && check_dca && check_vProb) good_dimuon = true; 
	if(!_applyDiMuoncuts) good_dimuon=true; 
	return good_dimuon;
}

int 
JPsiAnalyzerPAT::getJpsiVarType(const double jpsivar, vector<double> vectbin) {
	//std::cout << "getJpsiVarType." << endl;
  for(unsigned int i=0;i<vectbin.size()-1;i++) {
    if(jpsivar > vectbin[i] && jpsivar < vectbin[i+1]) return i+1;
  }

  return -999;
}

// reset the global DataSet variables
void
JPsiAnalyzerPAT::resetDSVariables(){
	//std::cout << "resetDSVariables. " << endl; 
    //reset J/psi RECO variables
    JpsiMuScleMassCorr=-9999.;
    JpsiMuScleMassErr=-9999.;
    JpsiMassErr=-9999.;
    sigmaPtPos=-9999.;
    sigmaPtNeg=-9999.;
    // JpsiPt=-9999.;
    // JpsiRap=-9999.;
//     JpsiCharge=-9999;
    // JpsiPx=-9999.;
    // JpsiPy=-9999.;
    // JpsiPz=-9999.;
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

//     JpsiType=-1;

    //reset MUON RECO variables
    /* muPosPx=-9999.;
    muPosPy=-9999.;
    muPosPz=-9999.;
    muNegPx=-9999.;
    muNegPy=-9999.;
    muNegPz=-9999.;*/
    if(_isMC){
        MCType=-1;

        //reset J/psi GEN variables
        /* JpsiMass_Gen=-9999.;
        JpsiPt_Gen=-9999.;
        JpsiRap_Gen=-9999.;
        JpsiPx_Gen=-9999.;
        JpsiPy_Gen=-9999.;
        JpsiPz_Gen=-9999.; */

        //reset MUON GEN variables
        /* muPosPx_Gen=-9999.;
        muPosPy_Gen=-9999.;
        muPosPz_Gen=-9999.;
        muNegPx_Gen=-9999.;
        muNegPy_Gen=-9999.;
        muNegPz_Gen=-9999.; */
    }

    //reset EVENT information
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

//! fill Generator Information
void
JPsiAnalyzerPAT::analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
	bool print=0; 
	if(print) std::cout << "Begin JPsiAnalyzerPAT::analyzeGenerator" << endl; 
    using namespace trigger;
	
    std::vector < const reco::Candidate* > genMuons;
    //bool genjpsi= false;
    reco::Candidate::size_type nrD;
	
    //int count= 0;
	
	if(print) cout << "Gen Particle Size: " << genParticles->size() << endl; 
	
	bool missing_muon=false; 
	bool ups_good=false; 
	bool muon_good_p=false; 
	bool isMuMatched=false; 
	if(genParticles->size()!=3) cout << "Gen Particle Size: " << genParticles->size() << endl; 
    for( size_t i = 0; i < genParticles->size(); ++ i )
    {
        // std::cout << "analyzeGenerator: " << i << std::endl;
        const reco::Candidate & cand = (*genParticles)[ i ];
        int Mc_particleID = cand.pdgId();
        if (abs(Mc_particleID) == _oniaPDG && cand.status()==2 )//&& cand.pt() >= 1)
        {
			ups_good=true; 
			if(print)std::cout << "------::analyzeGenerator:: gen " << Mc_particleID <<"'s: pt=" << cand.pt() << "; eta=" << cand.eta() << "; phi=" << cand.phi() << std::endl;
			
			Double_t enGen = sqrt(cand.px()*cand.px() + cand.py()*cand.py() + cand.pz()*cand.pz() + cand.mass()*cand.mass());
			JpsiP_Gen->SetPxPyPzE(cand.px(),cand.py(),cand.pz(),enGen);
			
            //Jpsict_Gen=10.*cand.userFloat("ppdlTrue"));
			
            nrD= cand.numberOfDaughters();
			
			if(nrD!=2) cout << "Number of daughters: " << nrD << endl; 
			
			if(print) cout << "Number of Daughters: " << nrD << endl; 
			
            int count_muon=0;
            for(reco::Candidate::size_type t=0; t < nrD; t++){
                const reco::Candidate* muon= cand.daughter(t);
                int pID = muon->pdgId();
				if(print) std::cout << "------::analyzeGenerator:: gen " << Mc_particleID <<"'s daughter pdgId: " << pID << std::endl;
				
				if(abs(pID) !=13 && cand.daughter(t)->status()!=1) missing_muon=true; 
				
                if (abs(pID) == 13 && cand.daughter(t)->status()==1)
                {
                    genMuons.push_back(muon);
				    if(print) std::cout << "------::analyzeGenerator:: gen " << Mc_particleID <<"'s daughter #: " << count_muon << std::endl;
                    count_muon++;
                }
            }
			
			if( genMuons.empty() && print) cout << "No Muons. break loop" << endl; 
			if(missing_muon) cout << "Missing Muon. PDG ID !=13, status false" << endl; 
            if ( genMuons.empty() || missing_muon ) break;
			
			
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
				
                // fill Polarization variables - gen muons
                Double_t muMass = 0.105658;
				
                Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
                // TLorentzVector *muPos = new TLorentzVector();
                muPosP_Gen->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);
				
                Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
                // TLorentzVector *muNeg = new TLorentzVector();
                muNegP_Gen->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);
				
			if(print) cout << "Gen Upsilon Pt: " << JpsiP_Gen->Perp() << endl; 
			if(print) cout << "Gen Muon+ Pt: " << muPosP_Gen->Perp() << endl;	
			if(print) cout << "Gen Muon- Pt: " << muNegP_Gen->Perp() << endl; 
				
            }
			if(muPosP_Gen->Rho()>2.5 && muNegP_Gen->Rho()>2.5){
				muon_good_p=true; 
			}
			else {
				if(print) cout << "Muon Momentum < 2.5GeV/C" << endl; 
			}
			
			if(muon_good_p && !missing_muon && ups_good) isMuMatched=true;  
			
        } // end loop over genParticles
    }
	
	 bool ups_kin=false; 
	 bool muon1_good=false; 
	 bool muon2_good=false; 
	 bool seagull=false; 
	 bool rapidity_check=false;
	 
	double etaP=TMath::Abs(muPosP_Gen->Eta());
	double etaM=TMath::Abs(muNegP_Gen->Eta());
	
	double PtP=muPosP_Gen->Perp();
	double PtM=muNegP_Gen->Perp();
	
//compare dimuon mass with gen mass
	//fix kinematic cuts
	if(JpsiP_Gen->Perp()>10 && JpsiP_Gen->Perp()<100 && JpsiP_Gen->M()>8.7 && JpsiP_Gen->M() <11.2)
	 ups_kin=true; 
    
	if( ((etaP<1.2 && PtP>4.5) || (etaP>1.2 && etaP<1.4 && PtP>3.5) || (etaP>1.4 && etaP<1.6 && PtP>3.0)) && etaP<1.6)
		muon1_good=true;
	if( ((etaM<1.2 && PtM>4.5) || (etaM>1.2 && etaM<1.4 && PtM>3.5) || (etaM>1.4 && etaM<1.6 && PtM>3.0)) && etaM<1.6)
		muon2_good=true;
	
	if(JpsiP_Gen->Rapidity()<=0.6)
		rapidity_check=true;
	
	if(muPosP_Gen->DeltaPhi(*muNegP_Gen)<0)
	 seagull=true; 
	 //include seaguall line
	 if(ups_kin && muon1_good && muon2_good  && isMuMatched && rapidity_check) pass_cuts=true; 
	 
	if(pass_cuts)gen_mass->Fill(JpsiP_Gen->M()); 
		
	if(print) std::cout << "End JPsiAnalyzerPAT::analyzeGenerator" << endl; 
}

void
JPsiAnalyzerPAT::hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup)
{
	bool print=0; 
    if(print)	std::cout << "Begin JPsiAnalyzerPAT::hltReport" << endl; 
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
		if(print) cout << "size of hltConfig_: " << n << endl; 
		
		for(unsigned int i=0; i<n; i++){
			//cout << hltConfig_.triggerName(i) << endl;
		}
		for (std::map<std::string, unsigned int>::iterator it = mapTriggernameToHLTbit.begin(); it != mapTriggernameToHLTbit.end(); it++) {
			unsigned int triggerIndex= hltConfig_.triggerIndex( it->first );
			if(print) cout << "Trigger Index: " << triggerIndex << endl; 
			if (triggerIndex >= n) {
				if(print) std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " not available in config!" << std::endl;
			}
			else {
				it->second= triggerIndex;
				if(print) std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " available in config!" << std::endl;
			}
		}
    }//if hltConfigInit_
    
    // Get Trigger Results
    try {
		iEvent.getByLabel( tagTriggerResults_, handleTriggerResults_ );
		if(print) cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResult is present in current event" << endl;
    }
    catch(...) {
		if(print) cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults NOT present in current event" << endl;
    }
	
	if(handleTriggerResults_.isValid()!=1) cout << "handleTriggerResults_.isValid: " << handleTriggerResults_.isValid() << endl; 
	if(hltConfigInit_!=1) cout  << "hltInit: " << hltConfigInit_ << endl; 

	if(print) std::cout << "Before check: " << std::endl; 
	
	if(hltConfigInit_){
		for(unsigned int i=0; i<hltConfig_.size(); i++){
			TString name(hltConfig_.triggerName(i)); 
			cout << "Trigger Config: " << name << endl; 
			bool interesting_trigg=false; 
			if (name.Contains("Upsilon") && name.Contains("Barrel")) interesting_trigg=true; 
			if (name.Contains("HLT_Mu40_v")) interesting_trigg=true; 
			if(interesting_trigg && print) cout << name << " Bit: " << handleTriggerResults_->accept( mapTriggernameToHLTbit[hltConfig_.triggerName(i)] ) << " map: " << mapTriggernameToHLTbit[hltConfig_.triggerName(i)] << endl; 
		}
	}	
	
    if ( handleTriggerResults_.isValid() ){
		if(print) cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults IS valid in current event" << endl;
		
		// loop over Trigger Results to check if paths was fired
		for(std::vector< std::string >::iterator itHLTNames= HLTBitNames_.begin(); itHLTNames != HLTBitNames_.end(); itHLTNames++){
			const std::string triggerPathName =  *itHLTNames;
			if(print) std::cout << "[FloJPsiAnalyzer::hltReport] --- TriggerName --- TriggerName LOOP" << std::endl;
			if(print) cout << "Trigger Names: " << triggerPathName << endl; 
			if ( mapTriggernameToHLTbit[triggerPathName] < 1000 ) {
				//the number mapTriggernameToHLTbit takes the trigger path to the trigger index. 
				if (handleTriggerResults_->accept( mapTriggernameToHLTbit[triggerPathName] ) ){
					if(print) cout << "triggerPathFired: " << triggerPathName << endl; 
					mapTriggerNameToIntFired_[triggerPathName] = 3;
				}//if accept

	 //-------prescale factor------------
				if (!_isMC) {
					const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerPathName));
					if(print) std::cout << "[FloJPsiAnalyzer::prescalvalues] --- TriggerName"<<triggerPathName<<" prescales first "<< prescales.first <<" prescales second "<< prescales.second <<std::endl;
					mapTriggerNameToPrescaleFac_[triggerPathName] = prescales.first * prescales.second;
				}//!isMC_
			}//if there is a trigger name
		}//loop 
		bool keep_trigger=false; 
	    if(print) cout << "After Check: " << endl;
		if(hltConfigInit_){
			for(unsigned int i=0; i<hltConfig_.size(); i++){
				TString name(hltConfig_.triggerName(i)); 
				bool interesting_trigg=false; 
				if (name.Contains("Upsilon") && name.Contains("Barrel")) interesting_trigg=true; 
				//if (name.Contains("HLT_Mu40_v")) interesting_trigg=true; 
				if(interesting_trigg){
					if(print) cout << name << " Bit: " << mapTriggerNameToIntFired_[hltConfig_.triggerName(i)] << endl; 
					if(mapTriggerNameToIntFired_[hltConfig_.triggerName(i)]==3){
						keep_trigger=true; 
					}
				}
			}
		}//second trigger check
		HLT_tot++; 
		if(keep_trigger ) {
			if(print) cout << "good trigger. " << endl; 
			HLT_trig++; 
		}
    } else cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerResults NOT valid in current event" << endl;
	
	if(print)	std::cout << "End JPsiAnalyzerPAT::hltReport" << endl << endl; 
}

void
JPsiAnalyzerPAT::matchMuonToHlt(const pat::Muon* muon1, const pat::Muon* muon2)
{
	bool print=0; 
	if(print) std::cout << "Begin JPsiAnalyzerPAT::matchMuonToHlt " << endl; 
    std::string HLTL3MuCollName = "hltL3MuonCandidates::" + tagTriggerResults_.process();
    std::string HLTL2MuCollName = "hltL2MuonCandidates::" + tagTriggerResults_.process();
    std::string HLTTrackCollName = "hltMuTrackJpsiCtfTrackCands::" + tagTriggerResults_.process();
    std::string HLTTkMuCollName = "hltMuTkMuJpsiTrackerMuonCands::" + tagTriggerResults_.process();
    
    //! Loop over Trigger Paths and match muons to last Filter/collection
    for ( std::map<std::string, int>::iterator it = mapTriggerNameToIntFired_.begin(); it != mapTriggerNameToIntFired_.end(); it ++ ) {

        std::string triggerName = it->first;
        std::string hltLastFilterName = mapTriggerToLastFilter_[triggerName];
	    if(print) cout << "Trigger Name: " << triggerName << " Last Filter: " << hltLastFilterName << endl; 
		if(print) cout << "Fired " << mapTriggerNameToIntFired_[triggerName] << " second: " << it->second << endl; 
		
        //! just use Triggers which are in TriggerResults; value == 3
        if ( it->second != 3 ){
			if(print) cout << "no trigger info. " << endl; 
			continue;
		}
			
        const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter( hltLastFilterName );
        const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter( hltLastFilterName );
        bool pass1 = mu1HLTMatches.size() > 0;
        bool pass2 = mu2HLTMatches.size() > 0; 
		
		if(print) cout << "Passes Last Filter: " << endl; 
		if(print) cout << "Muon1 HLT Match: " << pass1 << " Muon2 HLT Match: " << pass2 << endl; 

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
	/*	
	cout << "The Check: " << *theCheck << endl; 
	cout << "end: " << *HLTBitNames_DoubleMu.end() << endl; 
	cout << "Pass1: " << pass1 << " Pass2: " << pass2 << endl;	
	*/ 
  if(print){
	if(theCheck == HLTBitNames_DoubleMu.end()) cout << "Failed theCheck " << endl; 	
	if(pass1==false) cout << "Failed on Muon1 Check." << endl;
	if(pass2==false) cout << "Failed on Muon2 Check. " << endl; 	
	if(pass1==false && pass2==false) cout << "Failed both Muons." << endl; 	
	if ( it->second != 3 && pass1==true && pass2==true) cout << "trigger information does not exist, but pass1 and pass2 are true" << endl;
  }
	if (theCheck != HLTBitNames_DoubleMu.end() && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
	theCheck = std::find(HLTBitNames_SingleMu.begin(),HLTBitNames_SingleMu.end(),triggerName);
	if (theCheck != HLTBitNames_SingleMu.end() && (pass1 == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
    }
	
	if(print) std::cout << "End JPsiAnalyzerPAT::matchMuonToHlt " << endl << endl; 
	
}

void 
JPsiAnalyzerPAT::muonStationDistance (const pat::CompositeCandidate* aCand) 
{
	//std::cout << "muonStationDistance." << endl; 
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
