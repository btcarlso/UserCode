
// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.cc
//
/*

 Description: an example analyzer using SUSY Ntuples

 Implementation:
 The macro is driven by ana.C in this directory.
 Fills MET and diEMPt histograms for gg and ff events with >= 1 good jets. gg and ff
 skimming can be done simultaneously by setting CopyEvents to true in ana.C.

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.cc,v 1.9 2013/06/26 20:26:20 weinberg Exp $
//

#include <TH1F.h>
#include <TFile.h>

#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

// REMOVE THE LINES BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h" // to access the JEC scales
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" // to access the uncertainties

#include "SusyEventAnalyzer.h"



template<typename T>
bool
PtGreater(const T* p1, const T* p2) {
  return p1->momentum.Pt() > p2->momentum.Pt();
}



template<typename T1, typename T2>
bool
isSameObject(const T1& p1, const T2& p2, double dRCut = 0.1)
{
  float dEta = p1.momentum.Eta() - p2.momentum.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.momentum.Phi() - p2.momentum.Phi());
  return dEta*dEta + dPhi*dPhi < dRCut * dRCut;
}



template<typename T1, typename T2>
bool
isSameObject2(const T1& p1, const T2& p2, double dRCut = 0.1)
{
  float dEta =                      p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());

  return dEta * dEta + dPhi * dPhi < dRCut * dRCut;
}



template<typename T1, typename T2>
Float_t
getDeltaR(T1& p1, T2& p2)
{
  Float_t dEta =                      p1.Eta() - p2.Eta();
  Float_t dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}



Float_t
muonEffectiveAreas(Float_t _eta)
{
  // Source: SWGuideMuonId twiki

  if      (_eta < 1.0) return 0.674;
  else if (_eta < 1.5) return 0.565;
  else if (_eta < 2.0) return 0.442;
  else if (_eta < 2.2) return 0.515;
  else if (_eta < 2.3) return 0.821;
  else if (_eta < 2.4) return 0.660;
                       return 0.000;
}



Float_t
electronEffectiveAreas(Float_t _eta)
{
  // Source: EgammaEARhoCorrection twiki
  if      (_eta < 1.0)   return 0.13;
  else if (_eta < 1.479) return 0.137;
  else if (_eta < 2.0)   return 0.067;
  else if (_eta < 2.2)   return 0.089;
  else if (_eta < 2.3)   return 0.107;
  else if (_eta < 2.4)   return 0.11;
                         return 0.138;
}



void
photonEffectiveAreas(Float_t _eta, Float_t* _effA)
{
  Float_t& effACH(_effA[0]);
  Float_t& effANH(_effA[1]);
  Float_t& effAPh(_effA[2]);

  // Source: CutBasedPhotonID2012 twiki
  if (_eta < 1.) {
    effACH = 0.012;
    effANH = 0.03;
    effAPh = 0.148;
  }
  else if (_eta < 1.479) {
    effACH = 0.010;
    effANH = 0.057;
    effAPh = 0.13;
  }
  else if (_eta < 2.) {
    effACH = 0.014;
    effANH = 0.039;
    effAPh = 0.112;
  }
  else if (_eta < 2.2) {
    effACH = 0.012;
    effANH = 0.015;
    effAPh = 0.216;
  }
  else if (_eta < 2.3) {
    effACH = 0.016;
    effANH = 0.024;
    effAPh = 0.262;
  }
  else if (_eta < 2.4) {
    effACH = 0.02;
    effANH = 0.039;
    effAPh = 0.26;
  }
  else {
    effACH = 0.012;
    effANH = 0.072;
    effAPh = 0.266;
  }
}



Float_t
SusyEventAnalyzer::SetSusyXSec(Float_t mSquark)
{
  if      (mSquark ==  200.0) return 183390.0;
  else if (mSquark ==  250.0) return  55865.5;
  else if (mSquark ==  300.0) return  19828.3;
  else if (mSquark ==  350.0) return   7982.94;
  else if (mSquark ==  400.0) return   3543.38;
  else if (mSquark ==  450.0) return   1679.47;
  else if (mSquark ==  500.0) return    847.051;
  else if (mSquark ==  550.0) return    447.055;
  else if (mSquark ==  600.0) return    244.862;
  else if (mSquark ==  650.0) return    137.949;
  else if (mSquark ==  700.0) return     79.9667;
  else if (mSquark ==  750.0) return     47.3374;
  else if (mSquark ==  800.0) return     28.4146;
  else if (mSquark ==  850.0) return     17.2924;
  else if (mSquark ==  900.0) return     10.6744;
  else if (mSquark ==  950.0) return      6.68877;
  else if (mSquark == 1000.0) return      4.24173;
  else if (mSquark == 1050.0) return      2.71284;
  else if (mSquark == 1100.0) return      1.7516;
  else if (mSquark == 1150.0) return      1.13551;
  else if (mSquark == 1200.0) return      0.739516;
  else if (mSquark == 1250.0) return      0.484981;
  else if (mSquark == 1300.0) return      0.319675;
  else if (mSquark == 1350.0) return      0.212278;
  else if (mSquark == 1400.0) return      0.14128;
  else                        return     -1.0;
}



Float_t
SusyEventAnalyzer::SetSignalEventWeight(TString dsName, Float_t lumi, Float_t m)
{
  if (dsName.Contains("stealthT2")) {
    if (dsName.Contains("WW"))
      return lumi * SetSusyXSec(m) * 0.543 / 20000.0;  // Efficiency of WW -> at least 1 lepton (including tau) ~ 54.3% from PDG

    return lumi * SetSusyXSec(m) / 20000.0;
  }

  return 1.0;
}



void
SusyEventAnalyzer::createHistogram(const char* name,   const char* title,
                                   const char* xTitle, const char* yTitle,
                                   Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
  TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);

  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);

  h->Sumw2();

  hName[name] = h;
}



////////// MAIN ANALYSIS FUNCTION //////////
void
SusyEventAnalyzer::Run()
{
  float SumW_Top=0;
  float SumW_TopM=0;
  float SumW_TopP=0;

  int nCnt[20];
  for(int i(0); i<20; i++) nCnt[i] = 0;

  ////////// TEXT OUTPUT //////////

  std::ofstream outFile;
  if(logFileName != "cout"){
    outFile.open(logFileName);
    if(!outFile.is_open()){
      std::cerr << "Log output " << logFileName << " could not be opened." << std::endl;
      return;
    }
  }
  std::ostream& out(logFileName == "cout" ? std::cout : outFile);

  ////////// SWITCH FOR USING USER-PROVIDED JEC //////////
  bool useCustomJEC(false);

  // REMOVE THE LINES BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
  FactorizedJetCorrector* jetCorrection(0);
  JetCorrectionUncertainty* jecUncertainty(0);

  TFile* ggFile(0);
  TFile* ffFile(0);
  TFile* fout(0);

  try{

    ////////// INITIALIZE SKIMMED CLONE //////////
    TTree* ggTree(0);
    TTree* ffTree(0);
    if(copyEvents){

      if(printLevel > 0) out << "Open file for skim output" << std::endl;

      ggFile = TFile::Open("susyEvents_" + outputName + "_gg.root", "RECREATE");
      ffFile = TFile::Open("susyEvents_" + outputName + "_ff.root", "RECREATE");
      if(!ggFile || ggFile->IsZombie() || !ffFile || ffFile->IsZombie()){
        std::cerr << "Cannot open output file susyEvents_" << outputName << ".root" << std::endl;
        throw std::runtime_error("IOError");
      }
      else{
        ggFile->cd();
        ggTree = new TTree("susyTree", "SUSY Event");
        ggTree->SetAutoSave(10000000);

        ffFile->cd();
        ffTree = new TTree("susyTree", "SUSY Event");
        ffTree->SetAutoSave(10000000);

        /* Register the output trees to the Event object - the branches will be booked internally */
        event.addOutput(*ggTree);
        event.addOutput(*ffTree);
      }

    }

    // Ntuple vectors
    muon_e          = new std::vector<float>(); // tight muons
    muon_px         = new std::vector<float>();
    muon_py         = new std::vector<float>();
    muon_pz         = new std::vector<float>();
    muon_charge     = new std::vector<float>();
	  
	loose_muon_e          = new std::vector<float>(); //loose muons
	loose_muon_px         = new std::vector<float>();
	loose_muon_py         = new std::vector<float>();
	loose_muon_pz         = new std::vector<float>();
	loose_muon_charge     = new std::vector<float>();  
	  
    electron_e      = new std::vector<float>();
    electron_px     = new std::vector<float>();
    electron_py     = new std::vector<float>();
    electron_pz     = new std::vector<float>();
    electron_charge = new std::vector<float>();
	  
	loose_electron_e      = new std::vector<float>();
	loose_electron_px     = new std::vector<float>();
	loose_electron_py     = new std::vector<float>();
	loose_electron_pz     = new std::vector<float>();
	loose_electron_charge = new std::vector<float>(); 
	  
    photon_e        = new std::vector<float>();
    photon_px       = new std::vector<float>();
    photon_py       = new std::vector<float>();
    photon_pz       = new std::vector<float>();
    jet_e           = new std::vector<float>();
    jet_px          = new std::vector<float>();
    jet_py          = new std::vector<float>();
    jet_pz          = new std::vector<float>();
	jet_bTagL        = new std::vector<bool>();
    jet_bTagM        = new std::vector<bool>();
	jet_bTagT        = new std::vector<bool>();
	jet_algFlavor        = new std::vector<int>();
	jet_phyFlavor        = new std::vector<int>();
    jet_unc         = new std::vector<float>();

    w_px            = new std::vector<float>(); 
    w_py            = new std::vector<float>();
    w_pz            = new std::vector<float>();
    w_e             = new std::vector<float>();

    z_px            = new std::vector<float>();
    z_py            = new std::vector<float>();
    z_pz            = new std::vector<float>();
    z_e             = new std::vector<float>();

    top_px            = new std::vector<float>();
    top_py            = new std::vector<float>();
    top_pz            = new std::vector<float>();
    top_e             = new std::vector<float>();


    ////////// INITIALIZE HISTOGRAM / NTUPLE OUTPUT //////////
    if(printLevel > 0) out << "Open file for histograms / ntuples" << std::endl;

    fout = TFile::Open(outputDirectoryName + "/" + outputName + ".root", "RECREATE");
    if(!fout || fout->IsZombie()){
      std::cerr << "Cannot open output file hist_" << outputName << ".root" << std::endl;
      throw std::runtime_error("IOError");
    }

    if(printLevel > 0) out << "Define histograms / ntuples" << std::endl;

    fout->cd();

    TTree* tree = new TTree("tree", "Stealth SUSY tree");

    tree->Branch("runNo",   & runNo,   "runNo/I");
    tree->Branch("lumiNo",  & lumiNo,  "lumiNo/I");
    tree->Branch("eventNo", & eventNo, "eventNo/l");

    tree->Branch("squark_m",     & squark_m,     "squark_m/F");
    tree->Branch("neutralino_m", & neutralino_m, "neutralino_m/F");

    tree->Branch("vertices_n",  & vertices_n,  "vertices_n/I");
    tree->Branch("muons_n",     & muons_n,     "muons_n/I");
    tree->Branch("electrons_n", & electrons_n, "electrons_n/I");
    tree->Branch("photons_n",   & photons_n,   "photons_n/I");
    tree->Branch("jets_n",      & jets_n,      "jets_n/I");

    tree->Branch("eventWeight", & evtWt, "eventWeight/F");

	tree->Branch("puWeight_nom",  & puWt_nom,  "puWeight_nom/F");
	tree->Branch("puWeight_up",   & puWt_up,   "puWeight_up/F");
	tree->Branch("puWeight_down", & puWt_down, "puWeight_dowm/F");
	tree->Branch("NumInteractions",  & NumInteractions,  "NumInteractions/I");//number of true interactions for pile-up re-weighting

	  
    tree->Branch("met_et",  & met_et,  "met_et/F");
    tree->Branch("met_phi", & met_phi, "met_phi/F");
    tree->Branch("st",      & st,      "st/F");
	tree->Branch("gen_st",      & gen_st,      "gen_st/F");
	tree->Branch("gen_ht",      & gen_ht,      "gen_ht/F");

//tight muons
    tree->Branch("muon_e",          "vector<float>", & muon_e);
    tree->Branch("muon_px",         "vector<float>", & muon_px);
    tree->Branch("muon_py",         "vector<float>", & muon_py);
    tree->Branch("muon_pz",         "vector<float>", & muon_pz);
    tree->Branch("muon_charge",     "vector<float>", & muon_charge);
	  
	tree->Branch("loose_muon_e",          "vector<float>", & loose_muon_e);
	tree->Branch("loose_muon_px",         "vector<float>", & loose_muon_px);
	tree->Branch("loose_muon_py",         "vector<float>", & loose_muon_py);
	tree->Branch("loose_muon_pz",         "vector<float>", & loose_muon_pz);
	tree->Branch("loose_muon_charge",     "vector<float>", & loose_muon_charge); 
	  
    tree->Branch("electron_e",      "vector<float>", & electron_e);
    tree->Branch("electron_px",     "vector<float>", & electron_px);
    tree->Branch("electron_py",     "vector<float>", & electron_py);
    tree->Branch("electron_pz",     "vector<float>", & electron_pz);
    tree->Branch("electron_charge", "vector<float>", & electron_charge);
	
	tree->Branch("loose_electron_e",      "vector<float>", & loose_electron_e);
	tree->Branch("loose_electron_px",     "vector<float>", & loose_electron_px);
	tree->Branch("loose_electron_py",     "vector<float>", & loose_electron_py);
	tree->Branch("loose_electron_pz",     "vector<float>", & loose_electron_pz);
	tree->Branch("loose_electron_charge", "vector<float>", & loose_electron_charge);
	  
    tree->Branch("photon_e",        "vector<float>", & photon_e);
    tree->Branch("photon_px",       "vector<float>", & photon_px);
    tree->Branch("photon_py",       "vector<float>", & photon_py);
    tree->Branch("photon_pz",       "vector<float>", & photon_pz);
    tree->Branch("jet_e",           "vector<float>", & jet_e);
    tree->Branch("jet_px",          "vector<float>", & jet_px);
    tree->Branch("jet_py",          "vector<float>", & jet_py);
    tree->Branch("jet_pz",          "vector<float>", & jet_pz);
	tree->Branch("jet_bTagL",       "vector<bool>",  & jet_bTagL);//loose b-tags
    tree->Branch("jet_bTagM",       "vector<bool>",  & jet_bTagM);//medium wp b-tags
	tree->Branch("jet_bTagT",       "vector<bool>",  & jet_bTagT);//medium wp b-tags
	tree->Branch("jet_algFlavor",   "vector<int>",  & jet_algFlavor);//alg flavor
	tree->Branch("jet_phyFlavor",   "vector<int>",  & jet_phyFlavor);//phy flavor

    tree->Branch("jet_unc",         "vector<float>", & jet_unc);

    tree->Branch("w_px",            "vector<float>", & w_px); 
    tree->Branch("w_py",          "vector<float>", & w_py);
    tree->Branch("w_pz",          "vector<float>", & w_pz);
    tree->Branch("w_e",          "vector<float>", & w_e);

    tree->Branch("z_px",          "vector<float>", & z_px);
    tree->Branch("z_py",          "vector<float>", & z_py);
    tree->Branch("z_pz",          "vector<float>", & z_pz);
    tree->Branch("z_e",          "vector<float>", & z_e);

    tree->Branch("top_px",          "vector<float>", & top_px);
    tree->Branch("top_py",          "vector<float>", & top_py);
    tree->Branch("top_pz",          "vector<float>", & top_pz);
    tree->Branch("top_e",          "vector<float>", & top_e);

	//create histogram
	  createHistogram("h_Zmass_nGamma0","Zmass, N_{#gamma}=0","M_{Z} [GeV]", "Events", 100,0,200); 
	  createHistogram("h_Zmass","Zmass","M_{Z} [GeV]", "Events", 100,0,200); 
	  createHistogram("h_Nevents","Nevents","N_{Events}","",10,0,10); 
	  //createHistogram("h_mumu","dimuon mass","M_{#mu#mu} [GeV]", "Events", 100,0,200); 

	  
    ////////// INITIALIZE JEC //////////
    // REMOVE THE LINES BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
    if(useCustomJEC){
      if(printLevel > 0) out << "Initialize jet energy corrections" << std::endl;

      std::string jecSourcePrefix("../jec/FT_53_V21_AN3_");
      jetCorrection = new FactorizedJetCorrector("L1FastJet:L2Relative:L3Absolute:L2Relative",
                                       jecSourcePrefix + "L1FastJet_AK5PF.txt:" +
                                       jecSourcePrefix + "L2Relative_AK5PF.txt:" +
                                       jecSourcePrefix + "L3Absolute_AK5PF.txt:" +
                                       jecSourcePrefix + "L2RelativeL3AbsoluteResidual_AK5PF.txt");
      jecUncertainty = new JetCorrectionUncertainty(jecSourcePrefix + "Uncertainty_AK5PF.txt");
    }

    ////////// SET MET FILTER COMBINATION //////////
    if(printLevel > 0) out << "Set MET filter combination" << std::endl;

    /*
      Recommended set of MET filters depends on the dataset being used. Analysists should check the JetMET twiki
      MissingETOptionalFilters (for general info) and PdmVKnowFeatures (for specifics, especially on HcalLaser
      filter) and modify the configuration below accordingly. The example below follows the JetMET recommendation
      for 22Jan2013 rereco datasets.
    */
    /*
      There are two ways to apply a combined filter. One is to simply take the logical AND of the individual filters:

      event.passMetFilter(susy::kCSCBeamHalo) && event.passMetFilter(susy::kHcalNoise) && ... (for each event)

      and the other is to use the bit mask below:

      event.metFilterMask = (1 << susy::kCSCBeamHalo) | (1 << susy::kHcalNoise) | ...; (once before the event loop)
      event.passMetFilters(); (for each event)

      Both will yield exactly same results.
    */
    event.metFilterMask =
      (1 << susy::kCSCBeamHalo) | (1 << susy::kHcalNoise) | (1 << susy::kEcalDeadCellTP) | (1 << susy::kHcalLaserRECOUserStep) | (1 << susy::kTrackingFailure) |
      (1 << susy::kEEBadSC) | (1 << susy::kLogErrorTooManyClusters) | (1 << susy::kLogErrorTooManyTripletsPairs) | (1 << susy::kLogErrorTooManySeeds);


    /////////////////////////////////////
    ////////// MAIN EVENT LOOP //////////
    /////////////////////////////////////

    if(printLevel > 0) out << "Start event loop" << std::endl;

//     Long64_t nEntries = event.getEntries();
    Long64_t nEntries = fTree->GetEntries();

    long iEntry(0);
    while(iEntry != processNEvents && event.getEntry(iEntry++) != 0){

      if(printLevel > 1 || iEntry % printInterval == 0)
        out << "Begin processing event " << iEntry - 1 << ". Current: run=" << event.runNumber << ", event=" << event.eventNumber << std::endl;

      ////////// CLEAR FLAT NTUPLE VECTORS //////////
      muon_e->clear();     muon_px->clear();     muon_py->clear();     muon_pz->clear();     muon_charge->clear();
	  loose_muon_e->clear();     loose_muon_px->clear();     loose_muon_py->clear();     loose_muon_pz->clear();     loose_muon_charge->clear();
	  electron_e->clear(); electron_px->clear(); electron_py->clear(); electron_pz->clear(); electron_charge->clear();
	  loose_electron_e->clear(); loose_electron_px->clear(); loose_electron_py->clear(); loose_electron_pz->clear(); loose_electron_charge->clear();

	  photon_e->clear();   photon_px->clear();   photon_py->clear();   photon_pz->clear();
		jet_e->clear();      jet_px->clear();      jet_py->clear();      jet_pz->clear();  jet_bTagL->clear();    jet_bTagM->clear(); jet_unc->clear();jet_bTagT->clear(); 
		jet_algFlavor->clear(); jet_phyFlavor->clear();
      w_px->clear();       w_py->clear();        w_pz->clear();        w_e->clear();
      z_px->clear();       z_py->clear();        z_pz->clear();        z_e->clear();    
      top_px->clear();     top_py->clear();      top_pz->clear();      top_e->clear();  

      ////////// DUMP EVENT CONTENT //////////
      if(printLevel > 2) event.Print(out);

      /* total number of events */
      nCnt[0]++;

      ////////// GOOD LUMI FILTER //////////
      if(printLevel > 1) out << "Apply good lumi list." << std::endl;
      if(!IsGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;

      ////////// MET FILTER //////////
      if(event.isRealData){
        if(printLevel > 1) out << "Apply MET filter." << std::endl;
        if(!event.passMetFilters()) continue;
      }

      ////////// TRIGGER //////////
      if(printLevel > 1) out << "Apply HLT cut." << std::endl;
      if(!PassTriggers()) continue;

      ////////// REQUIRE AT LEAST ONE GOOD VERTEX //////////
      if(printLevel > 1) out << "Require at least one good vertex." << std::endl;
      unsigned nV(event.vertices.size());
      unsigned iV(0);
      for(; iV != nV; ++iV){
        susy::Vertex& vertex(event.vertices[iV]);
        if(vertex.ndof >= 4 && std::abs(vertex.position.Z()) < 24. && vertex.position.Perp() < 2.) break;
      }
      if(iV == nV) continue;

      susy::Vertex const& primVtx(event.vertices[iV]);

      if(printLevel > 1) out << "Event passes preselection." << std::endl;

      /* number of events passing preselection */
      nCnt[1]++;

      ////////// SELECT MUONS //////////
      if (printLevel > 1) out << "Find muons in the event" << std::endl;

      std::vector<susy::Muon const*> looseMuons;
      std::vector<susy::Muon const*> tightMuons;

      susy::MuonCollection const& muons(event.muons["muons"]);

      for (size_t iMuon(0); iMuon != muons.size(); ++iMuon) {
        susy::Muon const& muon(muons[iMuon]);

        if (!muon.isPFMuon())                              continue;
        if (!muon.isTrackerMuon() && !muon.isGlobalMuon()) continue;
        if (muon.momentum.Pt() < muon2_ptCut)              continue;

        Float_t effA          = muonEffectiveAreas(fabs(muon.momentum.Eta()));
        Float_t pfCombinedIso = (muon.sumChargedHadronPt04 + std::max(double (muon.sumNeutralHadronEt04 +
                                 muon.sumPhotonEt04 - event.rho25 * effA), 0.0)) / muon.momentum.Pt();

        bool normChi2Cut        = false;
        bool dxyCut             = false;
        bool dzCut              = false;
        bool nValidPixelHitsCut = false;

        if (muon.globalTrack != 0)
          normChi2Cut = muon.globalTrack->normChi2() < 10.0;

        if (muon.bestTrack != 0) {
          dxyCut = fabs(muon.bestTrack->dxyPv(primVtx.position)) < 0.2;
          dzCut  = fabs(muon.bestTrack->dzPv (primVtx.position)) < 0.5;
        }

        if (muon.innerTrack != 0)
          nValidPixelHitsCut = muon.innerTrack->numberOfValidPixelHits >= 1;

        bool nValidMuonHitsCut     = muon.nValidMuonHits              >= 1;
        bool nMatchedStationsCut   = muon.nMatchedStations            >= 2;
        bool nTrackerLayersCut     = muon.nPixelLayersWithMeasurement +
                                     muon.nStripLayersWithMeasurement >= 6;
        bool loosePfCombinedIsoCut = pfCombinedIso                    <  0.2;
        bool tightPfCombinedIsoCut = pfCombinedIso                    <  0.12;

        bool looseCut = loosePfCombinedIsoCut;
        bool tightCut = muon.isGlobalMuon() && normChi2Cut       && nValidMuonHitsCut &&
                        nMatchedStationsCut && dxyCut            && dzCut             &&
                        nValidPixelHitsCut  && nTrackerLayersCut && tightPfCombinedIsoCut;

        if (looseCut) looseMuons.push_back(& muon);
        if (tightCut) tightMuons.push_back(& muon);
      }

      ////////// SORT AND CLEAN SELECTED MUONS //////////
      std::sort(looseMuons.begin(), looseMuons.end(), PtGreater<susy::Muon>);
      std::sort(tightMuons.begin(), tightMuons.end(), PtGreater<susy::Muon>);

      // Leading muon must pass higher pT cut
      if (looseMuons.size() >= 1 && looseMuons.at(0)->momentum.Pt() < muon1_ptCut) looseMuons.clear();
      if (tightMuons.size() >= 1 && tightMuons.at(0)->momentum.Pt() < muon1_ptCut) tightMuons.clear();

      // Extra cleaning for loose muons: Same sign muons should have DeltaR > 0.02 to avoid split tracks
      if (looseMuons.size() >= 2) {
        for (size_t i = 0; i < looseMuons.size() - 1; i++) {
          for (size_t j = i + 1; j < looseMuons.size(); j++) {
            if (looseMuons.at(i)->bestTrack->charge == looseMuons.at(j)->bestTrack->charge &&
                getDeltaR(looseMuons.at(i)->momentum, looseMuons.at(j)->momentum) < 0.02)
              looseMuons.erase(looseMuons.begin() + j);
          }
        }
      }

      ////////// SELECT ELECTRONS //////////
      if (printLevel > 1) out << "Find electrons in the event" << std::endl;

      std::vector<susy::Electron const*> looseElectrons;
      std::vector<susy::Electron const*> tightElectrons;

      susy::ElectronCollection const& electrons(event.electrons["gsfElectrons"]);

      for(size_t iEle(0); iEle != electrons.size(); ++iEle) {
        susy::Electron const& electron(electrons[iEle]);

        bool sameAsTightMuon = false;

        for (size_t i = 0; i < tightMuons.size(); i++) {
          if (isSameObject2(electron.momentum, tightMuons.at(i)->momentum, 0.3)) {
            sameAsTightMuon = true;
            break;
          }
        }

        Float_t dInvEInvP = 1.0 / electron.ecalEnergy * (1.0 - electron.eSuperClusterOverP);

        if (sameAsTightMuon)                           continue;
        if (electron.momentum.Et()  < electron2_etCut) continue;
        if (!electron.isEB())                          continue;
        if (electron.sigmaIetaIeta  >  0.01)           continue;
        if (electron.hcalOverEcalBc >  0.12)           continue;
        if (electron.gsfTrack       == 0)              continue;
        if (fabs(dInvEInvP)         >  0.05)           continue;
        if (!electron.passConversionVeto)              continue;

        Float_t effA          = electronEffectiveAreas(fabs(electron.momentum.Eta()));
        Float_t pfCombinedIso = (electron.chargedHadronIso + std::max(double (electron.neutralHadronIso +
                                 electron.photonIso - event.rho25 * effA), 0.0)) / electron.momentum.Et();

        bool looseDEtaInCut       = fabs(electron.deltaEtaSuperClusterTrackAtVtx)   <  0.007;
        bool tightDEtaInCut       = fabs(electron.deltaEtaSuperClusterTrackAtVtx)   <  0.004;
        bool looseDPhiInCut       = fabs(electron.deltaPhiSuperClusterTrackAtVtx)   <  0.15;
        bool tightDPhiInCut       = fabs(electron.deltaPhiSuperClusterTrackAtVtx)   <  0.03;
        bool looseDzVtxCut        = fabs(electron.gsfTrack->dzPv(primVtx.position)) <  0.2;
        bool tightDzVtxCut        = fabs(electron.gsfTrack->dzPv(primVtx.position)) <  0.1;
        bool looseIsoCut          = pfCombinedIso                                   <  0.15;
        bool tightIsoCut          = pfCombinedIso                                   <  0.1;
        bool looseNMissingHitsCut = electron.nMissingHits                           <= 2;
        bool tightNMissingHitsCut = electron.nMissingHits                           <= 1;

        bool looseCut = looseDEtaInCut && looseDPhiInCut && looseDzVtxCut &&
                        looseIsoCut    && looseNMissingHitsCut;
        bool tightCut = tightDEtaInCut && tightDPhiInCut && tightDzVtxCut &&
                        tightIsoCut    && tightNMissingHitsCut;

        if (looseCut) looseElectrons.push_back(& electron);
        if (tightCut) tightElectrons.push_back(& electron);
      }

      ////////// SORT SELECTED ELECTRONS //////////
      std::sort(looseElectrons.begin(), looseElectrons.end(), PtGreater<susy::Electron>);
      std::sort(tightElectrons.begin(), tightElectrons.end(), PtGreater<susy::Electron>);

      // Leading electron must pass higher pT cut
      if (looseElectrons.size() >= 1 && looseElectrons.at(0)->momentum.Et() < electron1_etCut) looseElectrons.clear();
      if (tightElectrons.size() >= 1 && tightElectrons.at(0)->momentum.Et() < electron1_etCut) tightElectrons.clear();

      ////////// SELECT PHOTONS //////////
      if (printLevel > 1) out << "Find photons in the event" << std::endl;

      std::vector<susy::Photon const*> loosePhotons;
      std::vector<susy::Photon const*> mediumPhotons;
      std::vector<susy::Photon const*> tightPhotons;

      susy::PhotonCollection const& photons(event.photons["photons"]);

      for (size_t iP(0); iP != photons.size(); ++iP) {
        susy::Photon const& photon(photons[iP]);

        if (photon.momentum.Et() < photon2_etCut) continue;
        if (!photon.isEB())                       continue;
        if (!photon.passelectronveto)             continue;
        if (photon.hadTowOverEm  > 0.05)          continue;

        Float_t effA[3];

        photonEffectiveAreas(photon.momentum.Eta(), effA);

        Float_t pfChIso = photon.chargedHadronIso - event.rho25 * effA[0];
        Float_t pfNhIso = photon.neutralHadronIso - event.rho25 * effA[1];
        Float_t pfPhIso = photon.photonIso        - event.rho25 * effA[2];

        bool looseSigmaIEtaIEtaCut    = (photon.sigmaIetaIeta < 0.012);
        bool tightSigmaIEtaIEtaCut    = (photon.sigmaIetaIeta < 0.011);
        bool loosePfChargedHadronIso  = (pfChIso              < 2.6);
        bool mediumPfChargedHadronIso = (pfChIso              < 1.5);
        bool tightPfChargedHadronIso  = (pfChIso              < 0.7);
        bool loosePfNeutralHadronIso  = (pfNhIso              < 3.5 + 0.04  * photon.momentum.Et());
        bool mediumPfNeutralHadronIso = (pfNhIso              < 1.0 + 0.04  * photon.momentum.Et());
        bool tightPfNeutralHadronIso  = (pfNhIso              < 0.4 + 0.04  * photon.momentum.Et());
        bool loosePfPhotonIso         = (pfPhIso              < 1.3 + 0.005 * photon.momentum.Et());
        bool mediumPfPhotonIso        = (pfPhIso              < 0.7 + 0.005 * photon.momentum.Et());
        bool tightPfPhotonIso         = (pfPhIso              < 0.5 + 0.005 * photon.momentum.Et());

        bool looseCut  = looseSigmaIEtaIEtaCut && loosePfChargedHadronIso  && loosePfNeutralHadronIso  && loosePfPhotonIso;
        bool mediumCut = tightSigmaIEtaIEtaCut && mediumPfChargedHadronIso && mediumPfNeutralHadronIso && mediumPfPhotonIso;
        bool tightCut  = tightSigmaIEtaIEtaCut && tightPfChargedHadronIso  && tightPfNeutralHadronIso  && tightPfPhotonIso;

        if (looseCut)  loosePhotons.push_back (& photon);
        if (mediumCut) mediumPhotons.push_back(& photon);
        if (tightCut)  tightPhotons.push_back (& photon);
      }

      ////////// SORT SELECTED PHOTONS //////////
      std::sort(loosePhotons.begin(),  loosePhotons.end(),  PtGreater<susy::Photon>);
      std::sort(mediumPhotons.begin(), mediumPhotons.end(), PtGreater<susy::Photon>);
      std::sort(tightPhotons.begin(),  tightPhotons.end(),  PtGreater<susy::Photon>);

      // Leading photon must pass higher ET cut
      if (loosePhotons.size()  >= 1 && loosePhotons.at (0)->momentum.Et() < photon1_etCut) loosePhotons.clear();
      if (mediumPhotons.size() >= 1 && mediumPhotons.at(0)->momentum.Et() < photon1_etCut) mediumPhotons.clear();
      if (tightPhotons.size()  >= 1 && tightPhotons.at (0)->momentum.Et() < photon1_etCut) tightPhotons.clear();

      ////////// SELECT JETS //////////
      if (printLevel > 1) out << "Find PF jets in the event" << std::endl;

      std::vector<susy::PFJet*> goodPfJets;

      susy::PFJetCollection& pfJets(event.pfJets["ak5"]);

      for (size_t iJ(0); iJ != pfJets.size(); ++iJ) {
        susy::PFJet& jet(pfJets[iJ]);

        bool sameAsTightMuon     = false;
        bool sameAsTightElectron = false;
        bool sameAsMediumPhoton  = false;

        for (size_t i = 0; i < tightMuons.size(); i++) {
          if (isSameObject2(jet.momentum, tightMuons.at(i)->momentum, 0.5)) {
            sameAsTightMuon = true;
            break;
          }
        }

        for (size_t i = 0; i < tightElectrons.size(); i++) {
          if (isSameObject2(jet.momentum, tightElectrons.at(i)->momentum, 0.5)) {
            sameAsTightElectron = true;
            break;
          }
        }

        for (size_t i = 0; i < mediumPhotons.size(); i++) {
          if (isSameObject2(jet.momentum, mediumPhotons.at(i)->momentum, 0.5)) {
            sameAsMediumPhoton = true;
            break;
          }
        }

        if (sameAsTightMuon     ||
            sameAsTightElectron ||
            sameAsMediumPhoton)                                  continue;
        if (fabs(jet.momentum.Eta())                   >  2.4)   continue;
        if (jet.chargedHadronEnergy / jet.momentum.E() <  1.e-6) continue;
        if (jet.neutralHadronEnergy / jet.momentum.E() >  0.99)  continue;
        if (jet.chargedEmEnergy     / jet.momentum.E() >  0.99)  continue;
        if (jet.neutralEmEnergy     / jet.momentum.E() >  0.99)  continue;
        if (jet.nConstituents                          <= 1)     continue;
        if (jet.chargedMultiplicity                    == 0)     continue;
        if (jet.passPuJetIdLoose(susy::kPUJetIdFull)   != 1)     continue;

        Float_t jecScale, jecScaleUncertainty;

        if (useCustomJEC) {
          jetCorrection->setJetPt (jet.momentum.Pt());
          jetCorrection->setJetEta(jet.momentum.Eta());
          jetCorrection->setJetA  (jet.jetArea);
          jetCorrection->setRho   (event.rho);

          jecScale = jetCorrection->getCorrection();

          jecUncertainty->setJetPt (jet.momentum.Pt());
          jecUncertainty->setJetEta(jet.momentum.Eta());

          try {
            jecScaleUncertainty = jecUncertainty->getUncertainty(true);
          }
          catch (std::exception& e) {
            std::cerr << "Cannot get uncertainty for jet pT = " << jet.momentum.Pt() << " eta = " << jet.momentum.Eta() << ". Setting to -1.0" << std::endl;
            jecScaleUncertainty = -1.0;
          }
        }
        else {
          jecScale            = jet.jecScaleFactors.find("L1FastL2L3")->second;
          jecScaleUncertainty = jet.jecUncertainty;
        }

        TLorentzVector corrP4(jecScale * jet.momentum);

        if (corrP4.Pt() > jet_ptCut) {
          goodPfJets.push_back(& jet);

          goodPfJets.back()->momentum = corrP4 * (1.0 + jecSystematic * jecScaleUncertainty);

	  jet_unc->push_back(jecScaleUncertainty);
	  // cout << "jecScaleUncertainty: "<< jecScaleUncertainty << endl; 
	  //cout << corrP4.Perp() << " Eta: " << corrP4.Eta() << endl; 
	  //corrP4*=(1+jecScaleUncertainty);
	  //cout << corrP4.Perp() << " Eta: " << corrP4.Eta() << endl;
        }
      }

      

      ////////// W,Z,t pt ///////////////
      TLorentzVector pGentot;
      pGentot.SetPxPyPzE(0,0,0,0); 
		Float_t St_gen=0;
		Float_t St_Lept_gen=0;
		Int_t Ngamma=0;
		
		for(size_t i=0; i< event.genParticles.size(); i++){
			if(event.genParticles.at(i).pdgId==22) Ngamma++;
		}
		
        float Wtopi=1;
        float WtopiP=1;
        float WtopiM=1;

        
      for(size_t i=0; i< event.genParticles.size(); i++){
	if(event.genParticles.at(i).status==1){
	  pGentot+=event.genParticles.at(i).momentum;
	}
		  if(event.genParticles.at(i).pdgId==23) hName["h_Zmass"]->Fill(event.genParticles.at(i).momentum.M()); 

		  if(Ngamma==0){
		  /*
		  cout << "pdgId: " << event.genParticles.at(i).pdgId << " "; 
		  cout << " E: " << event.genParticles.at(i).momentum.E();
		  cout << " pt: " << event.genParticles.at(i).momentum.Perp(); 
		  cout << " eta: " << event.genParticles.at(i).momentum.Eta();
		  cout << " phi: " << event.genParticles.at(i).momentum.Phi();
		  */
			  if(event.genParticles.at(i).pdgId==23) hName["h_Zmass_nGamma0"]->Fill(event.genParticles.at(i).momentum.M()); 
			  //cout << " "<< event.genParticles.at(i).momentum.M(); 
		  
		 // cout << endl; 
		  }
		  
		  
		if(event.genParticles.at(i).status==1 ){
			/*
				cout << "pdgId: " << event.genParticles.at(i).pdgId << " "; 
				cout << " E: " << event.genParticles.at(i).momentum.E();
				cout << " pt: " << event.genParticles.at(i).momentum.Perp(); 
				cout << " eta: " << event.genParticles.at(i).momentum.Eta();
				cout << " phi: " << event.genParticles.at(i).momentum.Phi();
				cout << endl; 
			*/
			if(abs(event.genParticles.at(i).pdgId)>=11 && abs(event.genParticles.at(i).pdgId)<=16 )  {
				St_Lept_gen+=event.genParticles.at(i).momentum.Perp();
			}
			
			St_gen+=event.genParticles.at(i).momentum.Perp();
			
	  }
	
          if(abs(event.genParticles.at(i).pdgId)==14){
              TLorentzVector Ptau=event.genParticles.at(i).momentum;
            //  cout << "tau pt: " << Ptau.Perp() << endl;
          }
	
	if(abs(event.genParticles.at(i).pdgId)==24){
	  //W 4 vector
	  TLorentzVector Pw=event.genParticles.at(i).momentum; 
	  w_px->push_back(Pw.Px());
	  w_py->push_back(Pw.Py());
	  w_pz->push_back(Pw.Pz());
	  w_e->push_back(Pw.E()); 
	  //	  cout << "W pt: " << Pw.Perp() << " eta: " << Pw.Eta() << endl; 
	}
	if(abs(event.genParticles.at(i).pdgId)==23){
	  //Z 4 vector
	  TLorentzVector Pz=event.genParticles.at(i).momentum;
		  z_px->push_back(Pz.Px());
          z_py->push_back(Pz.Py());
          z_pz->push_back(Pz.Pz());
          z_e->push_back(Pz.E());

        }
	if(abs(event.genParticles.at(i).pdgId)==6){
	  //top 4 vector
          TLorentzVector Ptop=event.genParticles.at(i).momentum;
        float x=Ptop.Perp();
        if(x>400)x=400;
        float w=TMath::Exp(0.156-0.00137*x);

        Wtopi*=w;
     //   std::cout << "top pt: " << Ptop.Perp() << " weight: " << TMath::Exp(0.156-0.00137*x) << std::endl;
        WtopiP*=1.2*Wtopi;
        WtopiM*=0.8*Wtopi;

		  top_px->push_back(Ptop.Px());
          top_py->push_back(Ptop.Py());
          top_pz->push_back(Ptop.Pz());
          top_e->push_back(Ptop.E());
        }


      }//end gen particles loop 
        SumW_Top+=TMath::Sqrt(Wtopi);
        SumW_TopP+=TMath::Sqrt(WtopiP);
        SumW_TopM+=TMath::Sqrt(WtopiM);

		Float_t Ht_gen=St_gen-St_Lept_gen;
		
		gen_st=St_gen;//for filling trees 
		gen_ht=Ht_gen;

		//cout << "GEN St: " << St_gen << " Ht: " << Ht_gen << endl; 
		
      //      cout << "pGenTot: " << pGentot.Perp() << endl; 

      ////////// SORT SELECTED JETS //////////
      std::sort(goodPfJets.begin(), goodPfJets.end(), PtGreater<susy::PFJet>);
      /*
      cout << "reconstructed jets: " << endl; 
      for (size_t i = 0; i < goodPfJets.size(); i++) {
	cout << " E: " << goodPfJets.at(i)->momentum.E();
        cout << " pt: " << goodPfJets.at(i)->momentum.Pt();
	cout << " eta: " <<  goodPfJets.at(i)->momentum.Eta();
	cout << " phi: " << goodPfJets.at(i)->momentum.Phi(); 
	cout << endl; 
      }
      cout << endl; 
      */
      ////////// MET //////////
      TVector2 const& metV(event.metMap["pfType01CorrectedMet"].mEt);

      ////////// INCREMENT COUNTERS //////////
      if (tightMuons.size()     >= 1) nCnt[2]++;
      if (tightElectrons.size() >= 1) nCnt[3]++;
      if (mediumPhotons.size()  >= 1) nCnt[4]++;

      ////////// FILL HISTOGRAMS AND NTUPLE //////////
      if (looseMuons.size() + looseElectrons.size() + loosePhotons.size() == 0) continue;

      // Number of events with either a loose lepton or a loose photon
      nCnt[5]++;

      if (goodPfJets.size() == 0) nCnt [6]++;
      if (goodPfJets.size() == 1) nCnt [7]++;
      if (goodPfJets.size() == 2) nCnt [8]++;
      if (goodPfJets.size() == 3) nCnt [9]++;
      if (goodPfJets.size() == 4) nCnt[10]++;
      if (goodPfJets.size() == 5) nCnt[11]++;
      if (goodPfJets.size() == 6) nCnt[12]++;
      if (goodPfJets.size() == 7) nCnt[13]++;

		int npv = 0;
		//cout << "isRealData: " << event.isRealData << endl; 
		if(!event.isRealData){
			for (susy::PUSummaryInfoCollection::const_iterator pu = event.pu.begin(); pu != event.pu.end(); ++pu){
				int bx = (int)(*pu).BX;
				if (bx == 0){
					npv = (*pu).trueNumInteractions;
					continue;
				}  
			}
			//cout << "npv = " << npv << " " << endl; 
			NumInteractions=npv;
			std::vector<Float_t> PUWeights = GetPileUpWeight(npv);
			puWt_nom  =  PUWeights[0];
			puWt_up   =  PUWeights[1];
			puWt_down =  PUWeights[2];
		}
		else{
			NumInteractions=0;
			puWt_nom  =  1.;
			puWt_up   =  1.;
			puWt_down =  1.;
		}
		//cout << " npv = " << npv <<", nominal PU weight = " << puWt_nom << endl;
		
		
      runNo        = event.runNumber;
      lumiNo       = event.luminosityBlockNumber;
      eventNo      = event.eventNumber;
      squark_m     = event.gridParams["pole_squark_m"];
      neutralino_m = event.gridParams["pole_gaugino_m"];
      vertices_n   = event.vertices.size();
      muons_n      = tightMuons.size();
      electrons_n  = tightElectrons.size();
      photons_n    = mediumPhotons.size();
      jets_n       = goodPfJets.size();
      met_et       = metV.Mod();
      met_phi      = metV.Phi();

      if (datasetName.Contains("stealthT2"))
        evtWt = SetSignalEventWeight(datasetName, lumiCalc, event.gridParams["pole_squark_m"]);
      else if (nScaledEvents > 0.0)
        evtWt = nScaledEvents / (Float_t)nEntries;
      else
        evtWt = 1.0;

      Float_t tempSt = 0.0;
		Float_t tempHt=0.0;
//loose muons - do not fill St
		for (size_t i = 0; i < looseMuons.size(); i++) {
			
			loose_muon_e->push_back     (looseMuons.at(i)->momentum.E());
			loose_muon_px->push_back    (looseMuons.at(i)->momentum.Px());
			loose_muon_py->push_back    (looseMuons.at(i)->momentum.Py());
			loose_muon_pz->push_back    (looseMuons.at(i)->momentum.Pz());
			loose_muon_charge->push_back(looseMuons.at(i)->bestTrack->charge);
		}	
		
//tight muons		
      for (size_t i = 0; i < tightMuons.size(); i++) {
        tempSt += tightMuons.at(i)->momentum.Pt();

        muon_e->push_back     (tightMuons.at(i)->momentum.E());
        muon_px->push_back    (tightMuons.at(i)->momentum.Px());
        muon_py->push_back    (tightMuons.at(i)->momentum.Py());
        muon_pz->push_back    (tightMuons.at(i)->momentum.Pz());
        muon_charge->push_back(tightMuons.at(i)->bestTrack->charge);
      }
		
//loose electrons - do not fill st
	  for (size_t i = 0; i < looseElectrons.size(); i++) {
		
			loose_electron_e->push_back     (looseElectrons.at(i)->momentum.E());
			loose_electron_px->push_back    (looseElectrons.at(i)->momentum.Px());
			loose_electron_py->push_back    (looseElectrons.at(i)->momentum.Py());
			loose_electron_pz->push_back    (looseElectrons.at(i)->momentum.Pz());
			loose_electron_charge->push_back(looseElectrons.at(i)->gsfTrack->charge);
		}
		
//tight electrons - use for st definition		
      for (size_t i = 0; i < tightElectrons.size(); i++) {
        tempSt += tightElectrons.at(i)->momentum.Et();

        electron_e->push_back     (tightElectrons.at(i)->momentum.E());
        electron_px->push_back    (tightElectrons.at(i)->momentum.Px());
        electron_py->push_back    (tightElectrons.at(i)->momentum.Py());
        electron_pz->push_back    (tightElectrons.at(i)->momentum.Pz());
        electron_charge->push_back(tightElectrons.at(i)->gsfTrack->charge);
      }

      for (size_t i = 0; i < mediumPhotons.size(); i++) {
        tempSt += mediumPhotons.at(i)->momentum.Et();

        photon_e->push_back (mediumPhotons.at(i)->momentum.E());
        photon_px->push_back(mediumPhotons.at(i)->momentum.Px());
        photon_py->push_back(mediumPhotons.at(i)->momentum.Py());
        photon_pz->push_back(mediumPhotons.at(i)->momentum.Pz());
      }

      for (size_t i = 0; i < goodPfJets.size(); i++) {
        tempSt += goodPfJets.at(i)->momentum.Pt();
		tempHt += goodPfJets.at(i)->momentum.Pt();
		  
        jet_e->push_back   (goodPfJets.at(i)->momentum.E());
        jet_px->push_back  (goodPfJets.at(i)->momentum.Px());
        jet_py->push_back  (goodPfJets.at(i)->momentum.Py());
        jet_pz->push_back  (goodPfJets.at(i)->momentum.Pz());
		  
		/*
		 https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP#B_tagging_Operating_Points_for_5 
		 */
		jet_bTagL->push_back(goodPfJets.at(i)->bTagDiscriminators[susy::kCSV] > 0.244); // loose csv working point b-tag
        jet_bTagM->push_back(goodPfJets.at(i)->bTagDiscriminators[susy::kCSV] > 0.679); //medium csv working point b-tag
		jet_bTagT->push_back(goodPfJets.at(i)->bTagDiscriminators[susy::kCSV] > 0.898); //tight csv working point b-tag
		  jet_algFlavor->push_back(goodPfJets.at(i)->algDefFlavour); 
		  jet_phyFlavor->push_back(goodPfJets.at(i)->phyDefFlavour); 

		//cout << "alg flavor: " << goodPfJets.at(i)->algDefFlavour << endl; 
		//  bool TAG=goodPfJets.at(i)->bTagDiscriminators[susy::kCSV] > 0.679;
		//  cout << "b-tagged: "<< TAG << endl; 
		//cout << "phy flavor: " << goodPfJets.at(i)->phyDefFlavour << endl; 
      }

      if (metV.Mod() > met_etCut)
        tempSt += metV.Mod();

      st = tempSt;
/*
		if(muon_e->size()>=1 && electron_e->size()>=1 && loose_electron_e->size()!=electron_e->size() && loose_muon_e->size()!=muon_e->size()){
			cout << "Nmu: "<< muon_e->size() << " NEl: " << electron_e->size() << endl;  
			int nGenMu=0;
			int nGenEl=0;
			for(size_t i=0; i< event.genParticles.size(); i++){
				
				if(event.genParticles.at(i).status==1){
					cout << "pdgId: " << event.genParticles.at(i).pdgId << " "; 
					cout << " E: " << event.genParticles.at(i).momentum.E();
					cout << " pt: " << event.genParticles.at(i).momentum.Perp(); 
					cout << " eta: " << event.genParticles.at(i).momentum.Eta();
					cout << " phi: " << event.genParticles.at(i).momentum.Phi();
					cout << endl; 
					if(abs(event.genParticles.at(i).pdgId)==11) nGenEl++;
					if(abs(event.genParticles.at(i).pdgId)==13) nGenMu++;

				}
			}
			cout << "nGenEl: " << nGenEl << " nGenMu: " << nGenMu << endl; 
			cout << endl; 
			
		}			
 */
		//cout << "RECO st: " << st << " ht: " << tempHt << endl << endl; 	
		
      //      cout << "jet px: "<< jet_px->size() << " jet_unc: " << jet_unc->size() << endl; 
      if(jet_px->size()>=1) tree->Fill();

      ////////// FILL SKIMS //////////
      if(copyEvents){
        if(tightPhotons.size() >= 2) ggTree->Fill();
        if(loosePhotons.size() >= 2) ffTree->Fill();
      }
    }

    ////////// END OF EVENT LOOP //////////
	  hName["h_Nevents"]->SetBinContent(1,nCnt[0]); 
    out << " -------------------- Job Summary -------------------- "                                                                  << std::endl;
    out << " Total events                                                 : " << nCnt [0]                                             << std::endl;
    out << " Sum Top weights                                                : " << SumW_Top << " +/ " << SumW_TopP << " " << SumW_TopM <<  std::endl;

    if (nCnt[0] >= 1) {
      out << " passed preselection                                        : " << nCnt [1] << " (" << nCnt [1] / float(nCnt[0]) << ")" << std::endl;
      out << " tightMuons     >= 1                                        : " << nCnt [2] << " (" << nCnt [2] / float(nCnt[1]) << ")" << std::endl;
      out << " tightElectrons >= 1                                        : " << nCnt [3] << " (" << nCnt [3] / float(nCnt[1]) << ")" << std::endl;
      out << " mediumPhotons  >= 1                                        : " << nCnt [4] << " (" << nCnt [4] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1                    : " << nCnt [5] << " (" << nCnt [5] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 0 : " << nCnt [6] << " (" << nCnt [6] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 1 : " << nCnt [7] << " (" << nCnt [7] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 2 : " << nCnt [8] << " (" << nCnt [8] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 3 : " << nCnt [9] << " (" << nCnt [9] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 4 : " << nCnt[10] << " (" << nCnt[10] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 5 : " << nCnt[11] << " (" << nCnt[11] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 6 : " << nCnt[12] << " (" << nCnt[12] / float(nCnt[1]) << ")" << std::endl;
      out << " looseLeptons   >= 1 || loosePhoton >= 1 && goodPfJets == 7 : " << nCnt[13] << " (" << nCnt[13] / float(nCnt[1]) << ")" << std::endl;
    }

    if(printLevel > 0) out << "Save outputs" << std::endl;

    ////////// END OF EVENT LOOP //////////

    fout->cd();
    fout->Write();

    if(copyEvents){
      event.releaseTree(*ggTree);
      event.releaseTree(*ffTree);

      ggFile->cd();
      ggFile->Write();

      ffFile->cd();
      ffFile->Write();
    }
  }
  catch(std::exception& e){
    std::cerr << e.what() << std::endl;

    event.releaseTrees();
  }

  delete jetCorrection;
  delete jecUncertainty;

  delete ggFile;
  delete ffFile;
  delete fout;
}

