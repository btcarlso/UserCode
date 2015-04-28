// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.h
//
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.h,v 1.6 2013/06/21 20:02:11 weinberg Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TChain.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TArrayI.h>

#include <map>
#include <set>
#include <fstream>

#include "SusyEvent.h"

class SusyEventAnalyzer {
public:
  SusyEventAnalyzer(TTree&);
  virtual ~SusyEventAnalyzer();

  /* Main analyzer function to be defined */
  virtual void Run();

  /* parameter configuration functions */
  void IncludeAJson(TString const&);
  void SetOutputDirectory(TString const& v) { outputDirectoryName = v; }
  void SetOutput(TString const& v) { outputName = v; }
  void SetLogFile(TString const& v) { logFileName = v;}
  void SetPrintInterval(int v) { printInterval = v; }
  void SetPrintLevel(int v) { printLevel = v; }
  void SetProcessNEvents(int v) { processNEvents = v; }
  void AddHltName(TString const& v) { hltNames.push_back(v + "_v*"); }
  void CopyEvents(bool v) { copyEvents = v; }

  void SetNScaledEvents(Float_t lumi, Float_t xSec) { lumiCalc = lumi; nScaledEvents = lumi * xSec; }
  void SetDatasetName  (TString const& v)           { datasetName = v; }
  void SetCutComplement(TString const& v)           { cutComplement = v; SetObjectCuts(cutComplement); }
  void SetJecUncert    (Float_t jec)                { jecSystematic = jec; }

  void    SetObjectCuts       (TString);

  std::vector<Float_t> GetPileUpWeight(Int_t);
	
  Float_t SetSusyXSec         (Float_t);
  Float_t SetSignalEventWeight(TString, Float_t, Float_t);

  void createHistogram(const char*, const char*, const char*, const char*, Int_t, Double_t, Double_t);

protected:
  bool IsGoodLumi(UInt_t, UInt_t) const;
  bool PassTriggers() const;

  /* container of all event data */
  susy::Event event;
  /* input tree */
  TTree *fTree;
  /* directory of the output file */
  TString outputDirectoryName;
  /* suffix of the output file */
  TString outputName;
  /* log file name */
  TString logFileName;
  /* verbosity - 0 => no printout, 1 => print function control flow, 2 => print event processing flow, 3 => print event dump */
  int printLevel;
  /* print frequency */
  unsigned printInterval;
  /* maximum number of events */
  int processNEvents;
  /* HLT path names */
  std::vector<TString> hltNames;
  /* switch for saving skims */
  bool copyEvents;
  /* good lumi list */
  std::map<unsigned, std::set<unsigned> > goodLumiList;
  mutable std::pair<unsigned, unsigned> currentLumi;
  mutable bool currentLumiIsGood;

  std::map<TString, TH1F*> hName;  // Map for histograms
  Float_t lumiCalc;                // Luminosity used for scaling output
  Float_t nScaledEvents;           // Number of events to be used for scaling output histograms
  TString datasetName;
  TString cutComplement;
  Float_t jecSystematic;           // -1.0, 0.0, +1.0 for JEC down, nominal, JEC up

  // Cut variables
  Float_t muon1_ptCut,   muon2_ptCut,   electron1_etCut, electron2_etCut;
  Float_t photon1_etCut, photon2_etCut, jet_ptCut,       met_etCut;

  // Flat ntuple variables
  Int_t   runNo,      lumiNo;
  ULong_t eventNo;
  Float_t squark_m,   neutralino_m;
  Int_t   vertices_n, muons_n, electrons_n, photons_n, jets_n;
  Float_t evtWt;
  Float_t met_et,     met_phi, st;
  Float_t gen_st, gen_ht;
	
  Float_t puWt_nom,		puWt_up,	puWt_down;
  Int_t NumInteractions;

  std::vector<float>* muon_e; // tight muons 
  std::vector<float>* muon_px;
  std::vector<float>* muon_py;
  std::vector<float>* muon_pz;
  std::vector<float>* muon_charge;
	
  std::vector<float>* loose_muon_e; // loose muons 
  std::vector<float>* loose_muon_px;
  std::vector<float>* loose_muon_py;
  std::vector<float>* loose_muon_pz;
  std::vector<float>* loose_muon_charge;	
	
  std::vector<float>* electron_e; //tight electrons 
  std::vector<float>* electron_px;
  std::vector<float>* electron_py;
  std::vector<float>* electron_pz;
  std::vector<float>* electron_charge;
	
  std::vector<float>* loose_electron_e; //tight electrons 
  std::vector<float>* loose_electron_px;
  std::vector<float>* loose_electron_py;
  std::vector<float>* loose_electron_pz;
  std::vector<float>* loose_electron_charge;		
	
  std::vector<float>* photon_e;
  std::vector<float>* photon_px;
  std::vector<float>* photon_py;
  std::vector<float>* photon_pz;
  std::vector<float>* jet_e;
  std::vector<float>* jet_px;
  std::vector<float>* jet_py;
  std::vector<float>* jet_pz;
  std::vector<bool>*  jet_bTagL;
  std::vector<bool>*  jet_bTagM;
  std::vector<bool>*  jet_bTagT;
  std::vector<int>*  jet_algFlavor;
  std::vector<int>*  jet_phyFlavor;

  std::vector<float>* jet_unc;

  std::vector<float>* w_px;
  std::vector<float>* w_py;
  std::vector<float>* w_pz;
  std::vector<float>* w_e;

  std::vector<float>* z_px;
  std::vector<float>* z_py;
  std::vector<float>* z_pz;
  std::vector<float>* z_e;

  std::vector<float>* top_px;
  std::vector<float>* top_py;
  std::vector<float>* top_pz;
  std::vector<float>* top_e;
};

SusyEventAnalyzer::SusyEventAnalyzer(TTree& tree) :
  event(),
  fTree(&tree),
  outputDirectoryName("."),
  outputName("analysis"),
  logFileName(outputName + ".log"),
  printLevel(0),
  printInterval(1000),
  processNEvents(-1),
  hltNames(),
  copyEvents(false),
  goodLumiList(),
  currentLumi(0, 0),
  currentLumiIsGood(true),
  lumiCalc(1.0),
  nScaledEvents(-1.0),
  datasetName(""),
  cutComplement(""),
  jecSystematic(0.0)
{
  event.setInput(tree);

  muon1_ptCut     = 15.0;
  muon2_ptCut     = 15.0;
  electron1_etCut = 15.0;
  electron2_etCut = 15.0;
  photon1_etCut   = 15.0;
  photon2_etCut   = 15.0;
  jet_ptCut       = 30.0;
  met_etCut       = 15.0;
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
}

void
SusyEventAnalyzer::IncludeAJson(TString const& _fileName)
{
  if(_fileName == "") return;

  std::ifstream inputFile(_fileName);
  if(!inputFile.is_open()){
    std::cerr << "Cannot open JSON file " << _fileName << std::endl;
    return;
  }

  std::string line;
  TString jsonText;
  while(true){
    std::getline(inputFile, line);
    if(!inputFile.good()) break;
    jsonText += line;
  }
  inputFile.close();

  TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
  TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

  TArrayI positions(2);
  positions[1] = 0;
  while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
    TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
    TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

    unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
    std::set<unsigned>& lumis(goodLumiList[run]);

    TArrayI lumiPos(2);
    lumiPos[1] = 0;
    while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
      TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
      int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
      int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
      for(int lumi(begin); lumi <= end; ++lumi)
        lumis.insert(lumi);
    }
  }
}



void
SusyEventAnalyzer::SetObjectCuts(TString cuts)
{
  if (cuts == "singleMu")
    muon1_ptCut = 30.0;
  else if (cuts == "singleElectron")
    electron1_etCut = 30.0;
  else if (cuts == "doublePhoton") {
    photon1_etCut = 40.0;
    photon2_etCut = 25.0;
  }
  else if (cuts == "photonHad")
    photon1_etCut = 75.0;
}

std::vector<Float_t>
SusyEventAnalyzer::GetPileUpWeight(Int_t npv)
{
	std::vector < float > NominalWeights;   
	Double_t Nominal_S10[60] = {
		0.2681,
		0.4996,
		7.5583,
		0.3367,
		0.2822,
		0.6834,
		0.5136,
		0.5068,
		0.6935,
		1.0129,
		1.4198,
		1.7409,
		1.7484,
		1.5471,
		1.3143,
		1.1420,
		1.0512,
		1.0219,
		1.0392,
		1.0838,
		1.1256,
		1.1524,
		1.1685,
		1.1734,
		1.1584,
		1.1165,
		1.0473,
		0.9564,
		0.8501,
		0.7353,
		0.6188,
		0.5068,
		0.4033,
		0.3116,
		0.2335,
		0.1701,
		0.1217,
		0.0876,
		0.0654,
		0.0527,
		0.0470,
		0.0464,
		0.0495,
		0.0555,
		0.0641,
		0.0752,
		0.0892,
		0.1065,
		0.1277,
		0.0527,
		0.1851,
		0.2234,
		0.2701,
		0.3268,
		0.3954,
		0.4785,
		0.5789,
		0.7004,
		0.8473,
		2.1570
	};
	
	std::vector < float > PileupUpWeights;   
	Double_t PileupUp_S10[60] = {
		0.2523,
		0.2739,
		7.0764,
		0.4072,
		0.2458,
		0.3487,
		0.3141,
		0.3367,
		0.4788,
		0.7302,
		1.0466,
		1.3589,
		1.4819,
		1.3804,
		1.1956,
		1.0520,
		0.9845,
		0.9705,
		0.9936,
		1.0442,
		1.0945,
		1.1279,
		1.1493,
		1.1672,
		1.1816,
		1.1850,
		1.1683,
		1.1271,
		1.0621,
		0.9780,
		0.8797,
		0.7726,
		0.6621,
		0.5540,
		0.4523,
		0.3592,
		0.2778,
		0.2101,
		0.1570,
		0.1183,
		0.0928,
		0.0782,
		0.0724,
		0.0735,
		0.0804,
		0.0922,
		0.1089,
		0.1309,
		0.1591,
		0.1183,
		0.2398,
		0.2962,
		0.3670,
		0.4562,
		0.5678,
		0.7081,
		0.8841,
		1.1051,
		1.3818,
		3.6347
	};
	
	std::vector < float > PileupDownWeights;   
	Double_t PileupDown_S10[60] = {
		0.2859,
		1.1471,
		7.5903,
		0.3720,
		0.3833,
		1.3167,
		0.8238,
		0.7683,
		1.0014,
		1.4116,
		1.9071,
		2.1424,
		2.0002,
		1.7165,
		1.4381,
		1.2281,
		1.1117,
		1.0718,
		1.0812,
		1.1162,
		1.1511,
		1.1718,
		1.1713,
		1.1406,
		1.0744,
		0.9781,
		0.8625,
		0.7371,
		0.6102,
		0.4894,
		0.3804,
		0.2860,
		0.2075,
		0.1456,
		0.0998,
		0.0683,
		0.0484,
		0.0371,
		0.0319,
		0.0305,
		0.0318,
		0.0348,
		0.0392,
		0.0449,
		0.0519,
		0.0603,
		0.0703,
		0.0821,
		0.0961,
		0.0305,
		0.1319,
		0.1544,
		0.1806,
		0.2112,
		0.2466,
		0.2877,
		0.3357,
		0.3924,
		0.4602,
		1.1403
	};
	
	for (int i = 0; i < 60; ++i){
		NominalWeights.push_back(Nominal_S10[i]);
		PileupUpWeights.push_back(PileupUp_S10[i]);
		PileupDownWeights.push_back(PileupDown_S10[i]);
	}
	
	std::vector < float > PUWeights(3);
	PUWeights[0] =  NominalWeights[npv];
	PUWeights[1] =  PileupUpWeights[npv]; 
	PUWeights[2] =  PileupDownWeights[npv];
	
	return PUWeights;
}

bool
SusyEventAnalyzer::IsGoodLumi(UInt_t run, UInt_t lumi) const
{
  if(goodLumiList.size() == 0) return true;
  if(run == currentLumi.first && lumi == currentLumi.second) return currentLumiIsGood;
  currentLumi.first = run;
  currentLumi.second = lumi;
  currentLumiIsGood = false;

  std::map<unsigned, std::set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
  if(rItr != goodLumiList.end()){
    std::set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
    if(lItr != rItr->second.end()) currentLumiIsGood = true;
  }

  return currentLumiIsGood;
}

bool
SusyEventAnalyzer::PassTriggers() const
{
  unsigned nT(hltNames.size());
  if(nT == 0) return true;

  for(unsigned iT(0); iT != nT; ++iT)
    if(event.hltMap.pass(hltNames[iT])) return true;

  return false;
}

#endif

