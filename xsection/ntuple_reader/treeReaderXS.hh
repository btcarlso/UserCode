#ifndef TREEREADERXS_H
#define TREEREADERXS_H
                                                                       
#include <iostream>
  
#include <TROOT.h>
#include <TSystem.h>
#include "TUnixSystem.h"
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TStopwatch.h>

//#include "/afs/cern.ch/work/b/bcarlson/public/common_code/bins_new.card"

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

using namespace std;

#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567


class treeReaderXS {
public:
  treeReaderXS(TChain *tree, TString evtClassName);
  virtual      ~treeReaderXS();
  virtual void init(TString evtClassName);

  virtual void openHistFile(TString filename);
  virtual void openCopyFile1(TString filename); 
  virtual void openCopyFile2(TString filename);
  virtual void openCopyFile3(TString filename);
  virtual void openDiscardedFile(TString discarded_filename);
  virtual void initialize_counters();
  virtual void print_numbers(int nevents, int h);
  virtual void closeHistFile();
  virtual void closeCopyFile1();
  virtual void closeCopyFile2();
  virtual void closeCopyFile3();
  virtual void closeDiscardedFile();
  virtual void normalize();
  virtual void bookHist();
  virtual void readCuts(TString filename, int dump = 1);
  virtual bool check_evt(); 
  virtual void startAnalysis();
  virtual int  loop(int nevents = 1, int start = -1, bool s_a=1);
  virtual int  get_mem(long* vmrss_kb, long* vmsize_kb);	
  virtual void eventProcessing();
  // virtual bool eventProcessing();
  
  virtual bool goodRun();
  virtual void fillHist();

protected:
  //stuf for copying events 
  TTree *copy_tree1;
  TTree *copy_tree2;
  TTree *copy_tree3;
  TAna01Event *copy_evt1;
  TAna01Event *copy_evt2;
  TAna01Event *copy_evt3;
  TTree *copy_tree_discarded;
  TAna01Event *copy_evt_discarded;

  TChain     *fpChain;        // pointer to the analyzed TTree or TChain
  TFile      *fpHistFile;     // for output histograms and reduced trees
  TFile      *fWeightFile;
  TFile      *rho_file; 
  TFile      *mueff_file;   
  TFile      *functions;
   TFile      *fpCopyFile1;     // for skimmed trees
  TFile      *fpCopyFile2;  
  TFile      *fpCopyFile3;
  TFile      *fpDiscardedFile;   
  TString     fChainFileName; // the name of the chain file
  TString     fCutFile;       // contains file with the cut definitions
  int         fNentries;      // number of events in chain
  int         fEvent;         // current event number

  TAna01Event*fpEvt; 

  // -- Variables
  int        fRun; 
  int        fLS; 

  int NMuChk;
  int NScaled;
  int NFired;
  int Nnames; 
  int Nselected;  
  int Nsmeared;

  bool skim_analysis; 

  // counters for skimming

  int Nr1;
  int Nr2; 
  int Nr3; 
  int Nr5;
  int Nds; 


  // -- Histogram pointers 
  TTree *fTree;
  TTree *fTree1;
  TTree *fTree2;

  // -- Cut values
  double 
      PTLO
    , PTHI
    , ETALO
    , ETAHI   
    ;
  int TYPE;
  

};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
