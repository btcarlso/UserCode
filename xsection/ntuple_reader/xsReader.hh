#ifndef XSREADER_H
#define XSREADER_H

#include <iostream>

#include <map>
#include <memory>
#include <iostream>
#include <string>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TTree.h>

//#include "/afs/cern.ch/work/b/bcarlson/public/common_code/bins_new.card"

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaJet.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/JSON.hh"



#include "treeReaderXS.hh"


#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567




class xsReader: public treeReaderXS {
public:
  xsReader(TChain *tree, TString evtClassName);
  ~xsReader();

  void        bookHist();
  void        startAnalysis();
  void	      copy_event1();
  void	      copy_event2();
  void        copy_event3();
  void        copy_discarded();
  bool        MC_cuts(TLorentzVector P, TAnaTrack *pl1, TAnaTrack *pl2);
  int         numberOfPixLayers(TAnaTrack *pTrack) ;
  bool        cuts(TAnaCand *pCand, TAnaTrack *pTra1, TAnaTrack *pTra2, TAnaTrack *pl1, TAnaTrack *pl2, TAnaMuon *pMuon1, TAnaMuon *pMuon2);
  int         get_path();
  bool        check_label(int path);
  int	      check_pt();
  bool	      other_cuts();
  void        genMC();
  int         SingleMu_Match(int pTmin);
  void        listHLT(); 
  void        Trigger_Study();
  void		  MCefficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2);
  void        prob_efficiency(TLorentzVector P);
  void		  seagull_efficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2);
  void        eta_efficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2);
  void        impactPar_efficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2);
  bool        PS_cuts(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2,bool eta_cut);
  bool        getVectors();
  void        data_SG_eff();
  void	      select_dimuons(bool &PS, bool &muon_qual1, bool &muon_qual2, double &vprob);
  void        print_event_details();
  void        plot_kinematics();
  void        plot_errors();
  void		  plot_run_numbers(int RN);
  void        Pt_uncertainty();
  double      M_uncertainty();	
  void        eventProcessing();
  void        fillCandHist();
  void        UpsGun_acceptance();
  bool        isPathFired_Match(TString path, TString label);
  bool        isPathFired(TString path);
  void        fillCandVectors();
  void        cutvariable_histograms(TAnaCand *pCand, TAnaTrack *pTrack1, TAnaTrack *pTrack2, TAnaMuon *pMuon1, TAnaMuon *pMuon2);
  void        cutvariable_histograms_MC(TAnaTrack *pTrack1, TAnaTrack *pTrack2);	
  double      correction_functionJPsi(const double & pt, const double & eta, const double & phi, const int chg, double *par); // correction function JPsi
  double      correction_functionZ(const double & pt, const double & eta, const double & phi, const int chg,   double *parScale); //Z corr
  void        correct_Tracks(); 
  double      Pt_primeJPsi(TVector3 v1,int Q);
  double      Pt_primeZ(TVector3 v1,int Q); 	
  bool        CowboyVeto(TAnaCand *pCand);
  bool        CowboyVetoMC(TAnaTrack *pl1, TAnaTrack *pl2);
  int         eta_to_bin(double eta);	
  void        calculate_dimuon_efficiency();
  void        readCuts(TString filename, int dump = 1);
  void        candidateSelection(int mode);
  void        GetDataVectors(int best);
  void        freePointers();
  
  // -- Cut values
  int TYPE;
  int MODE;
  int MUTYPE1;
  int MUTYPE2;
  int RESTYPE;
  double PTHI;
  double PTMI;
  double PTLO;
  double ETAHI;
  double ETAMI;
  double ETABARREL;
  double PTBARREL;
  double PTCAND_MIN;
  double PTCAND;
  double RAPCAND;
  int BIN;
  double MASSLO;
  double MASSHI;
  int UPSTYPE;
  double DPHI;  
  double DETA;
  int BARREL;
  int PT_CORR;
  double JPSI_Z_CORR;
  double JPSI_Z_WEIGHT;
  double MU_ETA_MAX;
  double VPROB_MIN;	
  double VERTEX_SIG; 
  double DZ_MAX; 	
  double DXY_MAX;
  int    NV_HITS_MIN;
  double MUON_CHI2_DOF_MAX;

	
  TString HLTPATH, HLTPATH1, HLTPATH2, HLTPATH3, HLTPATH4, HLTPATH5, HLTPATH6, HLTPATH7, HLTPATH8, HLTPATH9, HLTPATH10, HLTPATH11, HLTPATH12, HLTPATH13;
  TString HLTPATH14, HLTPATH15;
  TString HLTLABEL, HLTLABEL1, HLTLABEL2, HLTLABEL3;
    
  // -- Variables
  TAnaCand    *fpCand; 
  double      fCandPt, fCandMass, fCandY;
  double      fGenMass;
  TLorentzVector fCand4V;
  TGenCand    *fgCand;
  double      fGenCandPt, fGenCandY;
  double      fGenMuon1Pt, fGenMuon1Eta, fMuon1Eta, fMuon1Pt;
  double      fGenMuon2Pt, fGenMuon2Eta, fMuon2Eta, fMuon2Pt;

  double fCandPtE; 	
  double fMuon1Phi, fMuon2Phi;
  double fMuonPhi12; 
  double      fMuon1PtE, fMuon2PtE, fMuon1EtaE, fMuon2EtaE, fMuon1PhiE, fMuon2PhiE;
  TVector3 fMuon1Vect, fMuon2Vect;
  TVector3 fMuon1gen, fMuon2gen;	
  double fGenPt; 
	
  PidTable    *fPidTableMuIDPos, *fPidTableMuIDNeg; 
  PidTable    *fPidTableTrigPos, *fPidTableTrigNeg;
  PidTable    *fPidTableTrckEff;
  PidTable    *fPidTableTrigFit, *fPidTableMuidFit;
  PidTable    *fPidTable2011SeagullPos, *fPidTable2011SeagullNeg; 
  PidTable    *fPidTable2011CowboyPos, *fPidTable2011CowboyNeg;  
  PidTable    *fPidTable1SLambdaThetaPos, *fPidTable1SLambdaThetaNeg, *fPidTable1SLambdaThetaPhiPos;
  PidTable    *fPidTable1SLambdaThetaPhiNeg, *fPidTable1SLambdaPhiNeg, *fPidTable1SLambdaPhiPos; 
  PidTable    *fPidTable2SLambdaThetaPos, *fPidTable2SLambdaThetaNeg, *fPidTable2SLambdaThetaPhiPos;
  PidTable    *fPidTable2SLambdaThetaPhiNeg, *fPidTable2SLambdaPhiNeg, *fPidTable2SLambdaPhiPos; 
  PidTable    *fPidTable3SLambdaThetaPos, *fPidTable3SLambdaThetaNeg, *fPidTable3SLambdaThetaPhiPos;
  PidTable    *fPidTable3SLambdaThetaPhiNeg, *fPidTable3SLambdaPhiNeg, *fPidTable3SLambdaPhiPos;   
  
 // TFile *rho_file;	
  //TFile *mueff_file; 	
	
  
  double      fWeight;
  double dimuon_efficiency; 
  double dimuon_efficiencyEp; 
  double dimuon_efficiencyEm; 
	/*
  static const int  fNpt = 24;
  static const int fNpt1=14;
  static const int fNpt2=3;
  static const int fNy2=6;
  static const int  fNy = 6;
  static const int fNmass=15;
  double      fPTbin[fNpt+1],fPTbin1[fNpt1+1], fYbin[fNy+1], fMassbin[fNmass+1];
  double fYbin2[fNy2+1], fPTbin2[fNpt2+1];
	 */
  int fBin;
  double fMassLow, fMassHigh;
  
  float um, uP, ue, up, gE, gP, ge, gp, dR, xbm;
  float m_um, m_uP, m_ue, m_up, m_gE, m_gP, m_ge, m_gp, m_dR, m_xbm;
  int m_xbid;
  float mbg_um, mbg_uP, mbg_ue, mbg_up, mbg_gE, mbg_gP, mbg_ge, mbg_gp, mbg_dR, mbg_xbm;
  
  vector<TAnaCand*> Cands;
  vector<TAnaCand*> Cands_ID;
  vector<TAnaCand*> Cands_TM;
  
  JSON   *fpJSON;
};


#endif
