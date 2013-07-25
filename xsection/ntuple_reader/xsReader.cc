

#include "xsReader.hh"
#include <vector>
#include "TRandom.h"
#include "TLorentzRotation.h"
#include "TDirectory.h"
#define MMUON 0.10566
#define MKAON 0.49368

double fYbin[]={0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6};
double fPTbin[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,25,30,50,70,100};
double fPTbin1[]={10.,12.0,14,16.0,18,20.0,22,24.0,26,28,30.0,32,34,36,38,40.0,43,46,50,55,60,70,100};
double fPTbin2[]={10.,12,14,16,18,20,22,24,26,28,30,35,50,70,100};
double fYbin2[]={0.0,0.6,1.2};
double fMassbin[]={8.7,8.9,8.95,9.2,9.25,9.435,9.485,9.71,9.76,9.76,10.0,10.05,10.33,10.38,10.9,10.95,11.2};

const int fNy=sizeof(fYbin)/sizeof(double)-1; 
const int fNpt=sizeof(fPTbin)/sizeof(double)-1;
const int fNpt1=sizeof(fPTbin1)/sizeof(double)-1;
const int fNpt2=sizeof(fPTbin2)/sizeof(double)-1; 
const int fNy2=sizeof(fYbin2)/sizeof(double)-1;
const int fNmass=sizeof(fMassbin)/sizeof(double)-1; 


// ----------------------------------------------------------------------
// Run with: ./runXSReaders -c chains/bg-test -D root 
//           ./runXSReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
xsReader::xsReader(TChain *tree, TString evtClassName): treeReaderXS(tree, evtClassName) {
  cout << "--> xsReader> This is the start ..." << endl;
  fpJSON = new JSON("Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt");//JSON used by polarization group
	
   ///// PidTable Tracking Efficiency for DATA
  fPidTableTrckEff = new PidTable("PidTables/DATA/Upsilon/PtTrackEff.dat");
  
  
  fPidTableTrigFit = new PidTable("PtTrigFit-jpsi.8ptbin.DATAv4.dat");
  fPidTableMuidFit = new PidTable("PtMuidFit-jpsi.8ptbin.DATAv4.dat");
  
  fPidTableTrigPos = new PidTable("PtMmbTrigBothv2-jpsi.8ptbin.DATA.dat");
  fPidTableMuIDPos = new PidTable("PtMmbMuidBothv2-jpsi.8ptbin.DATA.dat");
  
  fPidTable2011SeagullPos = new PidTable("Pt2011Seagull_PosErr.DATA.dat");
  fPidTable2011SeagullNeg = new PidTable("Pt2011Seagull_NegErr.DATA.dat");  
  
  fPidTable2011CowboyPos = new PidTable("Pt2011Cowboy_PosErr.DATA.dat");
  fPidTable2011CowboyNeg = new PidTable("Pt2011Cowboy_NegErr.DATA.dat");  
  
  fPidTable1SLambdaThetaPos = new PidTable("PidTables/Polarization/Ups1S/Pol_1S_LambdaTheta_PosErr.dat");
  fPidTable1SLambdaThetaNeg = new PidTable("PidTables/Polarization/Ups1S/Pol_1S_LambdaTheta_NegErr.dat");
  fPidTable1SLambdaPhiPos = new PidTable("PidTables/Polarization/Ups1S/Pol_1S_LambdaPhi_PosErr.dat");
  fPidTable1SLambdaPhiNeg = new PidTable("PidTables/Polarization/Ups1S/Pol_1S_LambdaPhi_NegErr.dat");  
  fPidTable1SLambdaThetaPhiPos = new PidTable("PidTables/Polarization/Ups1S/Pol_1S_LambdaThetaPhi_PosErr.dat");
  fPidTable1SLambdaThetaPhiNeg = new PidTable("PidTables/Polarization/Ups1S/Pol_1S_LambdaThetaPhi_NegErr.dat");
  
  fPidTable2SLambdaThetaPos = new PidTable("PidTables/Polarization/Ups2S/Pol_2S_LambdaTheta_PosErr.dat");
  fPidTable2SLambdaThetaNeg = new PidTable("PidTables/Polarization/Ups2S/Pol_2S_LambdaTheta_NegErr.dat");
  fPidTable2SLambdaPhiPos = new PidTable("PidTables/Polarization/Ups2S/Pol_2S_LambdaPhi_PosErr.dat");
  fPidTable2SLambdaPhiNeg = new PidTable("PidTables/Polarization/Ups2S/Pol_2S_LambdaPhi_NegErr.dat");  
  fPidTable2SLambdaThetaPhiPos = new PidTable("PidTables/Polarization/Ups2S/Pol_2S_LambdaThetaPhi_PosErr.dat");
  fPidTable2SLambdaThetaPhiNeg = new PidTable("PidTables/Polarization/Ups2S/Pol_2S_LambdaThetaPhi_NegErr.dat");  
  
  fPidTable3SLambdaThetaPos = new PidTable("PidTables/Polarization/Ups3S/Pol_3S_LambdaTheta_PosErr.dat");
  fPidTable3SLambdaThetaNeg = new PidTable("PidTables/Polarization/Ups3S/Pol_3S_LambdaTheta_NegErr.dat");
  fPidTable3SLambdaPhiPos = new PidTable("PidTables/Polarization/Ups3S/Pol_3S_LambdaPhi_PosErr.dat");
  fPidTable3SLambdaPhiNeg = new PidTable("PidTables/Polarization/Ups3S/Pol_3S_LambdaPhi_NegErr.dat");  
  fPidTable3SLambdaThetaPhiPos = new PidTable("PidTables/Polarization/Ups3S/Pol_3S_LambdaThetaPhi_PosErr.dat");
  fPidTable3SLambdaThetaPhiNeg = new PidTable("PidTables/Polarization/Ups3S/Pol_3S_LambdaThetaPhi_NegErr.dat");
   
//  rho_file = new TFile("rho_factor.root","READ"); 
 // mueff_file = new TFile("eff_all.root","READ"); 	
	
  
}
// ----------------------------------------------------------------------
xsReader::~xsReader() {
  cout << "--> xsReader> This is the end ..." << endl;
}

// ----------------------------------------------------------------------
void xsReader::startAnalysis() {
  cout << "--> xsReader> startAnalysis: ..." << endl;
}

void xsReader::copy_event1(){
	if(fpEvt){
		*copy_evt1 = *fpEvt;
		copy_tree1->Fill(); }
}

void xsReader::copy_event2(){
	if(fpEvt){
		*copy_evt2 = *fpEvt;
		copy_tree2->Fill(); }
}

void xsReader::copy_event3(){
	if(fpEvt){
		*copy_evt3 = *fpEvt;
		copy_tree3->Fill(); }
}

void xsReader::copy_discarded(){
	if(fpEvt){
		*copy_evt_discarded = *fpEvt;
		copy_tree_discarded->Fill();}
}

int xsReader::get_path(){
	int path(-9);
	
	if ( isPathFired(HLTPATH) )  path=0;
    if ( isPathFired(HLTPATH1) ) path=1;
    if ( isPathFired(HLTPATH2) ) path=2;
    if ( isPathFired(HLTPATH3) ) path=3;
    if ( isPathFired(HLTPATH4) ) path=4;
    if ( isPathFired(HLTPATH5) ) path=5;
    if ( isPathFired(HLTPATH6) ) path=6;  
    if ( isPathFired(HLTPATH7) ) path=7;
    if ( isPathFired(HLTPATH8) ) path=8;
    if ( isPathFired(HLTPATH9) ) path=9;
    if ( isPathFired(HLTPATH10) ) path=10; 
	
	//cout << "Path: " << path << endl; 
	
	return path; 
	
}

bool xsReader::check_label(int path){
	/*HLTPATH     0  HLT_Dimuon0_Barrel_Upsilon_v1 #HLTPATH goes with	HLTLABEL below
	 HLTLABEL    0  hltDimuon0BarrelUpsilonL3Filtered:HLT::
	 HLTPATH1     0  HLT_Dimuon5_Upsilon_Barrel_v1
	 HLTPATH2     0  HLT_Dimuon5_Upsilon_Barrel_v2
	 HLTPATH3     0  HLT_Dimuon5_Upsilon_Barrel_v3
	 HLTPATH4     0  HLT_Dimuon5_Upsilon_Barrel_v5
	 HLTLABEL1    0  hltBarrelUpsilonL3Filtered:HLT::
	 HLTPATH5     0  HLT_Dimuon7_Upsilon_Barrel_v1
	 HLTPATH6     0  HLT_Dimuon7_Upsilon_Barrel_v4
	 HLTLABEL2    0  hltBarrelDimuon7UpsilonL3Filtered:HLT::
	 HLTPATH7     0  HLT_Dimuon9_Upsilon_Barrel_v1
	 HLTPATH8     0  HLT_Dimuon9_Upsilon_Barrel_v4
	 HLTLABEL3    0  hltDimuon9BarrelUpsilonL3Filtered:HLT::
	 HLTPATH9     0  HLT_Dimuon7_Upsilon_Barrel_v5 #Goes with HLTLABEL2, dimuon7
	 HLTPATH10    0  HLT_Dimuon9_Upsilon_Barrel_v5
	 */
	bool tmp=1; 
	if(BARREL==0){
	//if (path ==0) if ( !isPathFired_Match(HLTPATH,HLTLABEL) ) tmp=0; 
    if (path ==1) if ( !isPathFired_Match(HLTPATH1,HLTLABEL) ) tmp=0;
    if (path ==2) if ( !isPathFired_Match(HLTPATH2,HLTLABEL1) ) tmp=0;
    if (path ==3) if ( !isPathFired_Match(HLTPATH3,HLTLABEL2) ) tmp=0;
    if (path ==4) if ( !isPathFired_Match(HLTPATH4,HLTLABEL2) ) tmp=0;
    if (path ==5) if ( !isPathFired_Match(HLTPATH5,HLTLABEL2) ) tmp=0;
    if (path ==6) if ( !isPathFired_Match(HLTPATH6,HLTLABEL2) ) tmp=0;
    if (path ==7) if ( !isPathFired_Match(HLTPATH7,HLTLABEL2) ) tmp=0;
    if (path ==8) if ( !isPathFired_Match(HLTPATH8,HLTLABEL2) ) tmp=0;
    if (path ==9) if ( !isPathFired_Match(HLTPATH9,HLTLABEL2) ) tmp=0;
	if (path ==10) if ( !isPathFired_Match(HLTPATH10,HLTLABEL2) ) tmp=0; }
	
	if(BARREL==1){
	
	//if (path ==0) if ( !isPathFired_Match(HLTPATH,HLTLABEL) ) tmp=0;
    if (path ==1) if ( !isPathFired_Match(HLTPATH1,HLTLABEL1) ) tmp=0;
    if (path ==2) if ( !isPathFired_Match(HLTPATH2,HLTLABEL1) ) tmp=0;
    if (path ==3) if ( !isPathFired_Match(HLTPATH3,HLTLABEL1) ) tmp=0;
    if (path ==4) if ( !isPathFired_Match(HLTPATH4,HLTLABEL1) ) tmp=0;
    if (path ==5) if ( !isPathFired_Match(HLTPATH5,HLTLABEL2) ) tmp=0;
    if (path ==6) if ( !isPathFired_Match(HLTPATH6,HLTLABEL2) ) tmp=0;
    if (path ==7) if ( !isPathFired_Match(HLTPATH7,HLTLABEL3) ) tmp=0;
    if (path ==8) if ( !isPathFired_Match(HLTPATH8,HLTLABEL3) ) tmp=0;
    if (path == 9)  if ( !isPathFired_Match(HLTPATH9,HLTLABEL2) ) tmp=0;
	if (path == 10) if ( !isPathFired_Match(HLTPATH10,HLTLABEL3) ) tmp=0;
	//HLTPATH9 goes with HLTLABEL2 since both correspond to dimuon7
	}
	
	//cout << "BARREL: " << BARREL << " PATH: " << path << " TMP: " << tmp << endl; 
	
	return tmp; 
}

int xsReader::check_pt(){
	//  cout << "Pt check: " << fpEvt->nCands() << endl; 
	TAnaCand *pCand(0);
	bool c1=0; 
	bool c2=0; 
	bool c3=0; 
	bool c4=0; 
	bool cpt1=0; 
	bool cpt2=0; 
	
	TAnaTrack *pl1(0); TAnaTrack *pl2(0);
	
	for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
		pCand = fpEvt->getCand(iC);

		double ptCand=pCand->fPlab.Perp();

		if(ptCand<9)
		  return 5;

		//cout <<"test" <<endl;  
		c1=0; 
		c2=0; 
		c3=0;
		c4=0; 
		cpt1=0; 
		cpt2=0; 
		pl1 = fpEvt->getSigTrack(pCand->fSig1);
		pl2 = fpEvt->getSigTrack(pCand->fSig2);
		
		double pt1=pl1->fPlab.Perp(); 
		double pt2=pl1->fPlab.Perp();
		double eta1=TMath::Abs(pl1->fPlab.Eta()); 
		double eta2=TMath::Abs(pl2->fPlab.Eta()); 
		// cout <<"pt1: " << pt1 << " pt2: "<< pt2<<endl; 
		//cout << "eta1: " << eta1 << " eta2 " << eta2 <<endl; 
		//PT cuts. If Pt outside of the given range, then go to the end of the function 
		
		if(pt1<3)
			cpt1=true;
		if(pt2<3)
			cpt2=true;
		
		if ( eta1 <1.2 ){
			if(pt1 < 3.5 )
				c1=true; }
		if ( eta2 <1.2 ){
			if(pt2 < 3.5 )
				c3=true; }
		
		if(eta2>1.2 && eta2<2.4)
			if(pt2<3)
				c4=true;
		if(eta1>1.2 && eta1<2.4)
			if(pt1<3)
				c2=true;
		
		if(cpt1==1 || cpt2==1 || c1==1 ||c2==1 || c3==1 || c4==1 ){
			/*
			 cout << "Discard:" <<endl;
			 cout << "c1: " << c1 << " c2: " << c2 << " c3: " << c3 << " c4: " << c4 << endl;
			 cout << "eta1: " << eta1 <<" pt1 " << pt1 <<endl;
			 cout << "eta2: "<< eta2 << " pt2 " << pt2 <<endl;
			 */ 
			return 0; }
		
		if(eta1<1.2 && eta2<1.2)
			if(pt1>3.5 && pt2>3.5)
				return 1;
		
		if(eta1>1.2 && eta2>1.2 && eta1<2.4 && eta2<2.4)
			if(pt1>3 && pt2>3)
				return 2;
		/*
		 cout << endl; 
		 cout << "eta1: " << eta1 <<" pt1 " << pt1 <<endl;
		 cout << "eta2: "<< eta2 << " pt2 " << pt2 <<endl;
		 */
	}//end for loop  
	
	return 3; 
	
}

bool xsReader::other_cuts(){
	bool good_event =0;  
	TAnaCand *pCand(0);
	
	TAnaTrack *pl1(0); TAnaTrack *pl2(0);
	
	for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
		good_event=0; 
		pCand = fpEvt->getCand(iC);
		
		pl1 = fpEvt->getSigTrack(pCand->fSig1);
		pl2 = fpEvt->getSigTrack(pCand->fSig2);
		
		TAnaTrack *pTrack1 = fpEvt->getRecTrack(pl1->fIndex);
		TAnaTrack *pTrack2 = fpEvt->getRecTrack(pl2->fIndex);
		
		int index1=0; 
		int index2=0; 
		
		double mass=pCand->fMass; //mass for each candidate
		double ptCand=pCand->fPlab.Perp();
		double candeta=pCand->fPlab.Eta();
		
		// Find the correct Muons. 
		for(int i=0; i<fpEvt->nMuons();i++){
			TAnaMuon *testMuon = fpEvt->getMuon(i);
			if(testMuon->fIndex==pl1->fIndex)
				index1=i; 
			if(testMuon->fIndex==pl2->fIndex)
				index2=i;
		}//end muon loop
		
		TAnaMuon *pMuon1 = fpEvt->getMuon(index1);
		TAnaMuon *pMuon2 = fpEvt->getMuon(index2);
		
		if (TYPE != pCand->fType) continue;
		
		if(pCand->fVtx.fProb<VPROB_MIN) continue; 
		
		if(VERTEX_SIG>0){
			if(pCand->fVtx.fD3d/pCand->fVtx.fD3dE > VERTEX_SIG) continue; 
		}
		
		if(TMath::Abs(pTrack1->fdz)>DZ_MAX || TMath::Abs(pTrack2->fdz)>DZ_MAX) continue; 
		
		if(TMath::Abs(pTrack1->fdxy)>DXY_MAX ||TMath::Abs(pTrack2->fdxy)>DXY_MAX) continue; 
		
		if(pTrack1->fValidHits<NV_HITS_MIN || pTrack2->fValidHits<NV_HITS_MIN) continue; 
		
		// if((pMuon1->fMuonChi2)/static_cast<double>(pMuon1->fTimeNdof)>1.8 || pMuon2->fMuonChi2/static_cast<double>(pMuon2->fTimeNdof)>1.8) continue; 
		
		if( (((pMuon1->fMuID)&(1<<4))>>4)!=1 || (((pMuon1->fMuID)&(1<<12))>>12)!=1) continue; 
		
		if( (((pMuon2->fMuID)&(1<<4))>>4)!=1 || (((pMuon2->fMuID)&(1<<12))>>12)!=1) continue; 
		
		if ( (pl1->fPlab.Perp() > PTHI) || (pl1->fPlab.Perp() < PTLO) ) continue;
		
		if(ETAHI>0){
			if ( TMath::Abs(pl1->fPlab.Eta()) > ETAHI || TMath::Abs(pl2->fPlab.Eta())>ETAHI)  {
								continue;}
		}
		
		///////Muon in the barrel 
		if ( ((TMath::Abs(pl1->fPlab.Eta()) <= ETABARREL) && (pl1->fPlab.Perp() < PTBARREL)) || ((TMath::Abs(pl2->fPlab.Eta()) <= ETABARREL) && (pl2->fPlab.Perp() < PTBARREL)) ){
			continue;
		}
		((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon |#eta|<%.1f, P_{T}<%.1fGeV",ETABARREL,PTBARREL), 1); 
		
		if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETABARREL) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAMI) && (pl1->fPlab.Perp() < PTMI))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETABARREL) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAMI) && (pl2->fPlab.Perp() < PTMI)) ){
			continue;
		}
		((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon 1 or 2 %.1f<|#eta|<%.1f, P_{T}>%.1fGeV",ETABARREL,ETAMI,PTMI), 1); 
		
		
		if ( pl1->fQ*pl2->fQ > 0 ) continue;
		
		if (pCand->fMass < MASSLO) continue;
		if (pCand->fMass > MASSHI) continue;
		good_event=1; 
	}
	
	return good_event; 
	
}

bool xsReader::MC_cuts(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2){
  	bool print=0; 
	
	//((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Cowboy Veto",1);
	if(pl1==0 || pl2 ==0){
    
    cout << "Bad track. MC_cuts." << endl;
    return 0; 
  }
	if(print)cout  << "xsReader::MC_cuts " << endl; 
    cutvariable_histograms_MC(pl1,pl2); 
  /*
  cout << "INFO FOR FIRST & SECOND TRACK." << endl; 
  cout << "Cand: ";
  P.Print();
  cout << "Track1: ";
  (pl1->fPlab).Print();
  cout << "Track2: ";
  (pl2->fPlab).Print();
  */

  double mass=P.M();
  double ptCand = P.Pt();
  double candeta = P.Eta();	
  double candphi = P.Phi();
	
  TAnaTrack *pTra1=pl1; 
  TAnaTrack *pTra2=pl2;
	
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Total Events",1);	  
	
	//Various muon cuts 
	
	if ( (pl1->fPlab.Perp() > PTHI) || (pl1->fPlab.Perp() < PTLO) ) return 0;
	if ( (pl2->fPlab.Perp() > PTHI) || (pl2->fPlab.Perp() < PTLO) ) return 0;
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon %.1f<P_{T}<%.1fGeV",PTLO,PTHI),1);
	
	if(ETAHI>0){
	if ( TMath::Abs(pl1->fPlab.Eta()) > ETAHI || TMath::Abs(pl2->fPlab.Eta())>ETAHI)  {
		((TH1D*)fpHistFile->FindObjectAny("meta"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("meta2"))->Fill(mass,ptCand);
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon |#eta|<%.1f",ETAHI),1); 
	}
		
	///////Muon in the barrel 
	if ( ((TMath::Abs(pl1->fPlab.Eta()) <= ETABARREL) && (pl1->fPlab.Perp() < PTBARREL)) || ((TMath::Abs(pl2->fPlab.Eta()) <= ETABARREL) && (pl2->fPlab.Perp() < PTBARREL)) ){
		((TH1D*)fpHistFile->FindObjectAny("mPt_region1"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mPt_region1_2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPt_region1"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPt_region1"))->Fill(pl2->fPlab.Perp());
		return 0; 
	}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon |#eta|<%.1f, P_{T}<%.1fGeV",ETABARREL,PTBARREL), 1); 
	
	if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETABARREL) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAMI) && (pl1->fPlab.Perp() < PTMI))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETABARREL) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAMI) && (pl2->fPlab.Perp() < PTMI)) ){
		((TH1D*)fpHistFile->FindObjectAny("mPt_region2"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mPt_region2_2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl2->fPlab.Perp());
		return 0; 
	}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon 1 or 2 %.1f<|#eta|<%.1f, P_{T}>%.1fGeV",ETABARREL,ETAMI,PTMI), 1); 
	
	
	
	if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETAMI) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAHI) && (pl1->fPlab.Perp() < PTLO))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETAMI) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAHI) && (pl2->fPlab.Perp() < PTLO)) ){
		((TH1D*)fpHistFile->FindObjectAny("mPt_region2"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mPt_region2_2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl2->fPlab.Perp());
		
		// cout << "Eta> " << ETAMI << "and Eta< " << ETAHI << " and PT< " << PTLO<< endl; 
		// cout << "pt1: " << pl1->fPlab.Perp() << " eta1: " << pl1->fPlab.Eta() << endl;
		// cout << "pt2: " << pl2->fPlab.Perp() << " eta2: " << pl2->fPlab.Eta() << endl;
		
		return 0; 
	}
	
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon 1 or 2 %.1f<|#eta|<%.1f, P_{T}>%.1fGeV",ETAMI,ETAHI,PTLO), 1); 
	
	
	//Various candidate kinematic cuts 
	if ( P.Rapidity() < -RAPCAND ){ 
		return 0;}
	
	if ( P.Rapidity() > RAPCAND ){
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Cand |y|<%.1f",RAPCAND),1);
	
	if ( P.Perp() < PTCAND_MIN ){ 
		return 0;}
	
	if ( P.Perp() > PTCAND ){
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("%.1f<CandP_{T}<%.1fGeV",PTCAND_MIN,PTCAND),1);
	
	double pt1=pl1->fPlab.Perp();
	double pt2=pl2->fPlab.Perp();	  
	
	/// Vertex Cuts
	TAnaCand *pCand(0);
	TGenCand *pGenCand(0);
	
	for(int iC=0; iC<fpEvt->nGenCands();++iC){
		pGenCand=fpEvt->getGenCand(iC);
		if(pGenCand->fID == RESTYPE){
			((TH2D*)fpHistFile->FindObjectAny("GenVPhidxy"))->Fill(pGenCand->fV.Phi(),pTra1->fdxy);
			((TH2D*)fpHistFile->FindObjectAny("GenVPhidxy"))->Fill(pGenCand->fV.Phi(),pTra2->fdxy);
		}
	}

	if(fpEvt->nCands()==0) return 0;
	pCand=fpEvt->getCand(0);
	
	
	((TH2D*)fpHistFile->FindObjectAny("CandVPhidxy"))->Fill((pCand->fVtx.fPoint.Phi()),pTra1->fdxy);
	((TH2D*)fpHistFile->FindObjectAny("CandVPhidxy"))->Fill((pCand->fVtx.fPoint.Phi()),pTra2->fdxy);

	/// Vertex Cuts
	//fProb0.01	
	if(pCand->fVtx.fProb<VPROB_MIN){
		((TH1D*)fpHistFile->FindObjectAny("mProb"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mProb2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPtProb"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPtProb"))->Fill(pl2->fPlab.Perp());
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Vert. Prob.>%.2f", VPROB_MIN),1);
    
	// Remove this cut for optimization
	//>vertexsig>3	  
	  
	
    //////////
    //fdz>15
	if( TMath::Abs(pTra1->fdz) > DZ_MAX || TMath::Abs(pTra2->fdz) > DZ_MAX){ 
		((TH1D*)fpHistFile->FindObjectAny("mdz"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mdz2"))->Fill(mass,ptCand);
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("|dz|<%.1fcm",DZ_MAX),1); 	  
	
	//dxy>0.3  
	if( TMath::Abs(pTra1->fdxy)> DXY_MAX ||TMath::Abs(pTra2->fdxy) > DXY_MAX){
		
		((TH1D*)fpHistFile->FindObjectAny("mdxy"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mdxy2"))->Fill(mass,ptCand);
		
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("|dxy|>%.1fcm",DXY_MAX*10),1);	  
	//Number of hits <10	  
	if( pTra1->fValidHits < NV_HITS_MIN || pTra2->fValidHits < NV_HITS_MIN){ 
		((TH1D*)fpHistFile->FindObjectAny("mHits"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mHits2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPtHits"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPtHits"))->Fill(pl2->fPlab.Perp());
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Valid Hits>%d",NV_HITS_MIN),1);

	if ( pl1->fQ*pl2->fQ > 0 ){
		((TH1D*)fpHistFile->FindObjectAny("mfQ"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mfQ"))->Fill(mass,ptCand);
		return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Muons have opposite chage",1);
	
    if (mass < MASSLO) return 0;
    if (mass > MASSHI) return 0;
    ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Mass Range: %.1f<M_{#mu#mu}<%.1fGeV",MASSLO,MASSHI),1);
	
	if(print) cout << "eff Num." << endl; 
	((TH1D*)fpHistFile->FindObjectAny("Eff_vertexsig_den"))->Fill(ptCand);
	
	
	//compute vertex efficiency (not in pt bins) by dividing this histogram by the number of entries above a given value. This will require a bit more code 
	
	
	if(VERTEX_SIG>0){
	if(pCand->fVtx.fD3d/pCand->fVtx.fD3dE > VERTEX_SIG){
		((TH1D*)fpHistFile->FindObjectAny("mVertexSig"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mVertexSig2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPtVertexSig"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPtVertexSig"))->Fill(pl2->fPlab.Perp());
		return 0;}
	}
	if(print) cout << "Eff Deno. " << endl; 
	((TH1D*)fpHistFile->FindObjectAny("Eff_vertexsig_num"))->Fill(ptCand);// for vertex cut efficiency calculation
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("|fD3d/fD3dE|<%.1f",VERTEX_SIG),1);	

	for (int iy=0; iy<fNy2; iy++) {
		if(TMath::Abs(P.Rapidity())>fYbin2[iy] && TMath::Abs(P.Rapidity())<fYbin2[iy+1])
			((TH2D*)fpHistFile->FindObjectAny(Form("eta1_eta2_aftercuts_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(TMath::Abs(pl1->fPlab.Eta() ),TMath::Abs(pl2->fPlab.Eta()) ); 
	}
	
	
  return 1; 
   
}

void xsReader::eventProcessing() {
  // DiTracks();
  
  int path(-9);
  
  if(MODE == 5){
   //MC analysis 
	if(PT_CORR==1) correct_Tracks();
	UpsGun_acceptance();
	genMC();   

    if(getVectors()==0){
      // cout << "Event does not meet cuts." << endl; 
      goto end;
    }
 
    plot_kinematics();
    fillCandHist();
	M_uncertainty();
	Pt_uncertainty();
	plot_errors();

	  
    goto end; 
  }

  
	//Default running mode for data starts here 
	
    if(PT_CORR==1) correct_Tracks();//apply momentum corrections 
	
	fillCandVectors(); // fills Cands
  
//Skimming Code: 
	
	if(skim_analysis==0){
		
		if ( BARREL == 1 ){
			//Upsilons with the Barrel trigger tag
			path = get_path();
			if (path <0) goto end;
			
			if(check_label(path)==0)
				goto end; 
			else 
				NFired++;
		}
		if(fRun<165088) cout << "Before Run2011A promtRecov4" << endl; 	
		candidateSelection(4); 
		
		if(fpCand!=0) 
			copy_event1();
		
		/*
		path=get_path();
		bool names_good;
		if(path<0)
			names_good=0; 
		 else 
			 names_good=check_label(path); 
		int pt_good=0;
		
		if(names_good){
			Nnames++;
			pt_good=check_pt();
			//!!! USE ONLY TO IMPOSE ALL CUTS. Comment next line out for most general skim.  
			names_good=other_cuts();
		}
		
		if(names_good && pt_good==5){
		  //Pt<10GeV discard event. 
		  Nr5++;
		}


		if(names_good && pt_good==1){
		copy_event1();
			Nr1++;
		}
		if(names_good && pt_good==2){
			copy_event2();
			Nr2++;
		}
		
		if(names_good && pt_good==3){
			copy_event3();      	
			Nr3++;
		}
		if(names_good && pt_good==0){
			copy_discarded();
			Nds++;
		}
		*/ 
	} //end skimming
else { 
 
  if(!fpJSON->good(fRun, fLS)) goto end; 	// Check run & LS with JSON file 
	// Cowboys are veto'd at the trigger level after run 170249 
//	if(fRun>170249) goto end; 	
// listHLT(); //make a histogram of all the triggers. 	
//	Trigger_Study();
// to do the rho factor study, all you have to do is uncomment the listHLT function. Then comment out everything else. 	
	
	
  // analysis section for data
  if ( BARREL == 0 ){
	  //Analysis for Upsilons without barrel trigger
	path=get_path();
	if(path<0) goto end; // Added by BTC. Not in Bora's code. Why? 
	
	if(check_label(path)==0)
		goto end; 
	  
  }//BARREL==0
	
  if ( BARREL == 1 ){
	  //Upsilons with the Barrel trigger tag
	path = get_path();
	if (path <0) goto end;
	
	if(check_label(path)==0)
		goto end; 
	else 
		NFired++;
  }
	 
 
  if(fRun<165088) cout << "Before Run2011A promtRecov4" << endl; 	
  candidateSelection(4); // select candidates - including application of the cuts 
	
  if ( 0 != fpCand  ){
	  /*
	  data_SG_eff();
	  dimuon_efficiency=1; 
	  dimuon_efficiencyEp=1; 
	  dimuon_efficiencyEm=1; 
	  */
	  
    calculate_dimuon_efficiency(); 
    fillCandHist(); 
 
   // ((TH1D*)fpHistFile->FindObjectAny("Candidates"))->Fill(Cands_TM.size());
    ((TH1D*)fpHistFile->FindObjectAny("Paths"))->Fill(path);
	  plot_run_numbers(fpEvt->fRunNumber);// plot the run number ranges 
	  ((TH1D*)fpHistFile->FindObjectAny("nPV"))->Fill(fpEvt->nPV()); 
	  M_uncertainty();
	  Pt_uncertainty();
	  plot_kinematics();
	  plot_errors();
	  Nselected++;
	 // cout << " pt: " << fCandPt << endl; 
	 // cout << " Y: " << fCandY << endl; 
	  //Mode 4 used to be Ups_iso
    
  }//there is a valid candidate
  

  fpHistFile->cd();
}//end of analysis portion of code 
  end:
  freePointers();  
  
}

void xsReader::plot_run_numbers(int RN){

	string data_set="INIT";
	if(RN<165088)
		data_set="Not-in-data-set";
	if(RN>=165088 && RN<=167913)
		data_set="Run2011A-PromtReco-v4";
	if(RN>=170722 && RN<=172619)
		data_set="Run2011A-PromptReco-v5";
	if(RN>=172620 && RN<=173692)
		data_set="Run2011A-PromptReco-v6";
	if(RN>175832 && RN<=180252)
		data_set="Run2011B-PromptReco-v1"; 
	
	((TH1I*)fpHistFile->FindObjectAny("DataSets"))->Fill(data_set.c_str(),1);  
	
}


void xsReader::genMC(){
  /*
    This function is for looking at MC. It first loops over all the Gen-level candidates, then looks for candidates that correspond to 
    Y1S,Y2S,Y3S. The type is set in the cuts file, as UPSTYPE. Once a cadidate is found, the momentum 4 vector is set. The program them loops
    over all the daughters, and gets the candidate for each daughter. We look at only particles with fID of 13, which correspond to muons. +13 to
    muon- and -13 to muon+. We want to loop over all the tracks for each muon candidate. If the track genIndex matches the candidate number, we 
    consider that a good track, and write the momentum into a 4-vector. 
  */

 
  bool print=0;
  if(print==1)  cout << "GEN INFO. " << endl;  
  TGenCand *gCand(0);//Candidate

  TGenCand *gDau1(0);//Daughter 1
  TGenCand *gDau2(0);//Daugther 2

  TLorentzVector genCand;//Upsilon gen level candidate
  TLorentzVector pMuon1reco; // muon tracks
  TLorentzVector pMuon2reco;
  TLorentzVector pMuon1gen;
  TLorentzVector pMuon2gen;
  TLorentzVector pGamma;
  TLorentzVector pUpsgen;

  TGenCand *g2Cand(0);//daughter cand

  TAnaTrack *pTrack(0);//variable for muon tracks

  double pt(-9);//temporary variables to get track momentum components
  double eta(-9);
  double phi(-9);

  int par_id=-9;
  int dau_id=13; //13 for -muon, and -13 for +muon.
  int track_counter=1; 
  bool photon_check=0; 
  int genmuon=1; 

  //Define parent ID's for upsilon 1S,2S,3S

  
  //Temp check of vertex understanding

  TGenCand *pGenCand(0);
  TAnaCand *pCand(0);

  for(int iC=0; iC<fpEvt->nGenCands();++iC){
    pGenCand=fpEvt->getGenCand(iC);
    if(pGenCand->fID == 553) {
      if(print) cout << "Upsilon Gen Position: ";
      if(print) (pGenCand->fV).Print();
      ((TH2D*)fpHistFile->FindObjectAny("GenVxy"))->Fill(pGenCand->fV.X(),pGenCand->fV.Y());
      ((TH1D*)fpHistFile->FindObjectAny("GenVz"))->Fill(pGenCand->fV.Z());
    }
    if(TMath::Abs(pGenCand->fID) == 13){
      if(print) cout << "Muon Gen Position: ";
      if(print) (pGenCand->fV).Print();
    }
    }

  for(int iC=0; iC<fpEvt->nCands();++iC){
    pCand=fpEvt->getCand(iC);
    if(print) cout << "Cand Position: ";
    if(print) (pCand->fVtx.fPoint).Print();
    if(print) cout <<"pcand->fVtx.fD3d: " << pCand->fVtx.fD3d << endl;
    if(print) cout <<"pCand->fVtx.fD3dE: " << pCand->fVtx.fD3dE<<endl; 
    if(print) cout <<"pCand->fVtx.fD3d/pCand->fVtx.fD3dE=" << pCand->fVtx.fD3d/pCand->fVtx.fD3dE<<endl; 
    ((TH2D*)fpHistFile->FindObjectAny("CandVxy"))->Fill(pCand->fVtx.fPoint.X(),pCand->fVtx.fPoint.Y());
    ((TH1D*)fpHistFile->FindObjectAny("CandVz"))->Fill(pCand->fVtx.fPoint.Z());

  }

  if(print)  cout << "-------------" << endl; 

  TGenCand *test(0);
  if(fpEvt->nRecTracks()!=2){
  for(int iG=0; iG<fpEvt->nGenCands();++iG){
    test = fpEvt->getGenCand(iG);
    if(print) cout << "iG: " << iG << " nRecTracks: " << fpEvt->nRecTracks() << " Cand ID: " << test->fID << endl; 
     }
  }
  TAnaTrack *track_temp(0);

  if(fpEvt->nRecTracks()!=2){
    for(int iR=0; iR<fpEvt->nRecTracks();++iR){
      track_temp=fpEvt->getRecTrack(iR);
      if(print) cout << "iR: " << iR << " MCID: " << track_temp->fMCID << endl; 
    }

  }

  if(UPSTYPE==1)
    par_id=553;//Y1S
  if(UPSTYPE==2)
    par_id=100553;//Y2S
  if(UPSTYPE==3)
    par_id=200553;//Y3S

  for(int iG=0; iG<fpEvt->nGenCands();++iG){
    //loop over GenCands. 
    gCand=fpEvt->getGenCand(iG);
    if(gCand->fID == par_id){
      genCand.SetPtEtaPhiE(gCand->fP.Perp(),gCand->fP.Eta(),gCand->fP.Phi(),gCand->fP.Energy());
      for(int i=gCand->fDau1; i<=gCand->fDau2;++i){
	g2Cand=fpEvt->getGenCand(i);
	if(13==TMath::Abs(g2Cand->fID)){
	  for(int iR=0; iR<fpEvt->nRecTracks();iR++){
	    pTrack=fpEvt->getRecTrack(iR);
	    if(print)    cout << "PvIdx: " << pTrack->fPvIdx << endl; 
	    if(print)    cout << "Track " <<iR <<" dxy: " << pTrack->fdxy << endl; 
	    if(print)    cout << "fTip: " << pTrack->fTip << " fBsTip: " << pTrack->fBsTip << endl; 
	    if(print)    cout << "fLip: " << pTrack->fLip << " fBsLip: " << pTrack->fBsLip << endl; 
	    if(pTrack->fGenIndex == g2Cand->fNumber ){
	      pt=pTrack->fPlab.Pt();
	      eta=pTrack->fPlab.Eta();
	      phi=pTrack->fPlab.Phi();
	      if(track_counter==1) {
		pMuon1reco.SetPtEtaPhiM(pt,eta,phi,MMUON);//set muon track 4 vector
		pt=-9;
		eta=-9;
		phi=-9;
		track_counter++;
	      }//first  muon set momentum 
	      if(track_counter==2){
		pMuon2reco.SetPtEtaPhiM(pt,eta,phi,MMUON);
		pt=-9;
		eta=-9;
		phi=-9;
	      }//second muon set momentum
	    }//track check
	  }//loop over tracks

	}//daugher is muon, +/-
      }//Loop over daughters
    }//if cand is Upsilon     
    if(gCand->fID == 22){
      pGamma.SetPtEtaPhiE(gCand->fP.Perp(),gCand->fP.Eta(),gCand->fP.Phi(),gCand->fP.Energy());
      photon_check=1; 
    }
    if(TMath::Abs(gCand->fID) == 13){
	 if(genmuon==1) {
	   pMuon1gen.SetPtEtaPhiE(gCand->fP.Perp(),gCand->fP.Eta(),gCand->fP.Phi(),gCand->fP.Energy());
	   genmuon++;
	 }
	 if(genmuon==2){
	   pMuon2gen.SetPtEtaPhiE(gCand->fP.Perp(),gCand->fP.Eta(),gCand->fP.Phi(),gCand->fP.Energy());
	 }
       }
  
  }//loop of cands
  pUpsgen=pMuon1gen+pMuon2gen;
  if(genmuon<2){
    cout << "Not 2 muons: " << pUpsgen.M() << endl;        
      }
  if(print){
  if(photon_check!=1){
    cout << "No photon: "<< pUpsgen.M() << endl;}
  if(photon_check==1){
    cout << "Photon: " <<  pUpsgen.M() <<endl;
  }

  cout << "Number of Reco tracks: " << fpEvt->nRecTracks() << endl; 
  cout << "Momentum: " << endl; 
  cout <<"Cand: ";
  genCand.Print();

  cout <<"Reco Track1: ";
  pMuon1reco.Print();
  cout << "Gen Track1: "; 
  pMuon1gen.Print();

  cout << "Reco Track2: ";
  pMuon2reco.Print();
  cout << "Gen Track2: "; 
  pMuon2gen.Print();

  if(photon_check){
     cout << "Photon: ";
     pGamma.Print();
  }
    cout << endl; 
    if(track_counter!=2) cout << "Issue because 2 tracks were not found." << endl; 
  }//print

  fMuon1gen.SetXYZ(pMuon1gen.X(),pMuon1gen.Y(),pMuon1gen.Z()); 
  fMuon2gen.SetXYZ(pMuon2gen.X(),pMuon2gen.Y(),pMuon2gen.Z());
  fGenPt=pUpsgen.Perp();
	
  ((TH1D*)fpHistFile->FindObjectAny("genUpsilonMass"))->Fill(pUpsgen.M());
  ((TH1D*)fpHistFile->FindObjectAny("gen_reco_ptdiff"))->Fill((pUpsgen.Perp()-genCand.Perp())*1000); 
  if(pGamma.E()>0.0005)((TH1D*)fpHistFile->FindObjectAny("Egamma"))->Fill(pGamma.E());
  fGenMass=pUpsgen.M();

  //  if(fGenMass>9.4 && fGenMass<9.51)cout << "m: " << fGenMass << endl; 

}

void xsReader::data_SG_eff(){
	
	
	if(fCandY<0.6){
		if(fMuon1Vect.DeltaPhi(fMuon2Vect)<=0){
			((TH1D*)fpHistFile->FindObjectAny("seagull_numerator"))->Fill(fCandPt);
			if(fCandMass>9.26 && fCandMass<=9.66) ((TH1D*)fpHistFile->FindObjectAny("seagull_numerator_SB"))->Fill(fCandPt);
		}
		if(fCandMass>9.26 && fCandMass<=9.66) ((TH1D*)fpHistFile->FindObjectAny("seagull_denominator_SB"))->Fill(fCandPt);
		((TH1D*)fpHistFile->FindObjectAny("seagull_denominator"))->Fill(fCandPt); 
	}
	
	
}

void xsReader::MCefficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2){
	if(fpEvt->nCands()==0) return;
	eta_efficiency(P,pl1,pl2);
	bool PS=PS_cuts(P,pl1,pl2,0);
	if(!PS) return; 
	
	seagull_efficiency(P, pl1,pl2);
	prob_efficiency(P);
	impactPar_efficiency(P, pl1,pl2); 
	
}

void xsReader::prob_efficiency(TLorentzVector P){
	//Compute the efficiency of the vertex probability cut in rapidity & Pt bins 

	TAnaCand *pCand(0);
	
	pCand=fpEvt->getCand(0);
	double Y=TMath::Abs(P.Rapidity());

	for (int iy=0; iy<fNy2; iy++) {
		if(Y>fYbin2[iy] && Y<fYbin2[iy+1]){
			((TH2D*)fpHistFile->FindObjectAny(Form("vertex_prob_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp(),pCand->fVtx.fProb);
		}//iy
	}	//for
		
	
}

void xsReader::eta_efficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2){
	
	double Y=TMath::Abs(P.Rapidity());
	double Pt=P.Perp();
	double eta1=TMath::Abs(pl1->fPlab.Eta());
	double eta2=TMath::Abs(pl2->fPlab.Eta());

	bool PS=PS_cuts(P,pl1,pl2,1);
	
	if(!PS) return; 
	
	for (int iy=0; iy<fNy2; iy++) {
		if(Y>fYbin2[iy] && Y<fYbin2[iy+1])
			((TH2D*)fpHistFile->FindObjectAny(Form("eta1_eta2_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(eta1,eta2); 
	}
	
	for (int iy=0; iy<fNy2; iy++) {
		for (int ipt=0; ipt<fNpt1; ipt++) {
			if(Y>fYbin2[iy] && Y<fYbin2[iy+1] && Pt>fPTbin1[ipt] && Pt<fPTbin1[ipt+1])
				((TH2D*)fpHistFile->FindObjectAny(Form("eta1_eta2_bin1_y%.1f_%.1f_pt%.1f_%.1f",fYbin2[iy],fYbin2[iy+1],fPTbin1[ipt],fPTbin1[ipt+1])))->Fill(eta1,eta2);
		}
		
	}//end of iy loop 
}

void xsReader::impactPar_efficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2){

	double Y=TMath::Abs(P.Rapidity());
	
	for (int iy=0; iy<fNy2; iy++) {
		if(Y>fYbin2[iy] && Y<fYbin2[iy+1]){
			((TH2D*)fpHistFile->FindObjectAny(Form("dxy_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp(),TMath::Abs(pl1->fdxy) );
			((TH2D*)fpHistFile->FindObjectAny(Form("dz_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp(),TMath::Abs(pl1->fdz) );
			((TH2D*)fpHistFile->FindObjectAny(Form("dxy_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp(),TMath::Abs(pl2->fdxy) );
			((TH2D*)fpHistFile->FindObjectAny(Form("dz_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp(),TMath::Abs(pl2->fdz) );

		}//iy
	}	//for
	
}


void xsReader::seagull_efficiency(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2){
	bool PS=PS_cuts(P,pl1,pl2,0);
	bool seagull=!CowboyVetoMC(pl1,pl2);
		
	double Y=TMath::Abs(fCand4V.Rapidity());
	
	for(int iy=0; iy<fNy2; iy++){
		if(Y>fYbin2[iy] && Y<fYbin2[iy+1]){
			if(PS && seagull)((TH1D*)fpHistFile->FindObjectAny(Form("seagull_num_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp()); 
			if(PS) ((TH1D*)fpHistFile->FindObjectAny(Form("seagull_den_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1])))->Fill(P.Perp()); 
			
		}//if
	}//for
	
}

bool xsReader::PS_cuts(TLorentzVector P,TAnaTrack *pl1, TAnaTrack *pl2, bool eta_cut){

	//Compute the seagull efficiency = Nmc(seagull && PS)/Nmc(PS)
	
	bool print=0; 
	
	//((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Cowboy Veto",1);
	if(pl1==0 || pl2 ==0){
		
		cout << "Bad track. MC_cuts." << endl;
		return 0; 
	}
	if(print)cout  << "xsReader::PS_cuts " << endl; 
   
	double mass=P.M();
	double ptCand = P.Pt();
	double candeta = P.Eta();	
	double candphi = P.Phi();
	
	TAnaTrack *pTra1=pl1; 
	TAnaTrack *pTra2=pl2;
	
	
	//Various muon cuts 
	
	if(print){
		cout << "mass: " << mass << endl;
		pl1->fPlab.Print(); 
		pl2->fPlab.Print();
		P.Print(); 
		
	}
	
	if ( (pl1->fPlab.Perp() > PTHI) || (pl1->fPlab.Perp() < PTLO) ) return 0;
	if ( (pl2->fPlab.Perp() > PTHI) || (pl2->fPlab.Perp() < PTLO) ) return 0;
	
	
	if(ETAHI>0 && !eta_cut){
		if ( TMath::Abs(pl1->fPlab.Eta()) > ETAHI || TMath::Abs(pl2->fPlab.Eta())>ETAHI)  {
				return 0;}
	}
	
	///////Muon in the barrel 
	if ( ((TMath::Abs(pl1->fPlab.Eta()) <= ETABARREL) && (pl1->fPlab.Perp() < PTBARREL)) || ((TMath::Abs(pl2->fPlab.Eta()) <= ETABARREL) && (pl2->fPlab.Perp() < PTBARREL)) ){

		return 0; 
	}
	
	if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETABARREL) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAMI) && (pl1->fPlab.Perp() < PTMI))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETABARREL) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAMI) && (pl2->fPlab.Perp() < PTMI)) ){
		return 0; 
	}
	
	if(!eta_cut){
		if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETAMI) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAHI) && (pl1->fPlab.Perp() < PTLO))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETAMI) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAHI) && (pl2->fPlab.Perp() < PTLO)) ){
			return 0; 
		}
	}
	else {
		if (( (TMath::Abs(pl1->fPlab.Eta()) >= ETAMI) && ((pl1->fPlab.Perp() < PTLO))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETAMI) && (pl2->fPlab.Perp() < PTLO)) ){
			return 0; 
		}
	}

		
	//Various candidate kinematic cuts 
	if ( P.Rapidity() < -RAPCAND ){ 
		return 0;}
	
	if ( P.Rapidity() > RAPCAND ){
		return 0;}
		
	if ( P.Perp() < PTCAND_MIN ){ 
		return 0;}
	
	if ( P.Perp() > PTCAND ){
		return 0;}

    if (mass < MASSLO) return 0;
    if (mass > MASSHI) return 0;
    
	return 1;
	
}

bool xsReader::getVectors(){
	bool print=0; 	
	
  TAnaTrack *pTrack1(0);
  TAnaTrack *pTrack2(0);

  TLorentzVector P1;
  TLorentzVector P2;
  TLorentzVector P_Ups;

  bool result=0;  

  ((TH1I*)fpHistFile->FindObjectAny("MCTracks"))->Fill(fpEvt->nRecTracks());
  
  if ( fpEvt->nRecTracks() != 2 ) return 0;
  else{
  for (int iR = 0; iR < fpEvt->nRecTracks(); ++iR) {
    
    if ( iR == 0 ){
      result =1; 
      pTrack1 = fpEvt->getRecTrack(iR);
      if(pTrack1->fPlab.Pt()<0)
	return 0;
      P1.SetPtEtaPhiM(pTrack1->fPlab.Pt(),pTrack1->fPlab.Eta(),pTrack1->fPlab.Phi(),MMUON);
    }
    if ( iR == 1 ){
      result =1; 
      pTrack2 = fpEvt->getRecTrack(iR);
      if(pTrack2->fPlab.Pt()<0)
	return 0;
      P2.SetPtEtaPhiM(pTrack2->fPlab.Pt(),pTrack2->fPlab.Eta(),pTrack2->fPlab.Phi(),MMUON);

    }
        
  }
  if(pTrack1==0 || pTrack2==0) return 0; 
  if(result!=1) return 0; // if it did not get tracks return 0. 
  //************************************************
  //  if(fGenMass>9.4605 || fGenMass<9.4599) return 0; //This is to look at just the gen mass peak
  //***********************************************
  P_Ups=P1+P2;  //Upsilon 4-Vector
  fCand4V = P_Ups;
  if(print)  fCand4V.Print();
  if(print) cout << endl;
  
  
  ((TH2D*)fpHistFile->FindObjectAny("dxy_candphi"))->Fill(pTrack1->fdxy,P_Ups.Phi());
  ((TH2D*)fpHistFile->FindObjectAny("dxy_candphi"))->Fill(pTrack2->fdxy,P_Ups.Phi());
  fpHistFile->cd(); 
  if(print) cout << "Before MC cuts. " << endl; 
  MCefficiency(fCand4V,pTrack1,pTrack2);	  
  if(MC_cuts(fCand4V, pTrack1, pTrack2)==0) return 0; 
  if(print) cout << "Passed MC cuts. " << endl; 
  if(CowboyVetoMC(pTrack1, pTrack2) )
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Cowboy Veto",1);
  fpHistFile->cd(); 
	 	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Cowboy Veto",1);
  if(print) cout << "Passed cuts up to result." << endl; 
  if(result==0) return 0; 
    
  	  fMuon1Pt = pTrack1->fPlab.Perp(); fMuon2Pt = pTrack2->fPlab.Perp(); 
	  fMuon1Eta = pTrack1->fPlab.Eta(); fMuon2Eta = pTrack2->fPlab.Eta();
	  fMuon1Phi = pTrack1->fPlab.Phi(); fMuon2Phi = pTrack2->fPlab.Phi();
	  
	  fMuon1Vect=pTrack1->fPlab;
	  fMuon2Vect=pTrack2->fPlab;
	  
	  fMuonPhi12 =pTrack1->fQ*(pTrack1->fPlab.DeltaPhi(pTrack2->fPlab));
	  
	  fMuon1PtE = pTrack1->fPtE;    fMuon2PtE = pTrack2->fPtE;
	  //fMuon1PtE=0.01*fMuon1Pt; fMuon2PtE=0.01*fMuon2Pt;
	  fMuon1EtaE = pTrack1->fEtaE;  fMuon2EtaE = pTrack2->fEtaE;
	  fMuon1PhiE = pTrack1->fPhiE;  fMuon2PhiE = pTrack2->fPhiE;


    fCandPt   = P_Ups.Pt();
    fCandMass = P_Ups.M();
    fCandY = TMath::Abs(P_Ups.Rapidity());//Take Abs of rapidity...
    
/*
    if(UPSTYPE==1){
    int index=((TH1D*)fWeightFile->FindObjectAny("Y1Pt"))->FindBin(fCandPt);
    fWeight=((TH1D*)fWeightFile->FindObjectAny("Y1Pt"))->FindObjectAnyBinContent(index);}
    if(UPSTYPE==2){
    int index=((TH1D*)fWeightFile->FindObjectAny("Y2Pt"))->FindBin(fCandPt);
    fWeight=((TH1D*)fWeightFile->FindObjectAny("Y2Pt"))->FindObjectAnyBinContent(index);}
    if(UPSTYPE==3){
    int index=((TH1D*)fWeightFile->FindObjectAny("Y3Pt"))->FindBin(fCandPt);
    fWeight=((TH1D*)fWeightFile->FindObjectAny("Y3Pt"))->FindObjectAnyBinContent(index);}
*/  


  }//else statement
  //if you get through the cuts, return 1.
	  fWeight=1;
	  return 1; 
}

void xsReader::freePointers(){
	/*
	 possibly faulty piece of pointer clearing code 
  while (!Cands.empty())
    {
      // remove it from the list
      Cands.erase(Cands.begin());
    }
  
  while (!Cands_ID.empty())
    {
      // remove it from the list
      Cands_ID.erase(Cands_ID.begin());
    } 
  
  while (!Cands_TM.empty())
    {
      // remove it from the list
      Cands_TM.erase(Cands_TM.begin());
    }  
  */

  Cands.clear();
  Cands_TM.clear();
  Cands_ID.clear();	 
}


void xsReader::plot_kinematics(){
	
  fpHistFile->cd("kinematics");
	((TH1D*)fpHistFile->FindObjectAny("dphi"))->Fill(fMuonPhi12);
	((TH1D*)fpHistFile->FindObjectAny("muonphi"))->Fill(fMuon1Vect.Phi());
	((TH1D*)fpHistFile->FindObjectAny("muonphi"))->Fill(fMuon2Vect.Phi());
	((TH2D*)fpHistFile->FindObjectAny("phi1_phi2"))->Fill(fMuon1Vect.Phi(),fMuon2Vect.Phi());
	((TH2D*)fpHistFile->FindObjectAny("pt_phi"))->Fill(fMuonPhi12, fCandPt);
	((TH2D*)fpHistFile->FindObjectAny("muonPxPy"))->Fill(fMuon1Vect.Px(),fMuon1Vect.Py());
	((TH2D*)fpHistFile->FindObjectAny("genMass_RecoMass"))->Fill(fGenMass,fCandMass);
	
	if(fCandMass>9.3 && fCandMass<10.3 && fCandY<fYbin2[1]){
		for(int iy =0; iy<fNy; ++iy){
			if(fMuon1Vect.Eta()>=fYbin[iy] && fMuon1Vect.Eta()<fYbin[iy+1]) {
				((TH1D*)fpHistFile->FindObjectAny(Form("muonPt_rapidity_%dS%.1f_%.1f", UPSTYPE, fYbin[iy], fYbin[iy+1])))->Fill(fMuon1Vect.Perp());
			}
			if(fMuon2Vect.Eta()>=fYbin[iy] && fMuon2Vect.Eta()<fYbin[iy+1]) {
				((TH1D*)fpHistFile->FindObjectAny(Form("muonPt_rapidity_%dS%.1f_%.1f", UPSTYPE, fYbin[iy], fYbin[iy+1])))->Fill(fMuon2Vect.Perp());
			}
			
		}//loop over pseudo-rapidity bins 
		
	}//if 9.3<M<10.3
	if(fCandY<fYbin2[1]){
		if(fCandMass>9.3 && fCandMass<10.3){
			((TH2D*)fpHistFile->FindObjectAny("MuonPt_Candpt"))->Fill(fCandPt,fMuon1Vect.Perp());
			((TH2D*)fpHistFile->FindObjectAny("MuonPt_Candpt"))->Fill(fCandPt,fMuon2Vect.Perp());
		}
		if(fCandMass>9.26 && fCandMass<9.66) ((TH2D*)fpHistFile->FindObjectAny("EventsY1_y_pt"))->Fill(fCandPt, fCandY);
		if(fCandMass>9.92 && fCandMass<10.12) ((TH2D*)fpHistFile->FindObjectAny("EventsY2_y_pt"))->Fill(fCandPt, fCandY);
		if(fCandMass>10.255 && fCandMass<10.455) ((TH2D*)fpHistFile->FindObjectAny("EventsY3_y_pt"))->Fill(fCandPt,fCandY);
	}
	((TH2D*)fpHistFile->FindObjectAny("muon_kinematics"))->Fill(fMuon1Vect.Perp(),fMuon1Vect.Eta());
	((TH2D*)fpHistFile->FindObjectAny("muon_kinematics"))->Fill(fMuon2Vect.Perp(),fMuon2Vect.Eta());
	if(fCandY<fYbin2[1]){
	  ((TH2D*)fpHistFile->FindObjectAny("muon_kinematics_y0.0_0.6"))->Fill(fMuon1Vect.Perp(),fMuon1Vect.Eta());
	  ((TH2D*)fpHistFile->FindObjectAny("muon_kinematics_y0.0_0.6"))->Fill(fMuon2Vect.Perp(),fMuon2Vect.Eta());
	}

	if(fCandMass>9.26 && fCandMass<9.66){
		((TH2D*)fpHistFile->FindObjectAny("muon_pt_upsilon1S_pt"))->Fill(fCandPt, fMuon1Vect.Perp());
		((TH2D*)fpHistFile->FindObjectAny("muon_pt_upsilon1S_pt"))->Fill(fCandPt, fMuon2Vect.Perp());
	}
	double angle=fMuon1Vect.Angle(fMuon2Vect);
	((TH1D*)fpHistFile->FindObjectAny("opening_angle"))->Fill(angle); 


	
   fpHistFile->cd();
}

void xsReader::plot_errors(){
        double dphi=fMuonPhi12;
	double dphiE=sqrt(fMuon1PhiE*fMuon1PhiE+fMuon2PhiE*fMuon2PhiE);

	fpHistFile->cd("errors");
	((TH2D*)fpHistFile->FindObjectAny("muonptE_pt"))->Fill(fMuon1Vect.Perp(),fMuon1PtE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_pt"))->Fill(fMuon2Vect.Perp(),fMuon2PtE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_eta"))->Fill(TMath::Abs(fMuon1Eta),fMuon1PtE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_eta"))->Fill(TMath::Abs(fMuon2Eta),fMuon2PtE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_dphi"))->Fill(TMath::Abs(dphi), fMuon1PtE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_dphi"))->Fill(TMath::Abs(dphi), fMuon2PtE*1000,fWeight);
	
	((TH2D*)fpHistFile->FindObjectAny("muonetaE_pt"))->Fill(fMuon1Vect.Perp(),fMuon1EtaE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonetaE_pt"))->Fill(fMuon2Vect.Perp(),fMuon2EtaE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonetaE_eta"))->Fill(TMath::Abs(fMuon1Eta),fMuon1EtaE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonetaE_eta"))->Fill(TMath::Abs(fMuon2Eta),fMuon2EtaE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonetaE_dphi"))->Fill(TMath::Abs(dphi), fMuon1EtaE*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muonetaE_dphi"))->Fill(TMath::Abs(dphi), fMuon2EtaE*1000,fWeight);


	((TH2D*)fpHistFile->FindObjectAny("muondphiE_pt"))->Fill(fMuon1Pt,TMath::Abs(dphiE)*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muondphiE_pt"))->Fill(fMuon2Pt,TMath::Abs(dphiE)*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muondphiE_eta"))->Fill(fMuon1Eta,TMath::Abs(dphiE)*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muondphiE_eta"))->Fill(fMuon2Eta,TMath::Abs(dphiE)*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("muondphiE_dphi"))->Fill(TMath::Abs(dphi),TMath::Abs(dphiE)*1000,fWeight);

	((TH1D*)fpHistFile->FindObjectAny("ptE_pt"))->Fill((fMuon1PtE/fMuon1Vect.Perp())*100,fWeight);
	((TH1D*)fpHistFile->FindObjectAny("ptE_pt"))->Fill((fMuon2PtE/fMuon2Vect.Perp())*100,fWeight);

	 for(int iy =0; iy<fNy; ++iy){
	   if(TMath::Abs(fMuon1Eta)>=fYbin[iy] && TMath::Abs(fMuon1Eta)<fYbin[iy+1])
	     ((TH1D*)fpHistFile->FindObjectAny(Form("ptE_pt_etabin%.1f_%.1f", fYbin[iy],fYbin[iy+1])))->Fill(fMuon1PtE*100/fMuon1Pt,fWeight);
	   if(TMath::Abs(fMuon2Eta)>=fYbin[iy] && TMath::Abs(fMuon2Eta)<fYbin[iy+1])
	     ((TH1D*)fpHistFile->FindObjectAny(Form("ptE_pt_etabin%.1f_%.1f", fYbin[iy],fYbin[iy+1])))->Fill(fMuon2PtE*100/fMuon2Pt,fWeight);

                  }
	 if(TMath::Abs(fMuon1Eta)<0.8)
	   ((TH1D*)fpHistFile->FindObjectAny("ptE_pt_etabin0_0.8"))->Fill(fMuon1PtE*100/fMuon1Pt,fWeight);
	 if(TMath::Abs(fMuon2Eta)<0.8)
	   ((TH1D*)fpHistFile->FindObjectAny("ptE_pt_etabin0_0.8"))->Fill(fMuon2PtE*100/fMuon2Pt,fWeight);
	 if(TMath::Abs(fMuon1Eta)>0.8 && TMath::Abs(fMuon1Eta<2.5))
	   ((TH1D*)fpHistFile->FindObjectAny("ptE_pt_etabin0.8_2.5"))->Fill(fMuon1PtE*100/fMuon1Pt,fWeight);
	 if(TMath::Abs(fMuon2Eta)>0.8 && TMath::Abs(fMuon2Eta)<2.5)
	   ((TH1D*)fpHistFile->FindObjectAny("ptE_pt_etabin0.8_2.5"))->Fill(fMuon2PtE*100/fMuon2Pt,fWeight);


	((TH2D*)fpHistFile->FindObjectAny("ptE_pt_eta"))->Fill(TMath::Abs(fMuon1Eta),fMuon1PtE*100/fMuon1Pt,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("ptE_pt_eta"))->Fill(TMath::Abs(fMuon2Eta),fMuon2PtE*100/fMuon2Pt,fWeight);

	((TH1D*)fpHistFile->FindObjectAny("etaE_eta"))->Fill(TMath::Abs((fMuon1EtaE/fMuon1Vect.Eta()*100)),fWeight);
	((TH1D*)fpHistFile->FindObjectAny("etaE_eta"))->Fill(TMath::Abs((fMuon2EtaE/fMuon2Vect.Eta()*100)),fWeight);
	((TH1D*)fpHistFile->FindObjectAny("dphiE_phi"))->Fill(TMath::Abs(dphiE/dphi)*100);

	if(fCandY<fYbin2[1])
		((TH2D*)fpHistFile->FindObjectAny("sigma_pt"))->Fill(fGenPt,fCandPtE); //sigmaPt vs Pt
	
}

void xsReader::print_event_details(){

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
	

	 cout << "candPt: " << fCandPt << " CandMass: " << fCandMass <<  endl;      
	 cout << "Pt1: " << pt1 << " +/- " << fMuon1PtE << endl;
	 cout << "Pt2: " << pt2 << " +/- " << fMuon2PtE << endl;
	 cout << "Eta1: " << fMuon1Eta << " Eta1 from Vector: " << fMuon1Vect.Eta() << " +/- " << fMuon1EtaE << endl;
	 cout << "Eta2: " << fMuon2Eta << " Eta2 from Vector: " << fMuon2Vect.Eta() << " +/- " << fMuon2EtaE << endl;
		 
	 cout << "dphi: " << dphi << " +/- " << dphiE << endl;
	 
	 cout << "P1: " << sqrt(pt1*pt1+Pl1*Pl1) << " " << fMuon1Vect.Mag()  << " P2: " << sqrt(pt2*pt2+Pl2*Pl2) << " " << fMuon2Vect.Mag() << endl;  
		 
	 cout << "E1: " << E1 << " +/- " << sigmaE1 << endl;
	 cout << "E2: " << E2 << " +/- " << sigmaE2 << endl;
	 

	 
	 //cout << endl; 
	
	 // cout << "dphi: " << dphi << endl;

	cout << "Polar Coordinates: " << endl; 

	cout << "P1 <Pt,Eta,Phi>: " << fMuon1Vect.Perp() << " " << fMuon1Vect.Eta() << " " << fMuon1Vect.Phi() <<endl;
        cout << "P2 <Pt,Eta,Phi>: " << fMuon2Vect.Perp() << " " << fMuon2Vect.Eta() << " " << fMuon2Vect.Phi() << endl;

	cout << "Cartesian coordintes:" << endl;
 
        cout << "P1 <x,y,z>: " << fMuon1Vect.Px() << " " << fMuon1Vect.Py() << " " << fMuon1Vect.Pz() << endl;   
	cout << "P2 <x,y,z>: " << fMuon2Vect.Px() << " " << fMuon2Vect.Py() << " " << fMuon2Vect.Pz() << endl;

	cout << endl; 

}

double xsReader::correction_functionJPsi(const double & pt, const double & eta, const double & phi, const int chg, double *par) {
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

double xsReader::correction_functionZ(const double & pt, const double & eta, const double & phi, const int chg,   double *parScale) {    
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

double xsReader::Pt_primeZ(TVector3 v1,int Q=1){
	//Z corrections on 2011 Data 
	double parScale[]={-0.000322,0.001403,0.98565,-0.000792,-2.1,0.000535,0.825151,-0.000234,-1.5,0.000148,-1.47645,-2.8e-05,1.5,0.000182,0.84186,-0.000425,2.1,0.000733,1.93751,-0.001493,4.8e-05,0.00048,2,1.79589,0.000586,2,-0.518118};
	double parScale_MC[]={0.00123,0.001113,0.168873,0.000124,-2.1,0.000336,-0.999195,0.000174,-1.5,0.000177,-1.76348,2.6e-05,1.5,0.000377,-1.08773,-0.000101,2.1,0.000879,-1.43165,6.4e-05,4.4e-05,0.00087,2,-1.94623,0.000651,2,2.54445};
	
	double muon_pt=v1.Perp();
	double muon_eta=v1.Eta();
	double muon_phi=v1.Phi();
	double corrected_pt=0; 
	if(MODE==4){
		//DATA
		corrected_pt=correction_functionZ(muon_pt,muon_eta,muon_phi,Q,parScale);
	}
	if(MODE==5){
		corrected_pt=correction_functionZ(muon_pt,muon_eta,muon_phi,Q,parScale_MC);
	}
	return corrected_pt; 
	
}

double xsReader::Pt_primeJPsi(TVector3 v1,int Q=1){
  //Momentum correction using JPsi, Pt<10GeV, maybe extrapolate to 20GeV? 
 
//Mode==4 (data)	
  double parScaleJPsi_it1[]={1.00074,0.000103333, -0.000947069, -0.00723064, 1.65256, 0., 0.000181535, 6}; // JPsi full 2011 data set, corrections for first iteration
  double parScaleJPsi_it2[]={1.00014, 2.2114e-05, -6.29854e-05,-0.000515906, 1.65164, 0., -7.40772e-06, 6}; // JPsi corrections for 2nd iteration
  
	double parScaleJPsi_MC[]={1.00006,-0.00042742,0.00011569,0.00393103,1.60458,0,5.87978e-05,6};

  double muon_pt=v1.Perp();
  double muon_eta=v1.Eta();
  double muon_phi=v1.Phi();
  double corrected_pt=0; 
  if(MODE==4){
	  //Data
		corrected_pt=correction_functionJPsi(muon_pt,muon_eta,muon_phi,Q,parScaleJPsi_it1);
		corrected_pt=correction_functionJPsi(corrected_pt,muon_eta,muon_phi,Q,parScaleJPsi_it2); 
	}
  if(MODE==5){
	  corrected_pt=correction_functionJPsi(muon_pt,muon_eta,muon_phi,Q,parScaleJPsi_MC);
	}
	
  return corrected_pt; 
}

void xsReader::correct_Tracks(){
	TAnaCand *pCand(0); 
	TAnaTrack *pTrack1(0); 
	TAnaTrack *pTrack2(0);
	for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
		pCand = fpEvt->getCand(iC);
		pTrack1 = fpEvt->getSigTrack(pCand->fSig1); 
		pTrack2 = fpEvt->getSigTrack(pCand->fSig2);
		
		double corrected_pt1=0;
		double corrected_pt2=0; 

		TVector3 p1=pTrack1->fPlab;
		int Q1=pTrack1->fQ; 
		double Pt1_before=p1.Perp();

		TVector3 p2=pTrack2->fPlab;
		int Q2=pTrack2->fQ;
		double Pt2_before=p2.Perp();

		TLorentzVector P1_4V;
		TLorentzVector P2_4V;
		TLorentzVector fPUps;
		P1_4V.SetVectM(p1,MMUON);
		P2_4V.SetVectM(p2,MMUON);
		fPUps=P1_4V+P2_4V;

		double mass_before=fPUps.M();
		if(p1.Perp()<JPSI_Z_CORR){
			corrected_pt1=Pt_primeJPsi(p1,Q1);
			
		}
		else {
			corrected_pt1=Pt_primeZ(p1,Q1);
			
		}

		if(p2.Perp()<JPSI_Z_CORR)
			corrected_pt2=Pt_primeJPsi(p2,Q2);
		else {
			corrected_pt2=Pt_primeZ(p2,Q2);
		}


		p1.SetPerp(corrected_pt1); 
		p2.SetPerp(corrected_pt2); 

		P1_4V.SetVectM(p1,MMUON);
		P2_4V.SetVectM(p2,MMUON);
		fPUps=P1_4V+P2_4V;
  
		double mass_after=fPUps.M();
		fpHistFile->cd("shifts");
		if(corrected_pt1 <20 && corrected_pt2<20){
			((TH1D*)fpHistFile->FindObjectAny("Mass_Shift_low"))->Fill((mass_after-mass_before)*1000); 
			((TH1D*)fpHistFile->FindObjectAny("Pt_Shift_low"))->Fill((corrected_pt1-Pt1_before)*1000);
			((TH1D*)fpHistFile->FindObjectAny("Pt_Shift_low"))->Fill((corrected_pt2-Pt2_before)*1000);
		}
		
		else {
			((TH1D*)fpHistFile->FindObjectAny("Mass_Shift_high"))->Fill((mass_after-mass_before)*1000); 
			((TH1D*)fpHistFile->FindObjectAny("Pt_Shift_high"))->Fill((corrected_pt1-Pt1_before)*1000);
			((TH1D*)fpHistFile->FindObjectAny("Pt_Shift_high"))->Fill((corrected_pt2-Pt2_before)*1000);
		}
		
		
		for(int ipt=1; ipt<fNmass; ipt++){
			if((corrected_pt1>fPTbin[ipt] && corrected_pt1<fPTbin[ipt+1]) || (corrected_pt2>fPTbin[ipt] && corrected_pt2<fPTbin[ipt+1])){
				 ((TH1D*)fpHistFile->FindObjectAny(Form("Mass_Shift_Pt_%.1f_%.1f",fPTbin[ipt],fPTbin[ipt+1])))->Fill((mass_after-mass_before)*1000);
				 ((TH1D*)fpHistFile->FindObjectAny(Form("Pt_Shift_Pt_%.1f_%.1f",fPTbin[ipt], fPTbin[ipt+1])))->Fill((corrected_pt1-Pt1_before)*1000);
				 ((TH1D*)fpHistFile->FindObjectAny(Form("Pt_Shift_Pt_%.1f_%.1f",fPTbin[ipt], fPTbin[ipt+1])))->Fill((corrected_pt2-Pt2_before)*1000);
				 }//if statement 
			}// ipt loop 
		fpHistFile->cd();
		
		//cout << "Delta Pt: " << (corrected_pt1-Pt1_before)*1000 << endl; 
		
		pTrack1->fPlab.SetPerp(corrected_pt1);
		pTrack2->fPlab.SetPerp(corrected_pt2);
		
		pCand->fPlab.SetPtEtaPhi(fPUps.Perp(),fPUps.Eta(),fPUps.Phi()); 
		pCand->fMass=fPUps.M();
		
	}//iC loop 
}

void xsReader::Pt_uncertainty(){
	//compute the error on the dimuon pt from track errors 
	
	double phi12=fMuon1gen.DeltaPhi(fMuon2gen);
	
	double D1=fMuon1gen.Perp()+fMuon2gen.Perp()*TMath::Cos(phi12);//Derivative 1
	double D2=fMuon2gen.Perp()+fMuon1gen.Perp()*TMath::Cos(phi12);
	double D3=fMuon1gen.Perp()*fMuon2gen.Perp()*TMath::Sin(phi12);
	
	double Phi12E=TMath::Sqrt(fMuon1PhiE*fMuon1PhiE+fMuon2PhiE*fMuon1PhiE); //Delta Phi error 
	
	double T1=TMath::Power(D1*fMuon1PtE,2); //Term 1
	double T2=TMath::Power(D2*fMuon2PtE,2); //Term 2
	double T3=TMath::Power(D3*Phi12E,2); //Term3
	
	fCandPtE=(1/fGenPt)*TMath::Sqrt(T1+T2+T3);
	
	//cout << "genPt: " << fGenPt << " fCandPtE: " << fCandPtE << endl; 
	
}

double xsReader::M_uncertainty(){
	bool print =0; 
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
	if(print) cout <<"M: " << fCandMass << " sigmaM: " << sigmaM << endl; 
	if(print) cout << "fCandPt: " << fCandPt << " y: " << fCandY << endl; 
	((TH2D*)fpHistFile->FindObjectAny("deta_dm"))->Fill(fMuon2Vect.Eta()-fMuon1Vect.Eta(),sigmaM*1000,fWeight);
	((TH2D*)fpHistFile->FindObjectAny("phi_dm"))->Fill(dphi,sigmaM*1000);
	((TH2D*)fpHistFile->FindObjectAny("DR"))->Fill(sqrt(pow((fMuon2Vect.Eta()-fMuon1Vect.Eta()),2)+dphi*dphi),sigmaM*1000,fWeight);
	
	   double Y_abs=TMath::Abs(fCandY);
	   double fCandEta = TMath::Abs(fCand4V.Eta());

     //specific dm binning


     // same thing, but 2D, m-dm histograms now. 

	
	for ( int iy = 0; iy < fNy2; ++iy ){
		for ( int ipt = 0; ipt < fNpt1; ++ipt ){
			if(fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1] && fCandPt>=fPTbin1[ipt] && fCandPt<fPTbin1[ipt+1]){
				if(print) cout << "Filling Histogram: " << endl; 
				((TH2D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->Fill(fCandMass,sigmaM*1000,fWeight);
			}//end if statement
		}//pt loop
	}//y loop
	
	if(MODE!=5){
		//only weight by 1/eff for data, not MC
		//First fill histograms with fPTbin1
		
		
		for ( int iy = 0; iy < fNy2; ++iy ){
			for ( int ipt = 0; ipt < fNpt1; ++ipt ){
				if(fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1] && fCandPt>=fPTbin1[ipt] && fCandPt<fPTbin1[ipt+1]){
					if(print) cout << "Filling Histogram: " << endl; 
					((TH2D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effw_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->Fill(fCandMass,sigmaM*1000,1./dimuon_efficiency);
				}//end if statement
			}//pt loop
		}//y loop
		
		
		for ( int iy = 0; iy < fNy2; ++iy ){
			for ( int ipt = 0; ipt < fNpt1; ++ipt ){
				if(fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1] && fCandPt>=fPTbin1[ipt] && fCandPt<fPTbin1[ipt+1]){
					if(print) cout << "Filling Histogram: " << endl; 
					((TH2D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effwEp_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->Fill(fCandMass,sigmaM*1000,1./dimuon_efficiencyEp);
				}//end if statement
			}//pt loop
		}//y loop	
		
		
		
		for ( int iy = 0; iy < fNy2; ++iy ){
			for ( int ipt = 0; ipt < fNpt1; ++ipt ){
				if(fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1] && fCandPt>=fPTbin1[ipt] && fCandPt<fPTbin1[ipt+1]){
					if(print) cout << "Filling Histogram: " << endl; 
					((TH2D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effwEm_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->Fill(fCandMass,sigmaM*1000,1./dimuon_efficiencyEm);
				}//end if statement
			}//pt loop
		}//y loop			
		

	}//Run over Data only 
	
	
	  for(int iy =0; iy<fNy2; ++iy){
	    if(fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1])
	      ((TH2D*)fpHistFile->FindObjectAny(Form("dm_pt_ybin%.1f_%.1f", fYbin2[iy],fYbin2[iy+1])))->Fill(fCandPt, sigmaM*1000,fWeight);

	  }

	for(int ipt=0; ipt<fNpt2; ++ipt){
	    if(fCandPt>=fPTbin2[ipt] && fCandPt<fPTbin2[ipt+1])
	      ((TH2D*)fpHistFile->FindObjectAny(Form("dm_y_ptbin_%.1f_%.1f",fPTbin2[ipt],fPTbin2[ipt+1])))->Fill(fCandY,sigmaM*1000,fWeight);
	  }
	
	for(int ipt=0; ipt<fNpt1; ++ipt){
	    if(fCandPt>=fPTbin1[ipt] && fCandPt<fPTbin1[ipt+1])
			((TH2D*)fpHistFile->FindObjectAny(Form("dm_y_ptbin1_%.1f_%.1f",fPTbin1[ipt],fPTbin1[ipt+1])))->Fill(fCandY,sigmaM*1000,fWeight);
	}
	

	
	for (int iy=0; iy<fNy2; iy++) {
		if(fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1]){
			((TH2D*)fpHistFile->FindObjectAny(Form("m_dm_y%d",iy)))->Fill(fCandMass,sigmaM*1000,fWeight);
			((TH2D*)fpHistFile->FindObjectAny(Form("m_dmw_y%d",iy)))->Fill(fCandMass,sigmaM*1000,1/dimuon_efficiency);
			((TH2D*)fpHistFile->FindObjectAny(Form("m_dmw_Ep_y%d",iy)))->Fill(fCandMass,sigmaM*1000,1/dimuon_efficiencyEp);
			((TH2D*)fpHistFile->FindObjectAny(Form("m_dmw_Em_y%d",iy)))->Fill(fCandMass,sigmaM*1000,1/dimuon_efficiencyEm);
		}
	}
	
	if(fCandY<fYbin2[1]){
		((TH1D*)fpHistFile->FindObjectAny("dm"))->Fill(sigmaM*1000);
		
		
	}
	
       ((TH2D*)fpHistFile->FindObjectAny("pt_dm"))->Fill(TMath::Abs(fCandPt),sigmaM*1000,fWeight);
       ((TH2D*)fpHistFile->FindObjectAny("y_dm"))->Fill(TMath::Abs(fCandY),sigmaM*1000,fWeight);
	   ((TProfile2D*)fpHistFile->FindObjectAny("pt_y_dm"))->Fill(fCandPt,fCandY,sigmaM*1000);
	
	for ( int ipt = 0; ipt < fNpt2; ++ipt ){
	      if (  fCandPt >= fPTbin2[ipt]  && fCandPt < fPTbin2[ipt+1] ) 
		((TH2D*)fpHistFile->FindObjectAny(Form("m_dm_pt%.1f_%.1f", fPTbin2[ipt], fPTbin2[ipt+1])))->Fill(fCandMass,sigmaM*1000,fWeight);
		
		}
		
	fpHistFile->cd();
	//if( fCandPt>10 && TMath::Abs(Pl1-Pl2)<12)
	((TH1D*)fpHistFile->FindObjectAny("UpsilonMasscut"))->Fill(fCandMass,fWeight);

	if(fCandY<fYbin2[1])
		((TH2D*)fpHistFile->FindObjectAny("dm_nPv"))->Fill(fpEvt->nPV(),sigmaM*1000); 

	return sigmaM; 
	
}

bool xsReader::CowboyVeto(TAnaCand *pCand){
  bool veto =false;
  TAnaTrack *pl1(0); TAnaTrack *pl2(0);
  pl1 = fpEvt->getSigTrack(pCand->fSig1); 
  pl2 = fpEvt->getSigTrack(pCand->fSig2);
  
  //if ( (pl1->fQ*( pl1->fPlab.Phi() - pl2->fPlab.Phi() )) > 0 ) veto=true;  //old method of computing delta phi
 
  if(  (pl1->fQ*(  pl1->fPlab.DeltaPhi(pl2->fPlab )))  >0)  veto = true;   // new method, which accounts for interval on which phi is defined
   return veto;
}

bool xsReader::CowboyVetoMC(TAnaTrack *pl1, TAnaTrack *pl2){
 bool veto =false;
 // cout << "CowboyVetoMC." << endl; 
  if(  (pl1->fQ*(  pl1->fPlab.DeltaPhi(pl2->fPlab )))  >0)  veto = true;   // new method, which accounts for interval on which phi is defined
   return veto;
}

void xsReader::UpsGun_acceptance(){
	  //Only run UpsGun_acceptance of MC - never Data 
	
	  if(MODE!=5) return; 
	  //acceptance calculating function, run over MC 
	
	  TGenCand *gCand(0); TAnaTrack *pTrack(0); 
	  TGenCand *gDau1(0); TGenCand *gDau2(0);  
	  double pt, rapidity; bool fill = false; 
	  bool ups = false; bool mu1 = false; bool mu2 = false; 
	  double pt1(-1.), pt2(-1.);
	  double eta1(-99.), eta2(-99);
	  int index1(-99), index2(-99);
	  double w1(-99); 
	  double w1Ep(-99); 
	  double w1Em(-99);
	  TLorentzVector genCand; TAnaMuon *pMuon;
	  TGenCand *g2Cand; TGenCand *g2Cand_; TGenCand *gUps; TGenCand *gMu1; TGenCand *gMu2;
	  int m(0);
	  int mp;  
	  TLorentzVector genMuPlus;
	  Float_t cosTheta; Float_t sinTheta; Float_t sin2Theta; Float_t cosPhi; Float_t cos2Phi;
	  TLorentzVector h1; h1.SetPxPyPzE(0,1,0,0);
	  double lambdaTheta(-99); double lambdaPhi(-99); double lambdaThetaPhi(-99); 
	  double lambdaTheta_Ep(-99); double lambdaPhi_Ep(-99); double lambdaThetaPhi_Ep(-99); //polarization parameter 1sigma Error +
	  double lambdaTheta_Em(-99); double lambdaPhi_Em(-99); double lambdaThetaPhi_Em(-99); // 1 sigma error -
	  
	  //fpEvt->dumpGenBlock();
	  //cout << " fpEvt->nGenCands() = " << fpEvt->nGenCands() << endl;
	  for (int iG = 0; iG < fpEvt->nGenCands(); ++iG) {
		  gCand = fpEvt->getGenCand(iG);
		  if ( gCand->fID == RESTYPE && gCand->fStatus == 2 ) {
			  ups = true;
			  gUps = fpEvt->getGenCand(iG);
		  }
		  if ( gCand->fID == -13 && gCand->fStatus == 1 ) {
			  mu1 = true;
			  gMu1 = fpEvt->getGenCand(iG);
		  }
		  if ( gCand->fID == 13 && gCand->fStatus == 1 ) {
			  mu2 = true;
			  gMu2 = fpEvt->getGenCand(iG);
			  
			  if ( ups ) {
				  genCand.SetPtEtaPhiE(gUps->fP.Perp(),gUps->fP.Eta(),gUps->fP.Phi(),gUps->fP.Energy());
				  TLorentzRotation boost(-genCand.BoostVector()); // For different Polarizations
				  
				  // calculate cosTheta Helicity
				  genMuPlus.SetPtEtaPhiM(gMu2->fP.Perp(), gMu2->fP.Eta(), gMu2->fP.Phi(), MMUON);
				  genMuPlus *= boost;
				  cosTheta = genMuPlus.Vect().Dot(genCand.Vect())/(genMuPlus.Vect().Mag()*genCand.Vect().Mag());
				  sinTheta = genMuPlus.Vect().Cross(genCand.Vect()).Mag()/(genMuPlus.Vect().Mag()*genCand.Vect().Mag());
				  sin2Theta = 2*sinTheta*cosTheta;
				  
				  h1 *= boost; //h2 *= boost; h3 *= boost;
				  cosPhi = genMuPlus.Vect().Dot(h1.Vect().Unit())/genMuPlus.Vect().Mag();
				  cos2Phi = 2*cosPhi*cosPhi - 1;
				  
				  double error_SF=1.0; 
				  if(gUps->fP.Perp()>50) error_SF=3.0; // if Pt>50 GeV, multiply errors by factor of 3 
				  
				  if ( UPSTYPE == 1  ){
					  //effD gets the nominal value in a given bin
					  lambdaTheta = fPidTable1SLambdaThetaPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi = fPidTable1SLambdaPhiPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi = fPidTable1SLambdaThetaPhiPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  //errD gets the 1 sigma error in a given bin, which means to get the actual value of lambda_1sigma we need
					  //lambda_1sigma=lambda+1sigma
					  lambdaTheta_Ep = lambdaTheta + error_SF*fPidTable1SLambdaThetaPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi_Ep = lambdaPhi + error_SF*fPidTable1SLambdaPhiPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi_Ep = lambdaThetaPhi + error_SF*fPidTable1SLambdaThetaPhiPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
					  //lembda is the same for pos / negative pidtables 
					  
					  lambdaTheta_Em = lambdaTheta - error_SF*fPidTable1SLambdaThetaNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi_Em = lambdaPhi - error_SF*fPidTable1SLambdaPhiNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi_Em = lambdaThetaPhi - error_SF*fPidTable1SLambdaThetaPhiNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);	    
					  
				  }
				  
				  if ( UPSTYPE == 2  ){
					  lambdaTheta = fPidTable2SLambdaThetaPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi = fPidTable2SLambdaPhiPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi = fPidTable2SLambdaThetaPhiPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
					  lambdaTheta_Ep = lambdaTheta + error_SF*fPidTable2SLambdaThetaPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi_Ep = lambdaPhi + error_SF*fPidTable2SLambdaPhiPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi_Ep = lambdaThetaPhi + error_SF*fPidTable2SLambdaThetaPhiPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
					  //lembda is the same for pos / negative pidtables 
					  
					  lambdaTheta_Em = lambdaTheta - error_SF*fPidTable2SLambdaThetaNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi_Em = lambdaPhi - error_SF*fPidTable2SLambdaPhiNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi_Em = lambdaThetaPhi - error_SF*fPidTable2SLambdaThetaPhiNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
				  }
				  
				  if ( UPSTYPE == 3  ){
					  lambdaTheta = fPidTable3SLambdaThetaPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi = fPidTable3SLambdaPhiPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi = fPidTable3SLambdaThetaPhiPos->effD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
					  lambdaTheta_Ep = lambdaTheta + error_SF*fPidTable3SLambdaThetaPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi_Ep = lambdaPhi + error_SF*fPidTable3SLambdaPhiPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi_Ep = lambdaThetaPhi + error_SF*fPidTable3SLambdaThetaPhiPos->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
					  //lembda is the same for pos / negative pidtables 
					  
					  lambdaTheta_Em = lambdaTheta - error_SF*fPidTable3SLambdaThetaNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaPhi_Em = lambdaPhi - error_SF*fPidTable3SLambdaPhiNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  lambdaThetaPhi_Em = lambdaThetaPhi - error_SF*fPidTable3SLambdaThetaPhiNeg->errD(gUps->fP.Perp(), genCand.Rapidity(), 0.);
					  
				  }
				  //Weight function From polarization parameters 
				  w1 = ((3/(4*TMath::Pi()))/(3+lambdaTheta))*(1+(lambdaTheta*cosTheta*cosTheta)+(lambdaPhi*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi*sin2Theta*cosPhi));
				  w1Ep = ((3/(4*TMath::Pi()))/(3+lambdaTheta_Ep))*(1+(lambdaTheta_Ep*cosTheta*cosTheta)+(lambdaPhi_Ep*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi_Ep*sin2Theta*cosPhi));
				  w1Em = ((3/(4*TMath::Pi()))/(3+lambdaTheta_Em))*(1+(lambdaTheta_Em*cosTheta*cosTheta)+(lambdaPhi_Em*sinTheta*sinTheta*cos2Phi)+(lambdaThetaPhi_Em*sin2Theta*cosPhi));
				  //weight, w+1sigma weight, w-1sigma weight 
				  
				  ((TH1D*)fpHistFile->FindObjectAny("Ups_phi_weighted"))->Fill(gUps->fP.Phi(),w1); 
				  ((TH1D*)fpHistFile->FindObjectAny("Ups_y_weighted"))->Fill(genCand.Rapidity(),w1); 
				  ((TH1D*)fpHistFile->FindObjectAny("Ups_pt_weighted"))->Fill(gUps->fP.Perp(),w1);
				   
				  if ( gUps->fP.Perp() < PTCAND_MIN || gUps->fP.Perp()>PTCAND ) return; // Upsilon kinematic cuts
				  if (TMath::Abs(genCand.Rapidity()>RAPCAND)) return; 	
				  if (genCand.M()<MASSLO || genCand.M()>MASSHI) return; 
				  
				  ((TH1D*)fpHistFile->FindObjectAny("PolarizationWeight"))->Fill(w1);
				  
			  }
		  }
	  }
	  
	  if ( ups ){
		  //Fill Histograms
		  genCand.SetPtEtaPhiE(gUps->fP.Perp(),gUps->fP.Eta(),gUps->fP.Phi(),gUps->fP.Energy());
		  
				  
		  if ( (gUps->fP.Perp() <= PTCAND) && (TMath::Abs(genCand.Rapidity()) <= RAPCAND) ){
			  ((TH2D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_y_pt_%dS",UPSTYPE)))->Fill(TMath::Abs(genCand.Rapidity()), gUps->fP.Perp(), w1);
			  if(TMath::Abs(genCand.Rapidity())<fYbin2[1]){
				   ((TH1D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_%dS",UPSTYPE)))->Fill(gUps->fP.Perp(),w1); // normally need to weight by w1
				   ((TH1D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_Ep_%dS",UPSTYPE)))->Fill(gUps->fP.Perp(),w1Ep); // + systematic 
				   ((TH1D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_Em_%dS",UPSTYPE)))->Fill(gUps->fP.Perp(),w1Em); // - negative systematic 
			  }//Fill weights histograms 
			  }

		  if ( mu1 && mu2 ){
			 			  
			  
			  pt1 = gMu1->fP.Perp();
			  eta1 = gMu1->fP.Eta();
			  pt2 = gMu2->fP.Perp();
			  eta2 = gMu2->fP.Eta();
		
			  if(TMath::Abs(genCand.Rapidity())<fYbin2[1] && gUps->fP.Perp()>fPTbin2[0] && gUps->fP.Perp()<fPTbin2[1]){
				  ((TH2D*)fpHistFile->FindObjectAny(Form("UG_pt1_pt2_%dS",UPSTYPE)))->Fill(pt2,pt1); 
				  if(pt1>PTBARREL && pt2>PTBARREL) ((TH2D*)fpHistFile->FindObjectAny(Form("UG_eta1_eta2_%dS",UPSTYPE)))->Fill(TMath::Abs(eta1),TMath::Abs(eta2));
			  }
			   //Done filling histograms 
			  fill=false; 
			  //check muon kinematics. muon pt>PTLO && muon |eta|<ETAHI
			  if ((pt1 >= PTLO) && (pt2 >= PTLO) && (TMath::Abs(eta1) <= ETAHI) && (TMath::Abs(eta2) <= ETAHI) && (pt1 <= PTHI) && (pt2 <= PTHI) ){
				  fill = true;
				  
				  if ( ((TMath::Abs(eta1) <= ETABARREL) && (pt1 < PTBARREL)) || ((TMath::Abs(eta2) <= ETABARREL) && (pt2 < PTBARREL)) ){
					  fill = false;
				  }
				  if ( ((TMath::Abs(eta1) >= ETABARREL) && ((TMath::Abs(eta1) <= ETAMI) && (pt1 < PTMI))) || ((TMath::Abs(eta2) >= ETABARREL) && (TMath::Abs(eta2) <= ETAMI) && (pt2 < PTMI)) ){
					  fill = false;
				  }
				  if ( ((TMath::Abs(eta1) >= ETAMI) && ((TMath::Abs(eta1) <= ETAHI) && (pt1 < PTLO))) || ((TMath::Abs(eta2) >= ETAMI) && (TMath::Abs(eta2) <= ETAHI) && (pt2 < PTLO)) ){
					  fill = false;
				  }
			  }
		  }
		 
		  if ( fill ){
			  if ( (gUps->fP.Perp() <= PTCAND) && (TMath::Abs(genCand.Rapidity()) <= RAPCAND) ){
				  ((TH2D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_y_pt_%dS",UPSTYPE)))->Fill(TMath::Abs(genCand.Rapidity()), gUps->fP.Perp(), w1);
				  if(TMath::Abs(genCand.Rapidity())<fYbin2[1]){
					  ((TH1D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_%dS",UPSTYPE)))->Fill(gUps->fP.Perp(), w1); 
					  ((TH1D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_Ep_%dS",UPSTYPE)))->Fill(gUps->fP.Perp(), w1Ep); 
					  ((TH1D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_Em_%dS",UPSTYPE)))->Fill(gUps->fP.Perp(), w1Em); 
				  }//Fill reco/passing events 
			  }
		  }//fill
	  }//ups

	
//	compute acceptance from these histograms UG_AllGenRes_1S->Divide(UG_RecoGenRes_1S)
	
}

void xsReader::Trigger_Study(){
	TTrgObj* pO;
	bool dimuon_trigger=false; 
	
	//cout << "---------------------------------------------------------------------- " << endl;
	
	//cout << "Found a total of " << fpEvt->nTrgObj() << " trigger objects, the trigger objects are: " << endl;
	
	for (int a=0; a<NHLT; ++a) {
		if(fpEvt->fHLTResult[a]==1){
			TString path=fpEvt->fHLTNames[a];
			if(path.Contains("Dimuon") && path.Contains("Upsilon_Barrel")) dimuon_trigger=true;
			
		}//loop over number of HLT paths
	}
			   
	bool L3filtered=false; 
	   if(dimuon_trigger){	   
		   for(int i=0; i<fpEvt->nTrgObj();i++){
			   pO=fpEvt->getTrgObj(i);
			   //if(pO->fLabel.Contains("L3Filtered") &&pO->fLabel.Contains("Upsilon") )cout <<pO->fLabel << endl; 
			   if(pO->fLabel.Contains("hltBarrelUpsilonL3Filtered") ) L3filtered=true; 
		   }//loop over trigger objects
	   }//dimuon trigger
	if(dimuon_trigger && !L3filtered) cout << "Trig: " << dimuon_trigger << " " << L3filtered << endl; 			
	if(!dimuon_trigger && L3filtered) cout << "Trig: " << dimuon_trigger << " " << L3filtered << endl; 
	/*
	cout << "----------Muons--------------" << endl;
	TAnaMuon *pM;
	for (int i = 0; i < fpEvt->nMuons(); ++i) {
		pM = fpEvt->getMuon(i);
		//if ((pM->fMuID & 2) ||  (pM->fMuID & 1) || (pM->fMuID & 4)){
			pM->dump();
		
	}//muon loop
	 */
	if(fpEvt->nCands()!=1 && fpEvt->nMuons()==2) cout << "nCand: " << fpEvt->nCands() << " nMuons: " << fpEvt->nMuons() << endl; 
	
}


int xsReader::SingleMu_Match(int pTmin){

		
	vector<int> fired;
	vector<bool> index1;
	vector<bool> index2;
	

	//cout << "SingleMu_Match with pT min: " << pTmin << endl; 	
	//cout << "look for labels: " << Form("hltSingleMu%dL3Filtered%d",pTmin,pTmin) << endl;
	//cout << "nTrig Objects: " << fpEvt->nTrgObj() << endl; 
	for(int i=0; i<fpEvt->nTrgObj();i++){
		
		TTrgObj *pTrig=fpEvt->getTrgObj(i);
		//cout << pTrig->fLabel << endl; 
		if(pTrig->fLabel.Contains(Form("hltSingleMu%dL3Filtered%d",pTmin,pTmin)) ) fired.push_back(i); 
		
		
	}//loop over trigger objects
	
	bool print=false; 
	//if(fired.size()>1) print=true;
	
	TAnaCand *pCand=fpEvt->getCand(0); 
	TAnaTrack *pl1=fpEvt->getSigTrack(pCand->fSig1); 
	TAnaTrack *pl2=fpEvt->getSigTrack(pCand->fSig2); 
	
	if(print) cout << "pTrack1: "<< endl; 
	if(print) pl1->fPlab.Print(); 
	if(print) cout << "pTrack2: "<< endl; 
	if(print) pl2->fPlab.Print(); 

	
	TVector3 pM1;
	pM1.SetXYZ(pl1->fPlab.X(),pl1->fPlab.Y(),pl1->fPlab.Z());
	TVector3 pM2;
	pM2.SetXYZ(pl2->fPlab.X(),pl2->fPlab.Y(),pl2->fPlab.Z());
	
	
	fMuon1Vect=pM1; 
	fMuon2Vect=pM2; 
	
	fMuon1Pt = fMuon1Vect.Perp(); fMuon2Pt = fMuon2Vect.Perp();
	fMuon1Eta = fMuon1Vect.Eta(); fMuon2Eta = fMuon2Vect.Eta();
	fMuon1Phi = fMuon1Vect.Phi(); fMuon2Phi = fMuon2Vect.Phi();
	
	
	TLorentzVector fP1;
	TLorentzVector fP2;
	
	fP1.SetVectM(fMuon1Vect,MMUON);
	fP2.SetVectM(fMuon2Vect,MMUON);
	
	TLorentzVector fPUps;
	fPUps=fP1+fP2;
	
	fCandY = TMath::Abs(fPUps.Rapidity()); // Rapidity can be +/-, so set the magnitude
	fCand4V=fPUps;
	fCandMass=fCand4V.M();
	fCandPt=fCand4V.Perp();
	
	
	double dRmin=0.1;
	
	for(int i=0; i<fired.size(); i++){
		TTrgObj *pTrig = fpEvt->getTrgObj(fired[i]); 
		if(print) pTrig->dump(); 
		TVector3 fP;
		fP.SetXYZ(pTrig->fP.X(),pTrig->fP.Y(),pTrig->fP.Z()); 
		
		double dR1=pM1.DeltaR(fP);
		double dEta1=TMath::Abs(pM1.Eta()-fP.Eta()); 
		double dPhi1=pM1.DeltaPhi(fP);
		
		double dR2=pM2.DeltaR(fP);
		double dEta2=TMath::Abs(pM2.Eta()-fP.Eta()); 
		double dPhi2=pM2.DeltaPhi(fP);

		
		if(print) cout << "Trig: " << fired[i] <<  " DR(1): " << dR1 << endl; 
		if(print) cout <<  "Trig: " << fired[i] << " DR(2): " << dR2 << endl; 
		
		if(dR1<dR2 && dR1<dRmin) {
			index1.push_back(true);  
			index2.push_back(false);
		}
		if(dR2<dR1 && dR2<dRmin){
			index2.push_back(true); 
			index1.push_back(false); 
		}
		
	}
	
	int RV=0; 
	int nTrg=0; 
	for(int i=0; i<index1.size(); i++){
		
		if(print)cout << "Trigger[" << i << "] Muon1: " << index1[i] << endl; 
		if(print)cout << "Trigger[" << i << "] Muon2: " << index2[i] << endl; 
		
		if(index1[i]==true) RV=1; 
		if(index2[i]==true) RV=2; 
		if(index1[i]==true|| index2[i]==true) nTrg++; 
	}
	
	if(nTrg>1) RV=3; 
	
	if(print)cout << "Return Value: "<< RV << endl;

	return RV;
	
}

void xsReader::listHLT(){

	fCandY = fCandPt = fCandMass = -1.; 
	fpCand = 0;
	fGenCandY = fGenCandPt = -1.;
	fgCand = 0;
	fMuon1Eta = fMuon1Pt = fMuon2Eta = fMuon2Pt = -1.;
	fMuon1Vect.SetXYZ(-99.,-99.,-99.);
	fMuon1Vect.SetXYZ(-99.,-99.,-99.);
	
	int pt_Tmax=100; //singule muon triggers with a pT threshold of up to 100
	int pt_Tmin=10; //but >10
	
	double pTmin=-9; 
	
	int dimuon_path=get_path();
	bool PS=false; 
	bool muon_qual=false; 
	bool muon_qual1=false; 
	bool muon_qual2=false; 
	double vProb=-1; 
		
	bool Tmu=false; 

	for (int a=0; a<NHLT; ++a) {
		if(fpEvt->fHLTResult[a]==1){
			TString path=fpEvt->fHLTNames[a];
			for(int i=pt_Tmax; i>=pt_Tmin; i--){
				if(path.Contains(Form("HLT_Mu%d_v",i))){
					Tmu=true;
					pTmin=static_cast<double>(i); 
				}//it path
				
			}//loop over triggers
					
		}//check fired
	}//loop over triggers
	int Trig_muon=0; 
	
	if(fpEvt->nMuons()!=2 ) return; 
	if(Tmu==false) return; 

	if(Tmu && fpEvt->nMuons()==2) Trig_muon=SingleMu_Match(static_cast<int>(pTmin));
	select_dimuons(PS, muon_qual1,muon_qual2,vProb); 
	
	if(Tmu) ((TH1D*)fpHistFile->FindObjectAny("fired_triggers"))->Fill(Form("HLT_Mu%f",pTmin), 1); 
	if(Tmu) ((TH1I*)fpHistFile->FindObjectAny("nMuons"))->Fill(fpEvt->nMuons()); 
		
	
	//if(Tmu && AC)cout << "Single Mu: " << Tmu << " dimuon trigger: " << dimuon_path << " muon qual: " << muon_qual << " Analysis cuts: " << AC << endl; 
	
	//if(Tmu) cout << "Vertex Probability: " << vProb << " PS: "<< PS << endl; 
	//if(Tmu && vProb>0) cout << "pT threshold: " << pTmin << endl; 
	
	pTmin-=pTmin*0.1; 
	
	if(fpEvt->nMuons()==2 && Tmu && Trig_muon!=0) {
		((TH1D*)fpHistFile->FindObjectAny("Mass_SingleMu_Trig"))->Fill(fCandMass); 
		((TH1D*)fpHistFile->FindObjectAny("dipT_SingleMu_Trig"))->Fill(fCandPt); 
		((TH1D*)fpHistFile->FindObjectAny("pT_SingleMu_Trig"))->Fill(fMuon1Vect.Perp()); 
		((TH1D*)fpHistFile->FindObjectAny("pT_SingleMu_Trig"))->Fill(fMuon2Vect.Perp()); 
		//cout << "Trigger Muon: "<< Trig_muon << endl; 
	}
	
	if(Trig_muon==1) muon_qual=muon_qual1; 
	if(Trig_muon==2) muon_qual=muon_qual2; 
	
	if(Trig_muon==3){	
		TRandom3 r(0);
		int PL=r.Integer(2); 
		if(PL==0) muon_qual=muon_qual1; 
		else muon_qual=muon_qual2; 
	}
		
	
	if(fpEvt->nMuons()==2 && Tmu && dimuon_path>0 && PS && muon_qual1 && muon_qual2 && vProb>0.01) {
		//cout << "Pt: " << fCandPt << endl; 
		double pt1=fMuon1Vect.Perp();
		double pt2=fMuon2Vect.Perp();
		double pt_notfired=0; 
		string name="INIT"; 
		if(Trig_muon==1) {
			name=Form("gEff_nominal_AETA%d",eta_to_bin(fMuon2Vect.Eta()) );
			pt_notfired=fMuon2Vect.Perp();
		}
		if(Trig_muon==2){
			name=Form("gEff_nominal_AETA%d",eta_to_bin(fMuon1Vect.Eta()) );
			pt_notfired=fMuon1Vect.Perp();
		}
		
		if(Trig_muon==3){	
				TRandom3 r(0);
				int PL=r.Integer(2); 
			if(PL==0){
				pt_notfired=fMuon1Vect.Perp();
				name=Form("gEff_nominal_AETA%d",eta_to_bin(fMuon1Vect.Eta()) );
			}
			else{
				pt_notfired=fMuon2Vect.Perp();
				name=Form("gEff_nominal_AETA%d",eta_to_bin(fMuon2Vect.Eta()) );
			}
			/*
			cout << "pT threshold: " << pTmin << endl; 
			cout << "name: " << name << endl; 
			cout << "pt1: " << pt1 << endl;
			cout << "pt2: " << pt2 << endl; 
			cout <<  "pT fired: " << pt_fired << endl; 
			cout << endl; 
			 */
			
		}//both muons fire the trigger, just randomly pick one. 		
		
		TGraph *grmu = (TGraph*)mueff_file->FindObjectAny(name.c_str()); 
		
		double W=grmu->Eval(pt_notfired);
		delete grmu;
		//cout << "Weight: " << W << endl; 
		if(fCandY<0.6) ((TH1D*)fpHistFile->FindObjectAny("rho_numerator"))->Fill(fCandPt,1./W); 
		if(fCandY<0.6) ((TH1D*)fpHistFile->FindObjectAny("rho_numerator_unW"))->Fill(fCandPt); 

	}
	if(fpEvt->nMuons()==2 && Tmu && PS && muon_qual && vProb>0.01){
		//cout << "Pt: " << fCandPt << endl; 
		if(fCandY<0.6) ((TH1D*)fpHistFile->FindObjectAny("rho_denominator"))->Fill(fCandPt); 
	}
	
	if(fpEvt->nMuons()==2 && Tmu && PS && muon_qual && vProb>=0 && fCandY<0.6){
		((TH2D*)fpHistFile->FindObjectAny("vertex_prob_pt"))->Fill(fCandPt,vProb); 
		((TH2D*)fpHistFile->FindObjectAny("vertex_prob_m"))->Fill(fCandMass,vProb); 
	}
	if(fpEvt->nMuons()==2 && Tmu && PS && dimuon_path>0 && muon_qual && vProb>=0 && fCandY<0.6){
		((TH2D*)fpHistFile->FindObjectAny("vertex_prob_pt_trig"))->Fill(fCandPt,vProb); 
		((TH2D*)fpHistFile->FindObjectAny("vertex_prob_m_trig"))->Fill(fCandMass,vProb); 
	}

	
}

bool xsReader::isPathFired( TString Path ){
  bool HLT_Path = false; TTrgObj *pTrig(0);
  for (int a = 0; a < NHLT ; ++a) {
  if ( fpEvt->fHLTNames[a] ==  Path  && fpEvt->fHLTResult[a] == 1  ) {
      HLT_Path = true;
      //cout << Path << " fired!!!! "  << endl;
      for (int s = 0; s < fpEvt->nTrgObj() ; ++s) {
	pTrig = fpEvt->getTrgObj(s);
	//cout << pTrig->fLabel << endl;
      }
      
    }
  }
  
  return HLT_Path;
}
  
void xsReader::select_dimuons(bool &PS, bool &muon_qual1, bool &muon_qual2, double &vProb){
	fpCand=0; 
	TAnaCand *pCand(0);
	const int N=Cands.size(); 
	
	vector<double> vertex_prob; 
	vector<int> good_cand;

	
	for(int iC=0; iC <Cands.size(); iC++){
		pCand=Cands[iC]; 
		double mass=pCand->fMass;
		double ptCand = pCand->fPlab.Perp();
		double candeta = pCand->fPlab.Eta();
		
		TAnaTrack *pTra1(0); 
		TAnaTrack *pTra2(0);
		TAnaTrack *pl1(0);
		TAnaTrack *pl2(0);
		
	
		
		pl1 = fpEvt->getSigTrack(pCand->fSig1); 
		pl2 = fpEvt->getSigTrack(pCand->fSig2);
		pTra1 = fpEvt->getRecTrack(pl1->fIndex);
		pTra2 = fpEvt->getRecTrack(pl2->fIndex);  
    	
		int index1; 
		int index2; 
		
		for(int i=0; i<fpEvt->nMuons();i++){
			TAnaMuon *Muon = fpEvt->getMuon(i);
			if(Muon->fIndex==pl1->fIndex)
				index1=i; 
			if(Muon->fIndex==pl2->fIndex)
				index2=i;
		}//Get muon loop 
		
		TAnaMuon *pMuon1 = fpEvt->getMuon(index1);
		TAnaMuon *pMuon2 = fpEvt->getMuon(index2);
		
		TLorentzVector tmp1;
		TLorentzVector tmp2;
		
		tmp1.SetVectM(pl1->fPlab,MMUON);
		tmp2.SetVectM(pl2->fPlab,MMUON);
		
		TLorentzVector tmp_ups;
		tmp_ups=tmp1+tmp2;

		if(pl1->fQ*pl2->fQ>0)
			continue;
		
		if(!PS_cuts(tmp_ups,pl1,pl2,0))
			continue;
		if ( CowboyVeto(pCand) )
			continue; 
			
		vertex_prob.push_back(pCand->fVtx.fProb);
		good_cand.push_back(iC); 

	}//for iC loop 
	int best=-1; 
	
	if(good_cand.size()>0){
		const int N=good_cand.size(); 
		double prob_array[N]; 
		double probMax=0; 
	
			for (int iC=0; iC<good_cand.size(); iC++) {
				if (vertex_prob[iC] > probMax) {
					best = iC;
					probMax = vertex_prob[iC];
				}
				
		}
	}//good candidates, get the best
	
	//This function gets vectors (fMuon1Vect) etc for the best candidate 
	if(best<0) return; 
	
	//cout << "best candidate index: " << best << endl; 
	
	fpCand = fpEvt->getCand(best); 		
	vProb=fpCand->fVtx.fProb; 
	//cout << "fpCand: " << fpCand << endl;
	
	TAnaTrack *pl1 = fpEvt->getSigTrack(fpCand->fSig1); 
	TAnaTrack *pl2 = fpEvt->getSigTrack(fpCand->fSig2);
	
	TAnaTrack *pTrack1 = fpEvt->getRecTrack(pl1->fIndex);
	TAnaTrack *pTrack2 = fpEvt->getRecTrack(pl2->fIndex);
	
		
	int index1=-9; 
	int index2=-9;
	
	for(int i=0; i<fpEvt->nMuons();i++){
		TAnaMuon *Muon = fpEvt->getMuon(i);
		if(Muon->fIndex==pl1->fIndex)
			index1=i; 
		if(Muon->fIndex==pl2->fIndex)
			index2=i;
	}//Get muon loop 
	
	TAnaMuon *pMuon1 = fpEvt->getMuon(index1);
	TAnaMuon *pMuon2 = fpEvt->getMuon(index2);
	
	TLorentzVector fP1;
	TLorentzVector fP2;
	
	fP1.SetVectM(fMuon1Vect,MMUON);
	fP2.SetVectM(fMuon2Vect,MMUON);
	
	TLorentzVector fPUps;
	fPUps=fP1+fP2;
	
	if(PS_cuts(fPUps,pl1,pl2,0))
		PS=true; 
	
	bool cut_dxy=true; 
	bool cut_dz=true; 
	bool cut_nv_hits=true; 
	bool chi2_ndof=true; 
	bool bit4=true;
	bool bit12=true; 
	bool pixL=true;
	
	if( TMath::Abs(pTrack1->fBsLip) > DZ_MAX) cut_dz=false; 
	if( TMath::Abs(pTrack1->fBsTip)> DXY_MAX) cut_dxy=false; 
	if( pTrack1->fValidHits < NV_HITS_MIN) cut_nv_hits=false;
	int pixelLayersWithMeasurement1 = pixelLayersWithMeasurement1=numberOfPixLayers(pTrack1); 
	int pixelLayersWithMeasurement2 = pixelLayersWithMeasurement2=numberOfPixLayers(pTrack2); 
	if( pixelLayersWithMeasurement1<=1) pixL=false; 
	if((pTrack1->fChi2)/static_cast<double>(pTrack1->fDof)>MUON_CHI2_DOF_MAX) chi2_ndof=false; 
	if( (((pMuon1->fMuID)&(1<<MUTYPE1))>>MUTYPE1)!=1) bit4=false; 	
	if( (((pMuon1->fMuID)&(1<<MUTYPE2))>>MUTYPE2)!=1) bit12=false;
		
	if(cut_dz && cut_dxy && cut_nv_hits && chi2_ndof && bit4 && bit12 && pixL) 
		muon_qual1=true; 
	
	cut_dxy=true; 
	cut_dz=true; 
	cut_nv_hits=true; 
	chi2_ndof=true; 
	bit4=true;
	bit12=true; 
	pixL=true;
	
	if( TMath::Abs(pTrack2->fBsLip) > DZ_MAX) cut_dz=false; 
	if( TMath::Abs(pTrack2->fBsTip)> DXY_MAX) cut_dxy=false; 
	if( pTrack2->fValidHits < NV_HITS_MIN) cut_nv_hits=false;
	if( pixelLayersWithMeasurement2<=1) pixL=false; 
	if((pTrack2->fChi2)/static_cast<double>(pTrack2->fDof)>MUON_CHI2_DOF_MAX) chi2_ndof=false; 
	if( (((pMuon2->fMuID)&(1<<MUTYPE1))>>MUTYPE1)!=1) bit4=false; 	
	if( (((pMuon2->fMuID)&(1<<MUTYPE2))>>MUTYPE2)!=1) bit12=false;
	
	if(cut_dz && cut_dxy && cut_nv_hits && chi2_ndof && bit4 && bit12 && pixL) 
		muon_qual2=true; 
	
	
}

bool xsReader::isPathFired_Match( TString Path, TString Label ){
	bool print=0; 
	if(print) cout << "Cand size at isPathFired_Match " << Cands.size() << endl; 
	//cout << "fpEvt->nCands() " << fpEvt->nCands() << endl; 
  bool HLT_Path = false;
  bool Leg1 = false;
  bool Leg2 = false;
  TAnaCand *pCand(0);
  TAnaCand *pCand_(0);
  double rapidity, pt;
  bool Match = false;
 
  for (int iC = 0; iC < Cands.size() ; ++iC) {
    pCand = Cands[iC];
    TLorentzVector Candi;
    Candi.SetPtEtaPhiM(pCand->fPlab.Perp(),pCand->fPlab.Eta(),pCand->fPlab.Phi(),pCand->fMass);
    if(skim_analysis==1){
    if ( Candi.Rapidity() >= 0 ) ((TH2D*)fpHistFile->FindObjectAny(Form("TrigCheck_before_%dS",UPSTYPE)))->Fill(Candi.Rapidity(), pCand->fPlab.Perp());
    if ( Candi.Rapidity() < 0 ) ((TH2D*)fpHistFile->FindObjectAny(Form("TrigCheck_before_%dS",UPSTYPE)))->Fill(-Candi.Rapidity(), pCand->fPlab.Perp());
    }
}
	
  for (int a = 0; a < NHLT ; ++a) {
    if ( fpEvt->fHLTNames[a] ==  Path  && fpEvt->fHLTResult[a] == 1  ) {
      HLT_Path = true;
	 ((TH1D*)fpHistFile->FindObjectAny("Labels"))->Fill(fpEvt->fHLTNames[a],1); 
      //cout << Path << " fired!!!! "  << endl;
    }
  }
  
  if ( HLT_Path  ){
	  for (int iC = 0; iC < Cands.size() ; ++iC) {
		  pCand_ = Cands[iC];
		  TTrgObj *pTrig(0); TTrgObj *pTrig_(0); 
		  int t(-1);
		  TLorentzVector tagD;
		  TAnaTrack *pTagD(0);
		  pTagD = fpEvt->getSigTrack(pCand_->fSig1);
		  tagD.SetPtEtaPhiM(pTagD->fPlab.Pt(), pTagD->fPlab.Eta(), pTagD->fPlab.Phi(), MMUON);
		  for (int s = 0; s < fpEvt->nTrgObj() ; ++s) {
			  pTrig = fpEvt->getTrgObj(s);
			  if ( !(Label.CompareTo(pTrig->fLabel)) ) {
				  double tagD_dR = tagD.DeltaR(pTrig->fP);
				  double tagD_dEta = TMath::Abs(pTagD->fPlab.Eta() - pTrig->fP.Eta());
				  // double tagD_dPhi = TMath::Abs(pTagD->fPlab.Phi() - pTrig->fP.Phi());
				  double tagD_dPhi = TMath::Abs(pTagD->fPlab.DeltaPhi(pTrig->fP.Vect()));
				  //cout << "Tag Phi: " << pTagD->fPlab.Phi() << " Trig Phi: " << pTrig->fP.Phi() << endl; 
				  //cout << "tagD_dphi: " << tagD_dPhi << " " << TMath::Abs(pTagD->fPlab.DeltaPhi(pTrig->fP.Vect())) << endl; 		  
				  if ( ( tagD_dPhi < DPHI ) && ( tagD_dEta < DETA )) {
					  Leg1 = true;				
					  //cout << " Leg1 matched to Double mu T.O.  " << endl;
					  t=s;
					  break;
				  } 
			  }//trigger label check
		  }//loop over trigger objects
		  
		  TLorentzVector probe;
		  TAnaTrack *pProbe(0);
		  pProbe = fpEvt->getSigTrack(pCand_->fSig2);
		  probe.SetPtEtaPhiM(pProbe->fPlab.Pt(), pProbe->fPlab.Eta(), pProbe->fPlab.Phi(), MMUON);
		  for (int i = 0; i < fpEvt->nTrgObj() ; ++i) {
			  if ( i == t ) continue;
			  pTrig_ = fpEvt->getTrgObj(i);
			  
			  if ( !(Label.CompareTo(pTrig_->fLabel)) ) {
				  double probe_dR = probe.DeltaR(pTrig_->fP);
				  double probe_dEta = TMath::Abs(pProbe->fPlab.Eta() - pTrig_->fP.Eta());
				  //double probe_dPhi = TMath::Abs(pProbe->fPlab.Phi() - pTrig_->fP.Phi());
				  double probe_dPhi = TMath::Abs(pProbe->fPlab.DeltaPhi(pTrig_->fP.Vect()));
				  if ( ( probe_dPhi < DPHI ) && ( probe_dEta < DETA )) {
					  Leg2 = true;
					  //cout << " Leg2 matched to Double mu T.O. " << endl;
					  break;
				  }//check dPhi
			  }//Trigger label check
		  }//loop over trigger objects
		  if(print) cout << "Leg1: " << Leg1 << " Leg2: " << Leg2 << endl; 
		  if(Leg1 && Leg2){
			  Cands_TM.push_back(pCand_);//set vector of candidates to be used later 
		  }//Leg1 && 2
	  }//loop over candidates 
  }//HLT_Path true 
	
	if ( Leg1 && Leg2 ){
		Match = true;
		TLorentzVector Cand;
		Cand.SetPtEtaPhiM(pCand_->fPlab.Perp(),pCand_->fPlab.Eta(),pCand_->fPlab.Phi(),pCand_->fMass);
		if(skim_analysis==1){
			if ( Cand.Rapidity() >= 0 ) ((TH2D*)fpHistFile->FindObjectAny(Form("TrigCheck_after_%dS",UPSTYPE)))->Fill(Cand.Rapidity(), pCand_->fPlab.Perp() );
			if ( Cand.Rapidity() < 0 ) ((TH2D*)fpHistFile->FindObjectAny(Form("TrigCheck_after_%dS",UPSTYPE)))->Fill(-Cand.Rapidity(), pCand_->fPlab.Perp() );
		}//fill histograms
		//cout << " Match " << endl;
	}//If Leg1 & Leg2 is satisfied for at least on candidate 
    
	return Match;
}

void xsReader::fillCandVectors(){
	//these are some general checks. 	
	
	TAnaCand *pCand(0);
	for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
		pCand = fpEvt->getCand(iC);
		Cands.push_back(pCand);
	}//for loop 
}						

void xsReader::cutvariable_histograms(TAnaCand *pCand, TAnaTrack *pTrack1, TAnaTrack *pTrack2, TAnaMuon *pMuon1, TAnaMuon *pMuon2){
	
	fpHistFile->cd("various_histograms"); 
	/*
	cout << "vertex position: " << endl;  
	pCand->fVtx.fPoint.Print(); 
	cout << endl; 
	
	cout << pTrack1->fdxy << " " << pTrack1->fTip << endl;
	cout << pTrack1->fdz << " " << pTrack1->fLip << endl; 
	
	cout << pTrack2->fdxy << " " << pTrack2->fTip << endl; 
	cout << pTrack2->fdz << " " << pTrack2->fLip << endl;
	*/
	 
    ((TH1D*)fpHistFile->FindObjectAny("vertex_prob"))->Fill(pCand->fVtx.fProb);
    ((TH1D*)fpHistFile->FindObjectAny("vertex_significance"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);    
	
    ((TH1D*)fpHistFile->FindObjectAny("vertexchi2_dof"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
	
	((TH1D*)fpHistFile->FindObjectAny("dxy"))->Fill(pTrack1->fdxy);
	((TH1D*)fpHistFile->FindObjectAny("dz"))->Fill(pTrack1->fdz);
	((TH1D*)fpHistFile->FindObjectAny("dxy"))->Fill(pTrack2->fdxy);
	((TH1D*)fpHistFile->FindObjectAny("dz"))->Fill(pTrack2->fdz);
	
	((TH1D*)fpHistFile->FindObjectAny("fTip"))->Fill(pTrack1->fTip);
	((TH1D*)fpHistFile->FindObjectAny("fBsTip"))->Fill(pTrack1->fBsTip);
	((TH1D*)fpHistFile->FindObjectAny("fLip"))->Fill(pTrack1->fLip);
	((TH1D*)fpHistFile->FindObjectAny("fBsLip"))->Fill(pTrack1->fBsLip);

	((TH1D*)fpHistFile->FindObjectAny("fTip"))->Fill(pTrack2->fTip);
	((TH1D*)fpHistFile->FindObjectAny("fBsTip"))->Fill(pTrack2->fBsTip);
	((TH1D*)fpHistFile->FindObjectAny("fLip"))->Fill(pTrack2->fLip);
	((TH1D*)fpHistFile->FindObjectAny("fBsLip"))->Fill(pTrack2->fBsLip);


	((TH1D*)fpHistFile->FindObjectAny("trackchi2_dof"))->Fill(pTrack1->fChi2/static_cast<double>(pTrack1->fDof));
	((TH1D*)fpHistFile->FindObjectAny("muonchi2_dof"))->Fill(pMuon1->fMuonChi2/static_cast<double>(pMuon1->fTimeNdof));
	((TH1D*)fpHistFile->FindObjectAny("trackchi2_dof"))->Fill(pTrack2->fChi2/static_cast<double>(pTrack2->fDof));
	((TH1D*)fpHistFile->FindObjectAny("muonchi2_dof"))->Fill(pMuon2->fMuonChi2/static_cast<double>(pMuon2->fTimeNdof));
	
	((TH1D*)fpHistFile->FindObjectAny("numberOfValidHits"))->Fill(pTrack1->fValidHits);
	((TH1D*)fpHistFile->FindObjectAny("numberOfValidHits"))->Fill(pTrack2->fValidHits);
	
	((TH1D*)fpHistFile->FindObjectAny("numberOfPixelLayersWithMeasurement"))->Fill(numberOfPixLayers(pTrack1));
	((TH1D*)fpHistFile->FindObjectAny("numberOfPixelLayersWithMeasurement"))->Fill(numberOfPixLayers(pTrack2));

	
	((TH1D*)fpHistFile->FindObjectAny("bit4"))->Fill((((pMuon1->fMuID)&(1<<4))>>4));
	((TH1D*)fpHistFile->FindObjectAny("bit12"))->Fill((((pMuon1->fMuID)&(1<<12))>>12));					      
	((TH1D*)fpHistFile->FindObjectAny("bit4"))->Fill((((pMuon2->fMuID)&(1<<4))>>4));
	((TH1D*)fpHistFile->FindObjectAny("bit12"))->Fill((((pMuon2->fMuID)&(1<<12))>>12));
		
	((TH2D*)fpHistFile->FindObjectAny("doca_prob"))->Fill(pCand->fVtx.fProb,pCand->fMaxDoca);
	((TH1D*)fpHistFile->FindObjectAny("max_doca"))->Fill(pCand->fMaxDoca); 
	
	fpHistFile->cd(); 
	
}

void xsReader::cutvariable_histograms_MC(TAnaTrack *pTrack1, TAnaTrack *pTrack2){
  //        cout << "Cut histos: " << endl; 
  //	cout << "dxy: " << pTrack1->fdxy << " " << pTrack2->fdxy << endl; 
  //	cout << "dz: "  << pTrack1->fdz << " " << pTrack2->fdz << endl;
	fpHistFile->cd("various_histograms"); 
	((TH1D*)fpHistFile->FindObjectAny("dxy"))->Fill(pTrack1->fdxy);
	((TH1D*)fpHistFile->FindObjectAny("dz"))->Fill(pTrack1->fdz);
	((TH1D*)fpHistFile->FindObjectAny("dxy"))->Fill(pTrack2->fdxy);
	((TH1D*)fpHistFile->FindObjectAny("dz"))->Fill(pTrack2->fdz);
	
	((TH1D*)fpHistFile->FindObjectAny("trackchi2_dof"))->Fill(pTrack1->fChi2/static_cast<double>(pTrack1->fDof));
	((TH1D*)fpHistFile->FindObjectAny("trackchi2_dof"))->Fill(pTrack2->fChi2/static_cast<double>(pTrack2->fDof));
	
	((TH1D*)fpHistFile->FindObjectAny("numberOfValidHits"))->Fill(pTrack1->fValidHits);
	((TH1D*)fpHistFile->FindObjectAny("numberOfValidHits"))->Fill(pTrack2->fValidHits);
	
	fpHistFile->cd();
	
}

int  xsReader::numberOfPixLayers(TAnaTrack *pTrack) {
	bool layer1 = false, layer2 = false, layer3 = false, disk1=false, disk2=false;
	int hits = pTrack->fValidHits;
	//cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
	if(hits>20) hits=20; // pattern has only 20 locations
	for(int i =0; i<hits; ++i){
		unsigned int pat = pTrack->fHitPattern[i];
		//cout<<pat<<" " << endl;
		if( pat == 0x488 ) layer1 = true;
		else if( pat == 0x490 ) layer2 = true;
		else if( pat == 0x498 ) layer3 = true;
		else if( pat == 0x508 ) disk1 = true;
		else if( pat == 0x510 ) disk2 = true;
	}
	//cout<<dec<<endl;
	//cout << "L1: " << layer1 << " L2: " << layer2 << " L3: " << layer3 << " D1: " << disk1 << " D2: "<< disk2 << endl; 
	int pixHits=0;
	if(layer1) {pixHits++;}
	if(layer2) {pixHits++;}
	if(layer3) {pixHits++;}
	if(disk1) {pixHits++;}
	if(disk2) {pixHits++;}
	
	return pixHits;
}

bool xsReader::cuts(TAnaCand *pCand, TAnaTrack *pTra1, TAnaTrack *pTra2, TAnaTrack *pl1, TAnaTrack *pl2, TAnaMuon *pMuon1, TAnaMuon *pMuon2){
	//these are some general checks. 	
   TLorentzVector fPUps;
   TLorentzVector fP1;
   TLorentzVector fP2;
   fP1.SetVectM(pl1->fPlab, MMUON);
   fP2.SetVectM(pl2->fPlab, MMUON);
   fPUps=fP1+fP2;
   double ptCand = fPUps.Perp();
   double candeta = fPUps.Eta();
   double mass=fPUps.M();

  if(MODE==4){
		
	  //bit 10 == 1 [Track]
	  //bit 9/8/7 must be = 1 for barrel 
	  //bit 6-3 must sum to decimal number 1-3 for barrel 
	  //bit 1-0 must be 0 for a valid hit 
	  
	  int pixelLayersWithMeasurement1=0; 
	  pixelLayersWithMeasurement1=numberOfPixLayers(pTra1); 
	  int pixelLayersWithMeasurement2=0; 
	  pixelLayersWithMeasurement2=numberOfPixLayers(pTra2); 
	  
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Total Events",1);	  
		  
	  //Various muon cuts 
		if ( CowboyVeto(pCand) ){
		  return 0;}
	  	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Cowboy Veto",1);
	  
	  if ( (pl1->fPlab.Perp() > PTHI) || (pl1->fPlab.Perp() < PTLO) ) return 0;
	  if ( (pl2->fPlab.Perp() > PTHI) || (pl2->fPlab.Perp() < PTLO) ) return 0;
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon %.1f<P_{T}<%.1fGeV",PTLO,PTHI),1);

	  
	  if(ETAHI>0){
		  if ( TMath::Abs(pl1->fPlab.Eta()) > ETAHI || TMath::Abs(pl2->fPlab.Eta())>ETAHI)  {
			  ((TH1D*)fpHistFile->FindObjectAny("meta"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("meta2"))->Fill(mass,ptCand);
			  return 0;}
	  }
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon |#eta|<%.1f",ETAHI),1); 
	  
	  
	  ///////Muon in the barrel 
	  if ( ((TMath::Abs(pl1->fPlab.Eta()) <= ETABARREL) && (pl1->fPlab.Perp() < PTBARREL)) || ((TMath::Abs(pl2->fPlab.Eta()) <= ETABARREL) && (pl2->fPlab.Perp() < PTBARREL)) ){
		  ((TH1D*)fpHistFile->FindObjectAny("mPt_region1"))->Fill(mass);
		  ((TH2D*)fpHistFile->FindObjectAny("mPt_region1_2"))->Fill(mass,ptCand);
		  ((TH1D*)fpHistFile->FindObjectAny("muonPt_region1"))->Fill(pl1->fPlab.Perp());
		  ((TH1D*)fpHistFile->FindObjectAny("muonPt_region1"))->Fill(pl2->fPlab.Perp());
		  
		 // cout << "eta< " << ETABARREL << " and PT< " << PTBARREL<< endl; 
		 // cout << "pt1: " << pl1->fPlab.Perp() << " eta1: " << pl1->fPlab.Eta() << endl;
		 // cout << "pt2: " << pl2->fPlab.Perp() << " eta2: " << pl2->fPlab.Eta() << endl;
		  
		  return 0; 
	  }
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon |#eta|<%.1f, P_{T}<%.1fGeV",ETABARREL,PTBARREL), 1); 
	  
	  if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETABARREL) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAMI) && (pl1->fPlab.Perp() < PTMI))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETABARREL) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAMI) && (pl2->fPlab.Perp() < PTMI)) ){
		  ((TH1D*)fpHistFile->FindObjectAny("mPt_region2"))->Fill(mass);
		  ((TH2D*)fpHistFile->FindObjectAny("mPt_region2_2"))->Fill(mass,ptCand);
		  ((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl1->fPlab.Perp());
		  ((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl2->fPlab.Perp());

		  //cout << "Eta> " << ETABARREL << "and Eta< " << ETAMI << " and PT< " << PTMI<< endl; 
		  //cout << "pt1: " << pl1->fPlab.Perp() << " eta1: " << pl1->fPlab.Eta() << endl;
		  //cout << "pt2: " << pl2->fPlab.Perp() << " eta2: " << pl2->fPlab.Eta() << endl;
		  
		  return 0; 
	  }
	   ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon 1 or 2 %.1f<|#eta|<%.1f, P_{T}>%.1fGeV",ETABARREL,ETAMI,PTMI), 1); 
	  
	  
	  
	  if ( ((TMath::Abs(pl1->fPlab.Eta()) >= ETAMI) && ((TMath::Abs(pl1->fPlab.Eta()) <= ETAHI) && (pl1->fPlab.Perp() < PTLO))) || ((TMath::Abs(pl2->fPlab.Eta()) >= ETAMI) && (TMath::Abs(pl2->fPlab.Eta()) <= ETAHI) && (pl2->fPlab.Perp() < PTLO)) ){
		  ((TH1D*)fpHistFile->FindObjectAny("mPt_region2"))->Fill(mass);
		  ((TH2D*)fpHistFile->FindObjectAny("mPt_region2_2"))->Fill(mass,ptCand);
		  ((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl1->fPlab.Perp());
		  ((TH1D*)fpHistFile->FindObjectAny("muonPt_region2"))->Fill(pl2->fPlab.Perp());
		  
		 // cout << "Eta> " << ETAMI << "and Eta< " << ETAHI << " and PT< " << PTLO<< endl; 
		 // cout << "pt1: " << pl1->fPlab.Perp() << " eta1: " << pl1->fPlab.Eta() << endl;
		 // cout << "pt2: " << pl2->fPlab.Perp() << " eta2: " << pl2->fPlab.Eta() << endl;
		  
		  return 0; 
	  }
	  
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon 1 or 2 %.1f<|#eta|<%.1f, P_{T}>%.1fGeV",ETAMI,ETAHI,PTLO), 1); 

	  if(pixelLayersWithMeasurement1<=1) return 0; 
	  if(pixelLayersWithMeasurement2<=1) return 0; 
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("pixelLayersWithMeasurement>1",1); 
	  
	  //Various candidate kinematic cuts 
	    if ( fPUps.Rapidity() < -RAPCAND ){ 
			return 0;}
	  
	  if ( fPUps.Rapidity() > RAPCAND ){
		  return 0;}
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Cand |y|<%.1f",RAPCAND),1);
	  
	  if ( fPUps.Perp() < PTCAND_MIN ){ 
			return 0;}
		  
	    if ( fPUps.Perp() > PTCAND ){
		  return 0;}
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("%.1f<CandP_{T}<%.1fGeV",PTCAND_MIN,PTCAND),1);
	  
	double pt1=pl1->fPlab.Perp();
	double pt2=pl2->fPlab.Perp();	  
		  
	/// Vertex Cuts
	//fProb0.01	
	if(pCand->fVtx.fProb<VPROB_MIN){
		((TH1D*)fpHistFile->FindObjectAny("mProb"))->Fill(mass);
		((TH2D*)fpHistFile->FindObjectAny("mProb2"))->Fill(mass,ptCand);
		((TH1D*)fpHistFile->FindObjectAny("muonPtProb"))->Fill(pl1->fPlab.Perp());
		((TH1D*)fpHistFile->FindObjectAny("muonPtProb"))->Fill(pl2->fPlab.Perp());
			  return 0;}
	  ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Vert. Prob.>%.2f", VPROB_MIN),1);
    
	// Remove this cut for optimization
	//>vertexsig>3	  
	  if(VERTEX_SIG>0){  
		  if(pCand->fVtx.fD3d/pCand->fVtx.fD3dE > VERTEX_SIG){
			  ((TH1D*)fpHistFile->FindObjectAny("mVertexSig"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("mVertexSig2"))->Fill(mass,ptCand);
			  ((TH1D*)fpHistFile->FindObjectAny("muonPtVertexSig"))->Fill(pl1->fPlab.Perp());
			  ((TH1D*)fpHistFile->FindObjectAny("muonPtVertexSig"))->Fill(pl2->fPlab.Perp());
			  return 0;}
	  }
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("|fD3d/fD3dE|<%.1f",VERTEX_SIG),1);	  
		  
    //////////
    //fdz>15
	if( TMath::Abs(pTra1->fBsLip) > DZ_MAX || TMath::Abs(pTra2->fBsLip) > DZ_MAX){ 
			  ((TH1D*)fpHistFile->FindObjectAny("mdz"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("mdz2"))->Fill(mass,ptCand);
			  return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("|dz|<%.1fcm",DZ_MAX),1); 	  
	
	//dxy>0.3  
	if( TMath::Abs(pTra1->fBsTip)> DXY_MAX ||TMath::Abs(pTra2->fBsTip) > DXY_MAX){
			  
			  ((TH1D*)fpHistFile->FindObjectAny("mdxy"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("mdxy2"))->Fill(mass,ptCand);
			  
			  return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("|dxy|>%.1fcm",DXY_MAX*10),1);	  
	//Number of hits <10	  
	if( pTra1->fValidHits < NV_HITS_MIN || pTra2->fValidHits < NV_HITS_MIN){ 
			  ((TH1D*)fpHistFile->FindObjectAny("mHits"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("mHits2"))->Fill(mass,ptCand);
			  ((TH1D*)fpHistFile->FindObjectAny("muonPtHits"))->Fill(pl1->fPlab.Perp());
			  ((TH1D*)fpHistFile->FindObjectAny("muonPtHits"))->Fill(pl2->fPlab.Perp());
			  return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Valid Hits>%d",NV_HITS_MIN),1);
    
    //////////
    
	 //muon chi2/dof>1.8  
	
	if((pTra1->fChi2)/static_cast<double>(pTra1->fDof)>MUON_CHI2_DOF_MAX || (pTra2->fChi2)/static_cast<double>(pTra2->fDof)>MUON_CHI2_DOF_MAX){
			  ((TH1D*)fpHistFile->FindObjectAny("mChi"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("mChi2"))->Fill(mass,ptCand);
			  ((TH1D*)fpHistFile->FindObjectAny("muonPtChi"))->Fill(pl1->fPlab.Perp());
			  ((TH1D*)fpHistFile->FindObjectAny("muonPtChi"))->Fill(pl2->fPlab.Perp());
			  return 0;}
		  
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Muon #chi^{2}/DOF >%.1f",MUON_CHI2_DOF_MAX),1); 

	  
	  if( (((pMuon1->fMuID)&(1<<MUTYPE1))>>MUTYPE1)!=1 || (((pMuon1->fMuID)&(1<<MUTYPE2))>>MUTYPE2)!=1) {
		  ((TH1D*)fpHistFile->FindObjectAny("mbits"))->Fill(mass);
		  ((TH2D*)fpHistFile->FindObjectAny("mbits2"))->Fill(mass,ptCand);
		  ((TH1D*)fpHistFile->FindObjectAny("muonPtbits"))->Fill(pl1->fPlab.Perp());
		  ((TH1D*)fpHistFile->FindObjectAny("muonPtbits"))->Fill(pl2->fPlab.Perp());
		  return 0;}
	  


	  if( (((pMuon2->fMuID)&(1<<MUTYPE1))>>MUTYPE1)!=1 || (((pMuon2->fMuID)&(1<<MUTYPE2))>>MUTYPE2)!=1) {
          ((TH1D*)fpHistFile->FindObjectAny("mbits"))->Fill(mass);
          ((TH2D*)fpHistFile->FindObjectAny("mbits2"))->Fill(mass,ptCand);
		  ((TH1D*)fpHistFile->FindObjectAny("muonPtbits"))->Fill(pl1->fPlab.Perp());
		  ((TH1D*)fpHistFile->FindObjectAny("muonPtbits"))->Fill(pl2->fPlab.Perp());
          return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Tight Muons",1);  
    
    //////////
    	if ( pl1->fQ*pl2->fQ > 0 ){
			  ((TH1D*)fpHistFile->FindObjectAny("mfQ"))->Fill(mass);
			  ((TH2D*)fpHistFile->FindObjectAny("mfQ"))->Fill(mass,ptCand);
			  return 0;}
	((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill("Muons have opposite chage",1);
		
    if (mass < MASSLO) return 0;
    if (mass > MASSHI) return 0;
    ((TH1D*)fpHistFile->FindObjectAny("cuts_sequence"))->Fill(Form("Mass Range: %.1f<M_{#mu#mu}<%.1fGeV",MASSLO,MASSHI),1);
   
	//cout << "pixelLayersWithMeasurement: " << pixelLayersWithMeasurement << endl;

  }//Mode 4

  return true; 

}

void xsReader::candidateSelection(int mode){

    bool print=0; 

  fCandY = fCandPt = fCandMass = -1.; 
  fpCand = 0;
  fGenCandY = fGenCandPt = -1.;
  fgCand = 0;
  fMuon1Eta = fMuon1Pt = fMuon2Eta = fMuon2Pt = -1.;
  fMuon1Vect.SetXYZ(-99.,-99.,-99.);
  fMuon1Vect.SetXYZ(-99.,-99.,-99.);
  TAnaCand *pCand(0);
  vector<int> lcands; //list of candidates that pass the cuts

  int index1=0; 
  int index2=0; 

  if(print) cout << "Candidate Selection N Cand: " << Cands_TM.size() << endl; 
//put in mass ranges


  for (int iC = 0; iC < Cands_TM.size() ; ++iC) {
	  
	  pCand = Cands_TM[iC];
	  double mass=pCand->fMass;
	  double ptCand = pCand->fPlab.Perp();
	  double candeta = pCand->fPlab.Eta();

	  if(iC==0){ 
		  //fill for first Candidate mass
		  if(mass>9.26 && mass<=9.66) ((TH1I*)fpHistFile->FindObjectAny("Candidates_pre_1S"))->Fill(Cands_TM.size());
		  if(mass>9.92 && mass<=10.12) ((TH1I*)fpHistFile->FindObjectAny("Candidates_pre_2S"))->Fill(Cands_TM.size());
		  if(mass>10.255 && mass<=10.455) ((TH1I*)fpHistFile->FindObjectAny("Candidates_pre_3S"))->Fill(Cands_TM.size());
	  }
	  
    TAnaTrack *pTra1(0); 
    TAnaTrack *pTra2(0);
    TAnaTrack *pl1(0);
    TAnaTrack *pl2(0);
   
    pl1 = fpEvt->getSigTrack(pCand->fSig1); 
    pl2 = fpEvt->getSigTrack(pCand->fSig2);
    pTra1 = fpEvt->getRecTrack(pl1->fIndex);
    pTra2 = fpEvt->getRecTrack(pl2->fIndex);  
    	
	fpHistFile->cd("various_histograms"); 
    ((TH2D*)fpHistFile->FindObjectAny("dxy_candphi"))->Fill(pTra1->fdxy,pCand->fPlab.Phi());
    ((TH2D*)fpHistFile->FindObjectAny("dxy_candphi"))->Fill(pTra2->fdxy,pCand->fPlab.Phi());
	fpHistFile->cd();
	  
	fpHistFile->cd("candidate_selection"); 
	  ((TH1D*)fpHistFile->FindObjectAny("vertex_prob_all"))->Fill(pCand->fVtx.fProb);
	fpHistFile->cd(); 
	  
    for(int i=0; i<fpEvt->nMuons();i++){
      TAnaMuon *Muon = fpEvt->getMuon(i);
      if(Muon->fIndex==pl1->fIndex)
        index1=i; 
      if(Muon->fIndex==pl2->fIndex)
        index2=i;
    }//Get muon loop 

   
    TAnaMuon *pMuon1 = fpEvt->getMuon(index1);
    TAnaMuon *pMuon2 = fpEvt->getMuon(index2);
    
    /////Make all the cut variable histograms of all candidates 
   
    cutvariable_histograms(pCand,pTra1, pTra2, pMuon1,pMuon2); 
   
    if(print){
      cout << "pTra1 fPvIdx: " << pTra1->fPvIdx << " pTra2 fPvIdx " << pTra2->fPvIdx << endl; 
      cout << "pTra1 fdxy: " << pTra1->fdxy << " pTra2: fdxy " << pTra2->fdxy << endl; 
      cout << "pTra1 fTip: " << pTra1->fTip << " pTra2: fTip " << pTra2->fTip << endl; 
      cout << "pTra1 fBsTip: " <<pTra1->fBsTip << " pTra2: fBsTip " << pTra2->fBsTip << endl; 
      cout << "pTra1 fLip: " << pTra1->fLip << " pTra2: fLip " << pTra2->fLip << endl; 
      cout << "pTra1 fBsLip: " << pTra1->fBsLip << " pTra2: fBsLip " << pTra2->fBsLip << endl; 
      cout << endl; 
    }
	
    ((TH2D*)fpHistFile->FindObjectAny("CandVxy"))->Fill(pCand->fVtx.fPoint.X(),pCand->fVtx.fPoint.Y());
    ((TH1D*)fpHistFile->FindObjectAny("CandVz"))->Fill(pCand->fVtx.fPoint.Z());

 if (TYPE != pCand->fType) continue;
   
	//remove for seagull study 
  if( cuts(pCand,pTra1,pTra2,pl1,pl2,pMuon1,pMuon2)==0){
			if(print) cout << "Cuts[" << iC << "] Failed. "<< endl; 
           continue;}
  else if(print)
		 cout << "Cuts[" << iC <<"] Passed." << endl; 
	/*
	  
	  //section for imposing PS cuts only 
	  TLorentzVector pMu1_cand; 
	  pMu1_cand.SetVectM(pl1->fPlab,MMUON);
	  
	  TLorentzVector pMu2_cand; 
	  pMu2_cand.SetVectM(pl2->fPlab,MMUON);
	  
	  TLorentzVector pUps_cand; 
	  pUps_cand=pMu1_cand + pMu2_cand; 
	  
	  if(PS_cuts(pUps_cand, pl1,pl2,0)==0){
		  continue; 
	   
	  }
	  */ 
	  
   //Fill optimization Histograms here. 

fpHistFile->cd("vertex_optimization");
 TLorentzVector p1;
 TLorentzVector p2;
 p1.SetVectM(pl1->fPlab,MMUON);
 p2.SetVectM(pl2->fPlab,MMUON);
 TLorentzVector v; 
 v=p1+p2;
 double candY = v.Rapidity();

 fpHistFile->cd();  

  lcands.push_back(iC); // candidates that pass the cuts 	 
	  
  }//iC Cand loop 
  
  
  if (0 == lcands.size()) return; // if no events pass the cuts, just return 
  int best(-1);
  if (lcands.size() > 1) {
    if(print) cout << "MORE THAN ONE CANDIDATE  " << lcands.size()   <<  endl;
    double ptMax(0.), pt(0.);
    double maxDocaMax(99.), maxDoca(0.);
    double chi2Max(99.), chi2(0.);
	double probMax(0), prob(0);
    
	  
	for (unsigned int iC = 0; iC < lcands.size(); ++iC) {
      //loop over all the candidates that were good and find 
	  pCand = fpEvt->getCand(lcands[iC]);       
      if ( mode == 1 ){
		  //max pt - probably not a good mode to run in this mode
	   pt = pCand->fPlab.Perp(); 
	   if (pt > ptMax) {
	   best = lcands[iC]; 
	   ptMax = pt;
		}
	}// mode 1
	else if ( mode == 2 ){
		//Get the event with the minimum distance of closest approach 
	maxDoca = pCand->fMaxDoca;
	if (maxDoca < maxDocaMax) {
	  best = lcands[iC]; 
	  maxDocaMax = maxDoca;
	  }
      }//else if mode 2 
	else if ( mode == 3 ){
		//Get the event with the lowest vertex Chi2
	chi2 = pCand->fVtx.fChi2;
	  if (chi2 < chi2Max) {
	    best = lcands[iC]; 
	    chi2Max = chi2;
	  }
	  }// else if mode 3	
	else if ( mode == 4 ){
		//Get the event with the highest vertex probability 
		prob = pCand->fVtx.fProb;
		if (prob > probMax) {
			best = lcands[iC]; 
			probMax = prob;
		}
		if(print) cout << "prob[" << iC << "]: " << prob << endl; 
	}// else if mode 4		
		
	}//iC loop over good candidates 
	  if(print) cout << "best: " << best << " max prob: " << probMax << endl; 
  }//if statement for dealing with more than one candidate 
  
  else if (lcands.size() == 1) {
	  //if there is just one candidate anyway, then just get the first one 
	  best = lcands[0];}
	
	if(print) cout << "Getting vectors[" << best << "]" << endl; 
	
	
	const int NCand=lcands.size(); 
	double vProb[NCand]; 
	int index[NCand]; 
	
	for (int iC=0;iC<NCand; iC++) {
		pCand = fpEvt->getCand(lcands[iC]);       
		vProb[iC]=pCand->fVtx.fProb; 
	}
	
	TMath::Sort(NCand, vProb,index); 
	if(NCand>1){
		double M1=fpEvt->getCand(lcands[index[0]])->fMass; 
		double M2=fpEvt->getCand(lcands[index[1]])->fMass; 
		((TH2D*)fpHistFile->FindObjectAny("m1_m2"))->Fill(M1,M2); 
	}
	if(NCand>0){
		double M=fpEvt->getCand(0)->fMass;
		((TH1I*)fpHistFile->FindObjectAny("Candidates"))->Fill(NCand);
		if(M>9.26 && M<=9.66) ((TH1I*)fpHistFile->FindObjectAny("Candidates_1S"))->Fill(NCand);
		if(M>9.92 && M<=10.12) ((TH1I*)fpHistFile->FindObjectAny("Candidates_2S"))->Fill(NCand);
		if(M>10.255 && M<=10.455) ((TH1I*)fpHistFile->FindObjectAny("Candidates_3S"))->Fill(NCand); 
	}
	GetDataVectors(best);
	
}

void xsReader::GetDataVectors(int best){
	//This function gets vectors (fMuon1Vect) etc for the best candidate 
	
	bool print=0; 
	if(print) cout << "best: " << best << endl; 
	if(best<=-1) return;
	
		fpCand = Cands_TM[best];
		fpHistFile->cd("candidate_selection");
		((TH1D*)fpHistFile->FindObjectAny("vertex_prob_best"))->Fill(fpCand->fVtx.fProb);
		fpHistFile->cd();						
		
		TAnaTrack *pl1 = fpEvt->getSigTrack(fpCand->fSig1); 
		TAnaTrack *pl2 = fpEvt->getSigTrack(fpCand->fSig2);
				
		TAnaTrack *pTrack1 = fpEvt->getRecTrack(pl1->fIndex);
		TAnaTrack *pTrack2 = fpEvt->getRecTrack(pl2->fIndex);
		
		//use signal tracks for momentum stuff
		if(pl1->fQ==1){
			//always make muon 1 the positive muon
			fMuon1Vect=pl1->fPlab;
			fMuon2Vect=pl2->fPlab;
		}
		else{
			fMuon1Vect=pl2->fPlab;
			fMuon2Vect=pl1->fPlab;
		}
		
		fMuon1Pt = fMuon1Vect.Perp(); fMuon2Pt = fMuon2Vect.Perp();
		fMuon1Eta = fMuon1Vect.Eta(); fMuon2Eta = fMuon2Vect.Eta();
		fMuon1Phi = fMuon1Vect.Phi(); fMuon2Phi = fMuon2Vect.Phi();
		
		TLorentzVector fP1;
		TLorentzVector fP2;
		
		fP1.SetVectM(fMuon1Vect,MMUON);
		fP2.SetVectM(fMuon2Vect,MMUON);
		
		TLorentzVector fPUps;
		fPUps=fP1+fP2;
		
		
		fCandY = TMath::Abs(fPUps.Rapidity()); // Rapidity can be +/-, so set the magnitude
		fCand4V=fPUps;
		fCandMass=fCand4V.M();
		fCandPt=fCand4V.Perp();
		
		fP1.SetVectM(fMuon1Vect,MMUON);
		fP2.SetVectM(fMuon2Vect,MMUON);
		
		fPUps=fP1+fP2;
		
		fMuonPhi12 =pl1->fQ*(pl1->fPlab.DeltaPhi(pl2->fPlab));
		
		fMuon1PtE = pTrack1->fPtE;    fMuon2PtE = pTrack2->fPtE;
		fMuon1EtaE = pTrack1->fEtaE;  fMuon2EtaE = pTrack2->fEtaE;
		fMuon1PhiE = pTrack1->fPhiE;  fMuon2PhiE = pTrack2->fPhiE;
		
		if(print) cout << endl; 
}

void xsReader::fillCandHist() {
    
	bool print =0; 
	
	if(print){
		cout << "xsReader::fillCandHist " << endl; 
		
		cout << "Dimuon Mass: " << fCandMass << " DiMuon Pt: " << fCandPt << " Dimuon y: " << fCandY << " Weight: " << fWeight << endl;  
		cout << "Muon1 Pt: " << fMuon1Pt << " Eta: " << fMuon1Eta << endl; 
		cout << "Muon2 Pt: " << fMuon2Pt << " Eta: " << fMuon2Eta << endl; 
		cout << "Number of Tracks: " << fpEvt->nRecTracks() << endl;
	}
    
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEta"))->Fill(fMuon1Eta);
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEta"))->Fill(fMuon2Eta);
    ((TH1D*)fpHistFile->FindObjectAny("SigMuPt"))->Fill(fMuon1Pt);
    ((TH1D*)fpHistFile->FindObjectAny("SigMuPt"))->Fill(fMuon2Pt);  
  
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEtaPt"))->Fill(fMuon1Eta,fMuon1Pt);  
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEtaPt"))->Fill(fMuon2Eta,fMuon2Pt);  
	//mode ==1

	  
    ((TH1D*)fpHistFile->FindObjectAny("weights"))->Fill(fWeight);
    ((TH1I*)fpHistFile->FindObjectAny("NTracks"))->Fill(fpEvt->nRecTracks());
    ((TH1D*)fpHistFile->FindObjectAny("CandPhi"))->Fill(fCand4V.Phi(),fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("CandMass"))->Fill(fCandMass,fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("CandPt"))->Fill(fCandPt,fWeight);
    ((TH2D*)fpHistFile->FindObjectAny("CandPxPy"))->Fill(fCand4V.Px(),fCand4V.Py(),fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("CandRapidity"))->Fill(fCandY,fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("CandEta"))->Fill(fCand4V.Eta(),fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("UpsilonMass"))->Fill(fCandMass,fWeight);
    
    if(fCandMass>9.26 && fCandMass<9.66 && fCandY<fYbin2[1])
      ((TH1D*)fpHistFile->FindObjectAny("Y1Pt"))->Fill(fCandPt,fWeight);

    if(fCandMass>9.92 && fCandMass<10.12 && fCandY<fYbin2[1])
      ((TH1D*)fpHistFile->FindObjectAny("Y2Pt"))->Fill(fCandPt,fWeight); 

    if(fCandMass>10.255 && fCandMass<10.455 && fCandY<fYbin2[1])
      ((TH1D*)fpHistFile->FindObjectAny("Y3Pt"))->Fill(fCandPt,fWeight); 
	
	if(fCandMass<9.2 && fCandY<fYbin2[1]) 
		((TH1D*)fpHistFile->FindObjectAny("SB_L_Pt"))->Fill(fCandPt,fWeight);
	if(fCandMass>10.6 && fCandY<fYbin2[1])
		((TH1D*)fpHistFile->FindObjectAny("SB_H_Pt"))->Fill(fCandPt,fWeight);


    ((TH1D*)fpHistFile->FindObjectAny("SigMuEta"))->Fill(fMuon1Eta,fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEta"))->Fill(fMuon2Eta,fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("SigMuPt"))->Fill(fMuon1Pt,fWeight);
    ((TH1D*)fpHistFile->FindObjectAny("SigMuPt"))->Fill(fMuon2Pt,fWeight);  
  
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEtaPt"))->Fill(fMuon1Eta,fMuon1Pt); // Not sure what these actually are supposed to be...  
    ((TH1D*)fpHistFile->FindObjectAny("SigMuEtaPt"))->Fill(fMuon2Eta,fMuon2Pt);  
    
	  if(fCandMass>9.26 && fCandMass<9.66){
		  ((TProfile*)fpHistFile->FindObjectAny("Y1mass_muonpt"))->Fill(fMuon1Pt,(fCandMass-9.4603)*1000); 
		  ((TProfile*)fpHistFile->FindObjectAny("Y1mass_muonpt"))->Fill(fMuon2Pt,(fCandMass-9.4603)*1000); 
	  }//mass with 100MeV of Y1 peak
	
	  
	for ( int ipt = 1; ipt < fNpt; ++ipt ){
		if(fMuon1Vect.Perp()>=fPTbin[ipt] && fMuon1Vect.Perp()<fPTbin[ipt+1]){
			((TH1D*)fpHistFile->FindObjectAny(Form("UpsilonMass_ptmu%.1f_%.1f", fPTbin[ipt], fPTbin[ipt+1])))->Fill(fCandMass,fWeight);
		}
		if(fMuon2Pt>=fPTbin[ipt] && fMuon2Pt<fPTbin[ipt+1]){
			((TH1D*)fpHistFile->FindObjectAny(Form("UpsilonMass_ptmu%.1f_%.1f", fPTbin[ipt], fPTbin[ipt+1])))->Fill(fCandMass,fWeight);
		}
	}//pt loop
	  

    for ( int iy = 0; iy < fNy2; ++iy ){
      for ( int ipt = 0; ipt < fNpt1; ++ipt ){
	    if(fCandPt>=fPTbin1[ipt] && fCandPt<fPTbin1[ipt+1] && fCandY>=fYbin2[iy] && fCandY<fYbin2[iy+1]){
	   ((TH1D*)fpHistFile->FindObjectAny(Form("UpsilonMass_y_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->Fill(fCandMass,1./dimuon_efficiency);
	   // cout << "Filling. " << endl; 
		}//pt & upsilon if statement 
	
      }//ipt loop
    }//iy loop 

  
} 

int xsReader::eta_to_bin(double eta){
	
	//binning for single muon efficiencies from AN2012-088
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

void xsReader::calculate_dimuon_efficiency(){
  double eff1(-99); double eff2(-99);
  double eff1_Ep(-99); double eff2_Ep; 
  double eff1_Em(-99); double eff2_Em; 
  
  double rho(-99); 
	
  double eff_dimuon(-99); 	
	
  TAnaCand *pCand;
  TLorentzVector Cand;
	
  bool print=0; 
  
  fpHistFile->cd();
  /*
	  //Run in mode 2. Not sure what mode 1 does. 
    TAnaTrack *pl1 = fpEvt->getSigTrack(fpCand->fSig1); 
    TAnaTrack *pl2 = fpEvt->getSigTrack(fpCand->fSig2);
    
    effID1 = fPidTable2011SeagullPos->effD(pl1->fPlab.Perp(), pl1->fPlab.Eta(), 0.);//there is only one bin in phi, so may as well set phi=0
    effID2 = fPidTable2011SeagullPos->effD(pl2->fPlab.Perp(), pl2->fPlab.Eta(), 0.);
	  
	effID1_Ep=effID1 + fPidTable2011SeagullPos->errD(pl1->fPlab.Perp(), pl1->fPlab.Eta(), 0.);
	effID2_Ep=effID2 + fPidTable2011SeagullPos->errD(pl2->fPlab.Perp(), pl2->fPlab.Eta(), 0.);
	  
	effID1_Em=effID1 - fPidTable2011SeagullNeg->errD(pl1->fPlab.Perp(), pl1->fPlab.Eta(), 0.);
	effID2_Em=effID2 - fPidTable2011SeagullNeg->errD(pl2->fPlab.Perp(), pl2->fPlab.Eta(), 0.);  
	*/  
	
	//Nominal efficiencies: 
	
	if(TMath::Abs(fMuon1Vect.Eta()) >1.6 || TMath::Abs(fMuon2Vect.Eta())>1.6 ){
		cout << "Single muon efficiencies comuputed to eta of 1.6 only." << endl; 
		dimuon_efficiency=1; 
		dimuon_efficiencyEp=1;
		dimuon_efficiencyEm=1;
		return;
	}
	
	string name1= Form("gEff_nominal_AETA%d",eta_to_bin(fMuon1Vect.Eta()) );
	string name2= Form("gEff_nominal_AETA%d",eta_to_bin(fMuon2Vect.Eta()) );

	//Sys +
	
	string name1P=Form("gEff_totsys_p_AETA%d",eta_to_bin(fMuon1Vect.Eta()));
	string name2P=Form("gEff_totsys_p_AETA%d",eta_to_bin(fMuon2Vect.Eta()));
	
	//Sys-
	
	string name1M=Form("gEff_totsys_m_AETA%d",eta_to_bin(fMuon1Vect.Eta()));
	string name2M=Form("gEff_totsys_m_AETA%d",eta_to_bin(fMuon2Vect.Eta()));

	
	TGraph *grmu1 = (TGraph*)mueff_file->FindObjectAny(name1.c_str()); 
	TGraph *grmu2 = (TGraph*)mueff_file->FindObjectAny(name2.c_str()); 
	
	TGraph *grmu1P = (TGraph*)mueff_file->FindObjectAny(name1P.c_str()); 
	TGraph *grmu2P = (TGraph*)mueff_file->FindObjectAny(name2P.c_str()); 
	
	TGraph *grmu1M = (TGraph*)mueff_file->FindObjectAny(name1M.c_str()); 
	TGraph *grmu2M = (TGraph*)mueff_file->FindObjectAny(name2M.c_str()); 
	
	TH1D *rho_pt = (TH1D*)rho_file->FindObjectAny("rho_pt"); 
	
	eff1=grmu1->Eval(fMuon1Vect.Perp()); 
	eff2=grmu2->Eval(fMuon2Vect.Perp()); 
	
	eff1_Ep=grmu1P->Eval(fMuon1Vect.Perp()); 
	eff2_Ep=grmu2P->Eval(fMuon2Vect.Perp()); 
	
	eff1_Em=grmu1M->Eval(fMuon1Vect.Perp()); 
	eff2_Em=grmu2M->Eval(fMuon2Vect.Perp()); 
	
	//rho=rho_pt->GetBinContent(rho_pt->FindBin(fCandPt)); 
	rho =rho_pt->Interpolate(fCandPt); 

	
	if(fCandPt<50) rho=1; // just use 1 for dimuon pt<50. 
	
    fWeight = 1;
    dimuon_efficiency = eff1*eff2*rho; // Two data driven muon efficiencies * rho factor gives total dimuon efficiency. 
	dimuon_efficiencyEp = eff1_Ep*eff2_Ep*rho; 
	dimuon_efficiencyEm = eff1_Em*eff2_Em*rho;
	
	  //cout << "pt1: " << pl1->fPlab.Perp() << " pt2: " << pl2->fPlab.Perp() << " eff: " << MuIdWeight << endl; 
	  if(eff1<0) cout << "eff1: " << eff1 << " pt1: " << fMuon1Vect.Perp() << " eta: " << fMuon1Vect.Eta() << " phi: " << fMuon1Vect.Phi() << endl; 
	  if(eff2<0) cout << "eff2: " << eff2 << " pt2: " << fMuon2Vect.Perp() << " eta: " << fMuon2Vect.Eta() << " phi: " << fMuon2Vect.Phi() << endl; 
	  
	  if(print) cout << "eff1: " << eff1 << " pt1: " << fMuon1Vect.Perp() << " eta: " << fMuon1Vect.Eta() << " phi: " << fMuon1Vect.Phi() << endl; 
	  if(print) cout << "eff2: " << eff2 << " pt2: " << fMuon2Vect.Perp() << " eta: " << fMuon2Vect.Eta() << " phi: " << fMuon2Vect.Phi() << endl; 
	  if(print) cout << "effmumu: " << dimuon_efficiency << " dimuon pt: " << fCandPt << endl; 
	  
	  if(dimuon_efficiency<0) cout << "Dimuon efficiency: " << dimuon_efficiency << " +/- " << dimuon_efficiencyEp << " " << dimuon_efficiencyEm << endl; 
	
	  if(eff1_Ep<0 || eff1_Em<0) cout << "less than 0." << endl; 
	  if(eff2_Ep<0 || eff2_Em<0) cout << "less than 0." << endl; 
	  
    fpHistFile->cd("efficiency");
	
   	if(TMath::Abs(fCandY)<fYbin2[1]){
		((TH1D*)fpHistFile->FindObjectAny("Eff_denominator"))->Fill(fCandPt,1); 
		((TH1D*)fpHistFile->FindObjectAny("Eff_numerator"))->Fill(fCandPt,dimuon_efficiency);   
		((TH1D*)fpHistFile->FindObjectAny("Eff_numerator_Ep"))->Fill(fCandPt,dimuon_efficiencyEp);   
		((TH1D*)fpHistFile->FindObjectAny("Eff_numerator_Em"))->Fill(fCandPt,dimuon_efficiencyEm);  
	}
	
    fpHistFile->cd();
  
	delete grmu1;
	delete grmu2;
	delete grmu1P;
	delete grmu2P;
	delete grmu1M;
	delete grmu2M;
	delete rho_pt; 
  
}

void xsReader::bookHist() {
  cout << "--> xsReader> bookHist> " << endl;
  fBin = BIN;
  fMassLow = MASSLO;
  fMassHigh = MASSHI;
	
	
  const double dm_max=1000;
  const int ndm_bins=static_cast<int>(dm_max/.25); 
  const int nm_bins=static_cast<int> ((fMassHigh-fMassLow)/0.01);
  
  TH1 *h;
  TH2 *k;
  TProfile2D *prof;
  TProfile *prof1D; 
		
	
	fpHistFile->cd("shifts");
	
	prof1D = new TProfile("Y1mass_muonpt","(M_{#mu#mu}-M_{PDG}) MeV vs Muon P_{T};Muon P_{T} GeV;Deviation from PDG value MeV",25,0,50); 
	
	h = new TH1D("Mass_Shift_low", "#DeltaM=M'-M, P_{T}<20GeV",100,-2,30);
	((TH1D*)fpHistFile->FindObjectAny("Mass_Shift_low"))->GetXaxis()->SetTitle("#DeltaM MeV");
	
	h = new TH1D("Pt_Shift_low","#DeltaP_{T}=P'_{T}-P_{T}",100,-100,100);
	((TH1D*)fpHistFile->FindObjectAny("Pt_Shift_low"))->GetXaxis()->SetTitle("#DeltaP_{T} MeV");

	 
	 h = new TH1D("Mass_Shift_high", "#DeltaM=M'-M, P_{T}>20GeV",100,-100,100);
	 ((TH1D*)fpHistFile->FindObjectAny("Mass_Shift_high"))->GetXaxis()->SetTitle("#DeltaM MeV");
	 
	 h = new TH1D("Pt_Shift_high","#DeltaP_{T}=P'_{T}-P_{T}",100,-100,100);
	 ((TH1D*)fpHistFile->FindObjectAny("Pt_Shift_low"))->GetXaxis()->SetTitle("#DeltaP_{T} MeV");
	  
	for(int ipt=1; ipt<fNpt; ipt++){
		
		h = new TH1D(Form("Mass_Shift_Pt_%.1f_%.1f",fPTbin[ipt],fPTbin[ipt+1]),Form("#DeltaM=M'-M, %.1f<P_{T}<%.1fGeV",fPTbin[ipt],fPTbin[ipt+1]),1000,-100,100);
		((TH1D*)fpHistFile->FindObjectAny(Form("Mass_Shift_Pt_%.1f_%.1f",fPTbin[ipt],fPTbin[ipt+1])))->GetXaxis()->SetTitle("#DeltaM MeV");
		h = new TH1D(Form("Pt_Shift_Pt_%.1f_%.1f",fPTbin[ipt], fPTbin[ipt+1]),Form("#DeltaP_{T}=P'_{T}-P_{T},  %.1f<P_{T}<%.1fGeV",fPTbin[ipt],fPTbin[ipt+1]),1000,-300,300);
		((TH1D*)fpHistFile->FindObjectAny(Form("Pt_Shift_Pt_%.1f_%.1f",fPTbin[ipt], fPTbin[ipt+1])))->GetXaxis()->SetTitle("#DeltaP_{T} MeV");
	}//for loop for muon Pt
		  
  fpHistFile->cd();
	
	// Histogram for paths 
	h = new TH1D("Paths","Paths",100,-1,12);//bin for version numbers
	//((TH1D*)fpHistFile->FindObjectAny("Paths"))->SetTitle("Paths");
	h = new TH1D("fired_triggers","Fired Triggers", 1000, 0, 1000); 
	h = new TH1I("nMuons", "Number of Muons", 20, 0, 20); 
	((TH1D*)fpHistFile->FindObjectAny("Paths"))->SetStats(0);
	((TH1D*)fpHistFile->FindObjectAny("Paths"))->SetBit(TH1::kCanRebin);
	
	h = new TH1D("Mass_SingleMu_Trig", "Mass_SinlgeMu_Trig", 2000,0.5,12.5); 
	h = new TH1D("dipT_SingleMu_Trig", "p_{T}(#mu#mu) for Single Muon Triggered Events", 100,0,100); 
	h = new TH1D("pT_SingleMu_Trig", "Single #mu p_{T}", 100, 0,100); 
	
	h = new TH1D("rho_numerator","Fired #mu#mu trigger & #mu trigger", fNpt1, fPTbin1); 
	h = new TH1D("rho_numerator_unW","Fired #mu#mu trigger & #mu trigger un-Weighted", fNpt1, fPTbin1); 

	h = new TH1D("rho_denominator","Fired #mu trigger", fNpt1, fPTbin1); 
	
	h = new TH1D("seagull_numerator","Fired #mu#mu trigger & seagull", fNpt1, fPTbin1); 
	h = new TH1D("seagull_denominator","Fired #mu#mu trigger", fNpt1, fPTbin1); 
	
	h = new TH1D("seagull_numerator_SB","Fired #mu#mu trigger & seagull Sideband Regions only", fNpt1, fPTbin1); 
	h = new TH1D("seagull_denominator_SB","Fired #mu#mu trigger Sideband regions only", fNpt1, fPTbin1); 
	
	//Histogram for labels
	h = new TH1D("Labels","Labels",50,0,3);//bin for labels
	// ((TH1D*)fpHistFile->FindObjectAny("Labels"))->SetTitle("Labels");
	((TH1D*)fpHistFile->FindObjectAny("Labels"))->SetStats(0);
	((TH1D*)fpHistFile->FindObjectAny("Labels"))->SetBit(TH1::kCanRebin);	
	

  h = new TH1D("weights", "Histogram of Weights", 1000,0,0.25);	
  h = new TH1I("DataSets","Run Numbers in File",5,0,5); 
  h = new TH1D("nPV", "Number of Primary Vertices in good events", 30,0,30); 	
  h = new TH1I("Candidates", "Candidates meeting cuts per Event", 100, 0,10);	
  h = new TH1I("Candidates_1S", "Candidates meeting cuts per Event", 100, 0,10);	
  h = new TH1I("Candidates_2S", "Candidates meeting cuts per Event", 100, 0,10);	
  h = new TH1I("Candidates_3S", "Candidates meeting cuts per Event", 100, 0,10);
	
  h = new TH1I("Candidates_pre_1S","Candidates per Event pre cut",100,0,10);
  h = new TH1I("Candidates_pre_2S","Candidates per Event pre cut",100,0,10);	
  h = new TH1I("Candidates_pre_3S","Candidates per Event pre cut",100,0,10);
	
  h = new TH1I("NTracks","Reco Tracks per Event; Number of reconstructed tracks per event", 100,0, 1000);
  h = new TH1I("MCTracks", "Number of Reco Tracks in MC", 100,0,5);
	const int  NPTbins=100;
	const int  PT0=0; 
	const int PTmax=100;
	
	//mass histograms of cut events 
	
	fpHistFile->cd("cuts");
	
	k = new TH2D("UpsilonMass2","UpsilonMass-PT",BIN,fMassLow,fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("UpsilonMass2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("UpsilonMass2"))->GetYaxis()->SetTitle("Pt,GeV");
	
	k = new TH2D("mProb2","Probability cut-PT",BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mProb2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mProb2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mVertexSig2","Vertex Distance Selection-PT",BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mVertexSig2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mVertexSig2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mIndex2","cut on index numbers",BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mIndex2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mIndex2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mdz2","cut on dz", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mdz2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mdz2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mdxy2","cut on dxy", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mdxy2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mdxy2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mHits2","cuts on hits",BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mHits2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mHits2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mChi2","cut on Chi2",BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mChi2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mChi2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mbits2", "cut on bits", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mbits2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mbits2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mPt_region1_2", "cut on Pt Region1", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mPt_region1_2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mPt_region1_2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mPt_region2_2", "cut on Pt Region2", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mPt_region2_2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mPt_region2_2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("meta2", "cut on eta", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("meta2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("meta2"))->GetYaxis()->SetTitle("Pt,GeV");
	k = new TH2D("mfQ2", "cut on fQ", BIN, fMassLow, fMassHigh,NPTbins,PT0,PTmax);
	((TH2D*)fpHistFile->FindObjectAny("mfQ2"))->GetXaxis()->SetTitle("Mass,GeV");
	((TH2D*)fpHistFile->FindObjectAny("mfQ2"))->GetYaxis()->SetTitle("Pt,GeV");

	fpHistFile->cd("candidate_selection");
	h =new TH1D("vertex_prob_all","Vertex Probability of all Candidates",100,0,1);
	h = new TH1D("vertex_prob_best","Vertex Probability of the Best Candidate",100,0,1); 
	
	k = new TH2D("m1_m2", "M_{#mu#mu}^{1} vs M_{#mu#mu}^{2} for the two best candidates", nm_bins, MASSLO, MASSHI, nm_bins, MASSLO,MASSHI); 
	
	fpHistFile->cd();
	
	fpHistFile->cd("dm_hist");

	
	prof1D = new TProfile("dm_prof","Profile of dm vs M_{#mu#mu} ;M_{#mu#mu} GeV;dm MeV",fNpt2,fPTbin2); 
	
	//bin 1
	
	fpHistFile->cd("dm_hist1"); 
	for ( int iy = 0; iy < fNy2; ++iy ){
		for ( int ipt = 0; ipt < fNpt1; ++ipt ){
			k = new TH2D(Form("m_dm_%dS_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
						 Form("m_dm_%dS,rapidity_bin1_%.1f_%.1f,pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
						 nm_bins,fMassLow, fMassHigh,ndm_bins,0,dm_max);  
			((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetXaxis()->SetTitle("Mass,GeV");
			((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetYaxis()->SetTitle("dm,MeV");
		}//pt loop
	}//y loop
	
	if(MODE!=5){
		//only create efficiency weightedhistograms if we aren't running MC
		//same histogram with efficiency weighting 
		for ( int iy = 0; iy < fNy2; ++iy ){
			for ( int ipt = 0; ipt < fNpt1; ++ipt ){
				k = new TH2D(Form("m_dm_%dS_effw_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
							 Form("m_dm_%dS,bin1 Efficiency weighted rapidity_%.1f_%.1f,pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin[ipt], fPTbin1[ipt+1]),
							 nm_bins,fMassLow, fMassHigh,ndm_bins,0,dm_max);  
				((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effw_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetXaxis()->SetTitle("Mass,GeV");
				((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effw_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetYaxis()->SetTitle("dm,MeV");
			}//pt loop
		}//y loop
		
		//same histogram with efficiency weighting, sys+
		for ( int iy = 0; iy < fNy2; ++iy ){
			for ( int ipt = 0; ipt < fNpt1; ++ipt ){
				k = new TH2D(Form("m_dm_%dS_effwEp_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
							 Form("m_dm_%dS,bin1 Efficiency, sys+ weighted rapidity_%.1f_%.1f,pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
							 nm_bins,fMassLow, fMassHigh,ndm_bins,0,dm_max);  
				((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effwEp_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetXaxis()->SetTitle("Mass,GeV");
				((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effwEp_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetYaxis()->SetTitle("dm,MeV");
			}//pt loop
		}//y loop
		
		//same histogram with efficiency weighting, sys-
		for ( int iy = 0; iy < fNy2; ++iy ){
			for ( int ipt = 0; ipt < fNpt1; ++ipt ){
				k = new TH2D(Form("m_dm_%dS_effwEm_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
							 Form("m_dm_%dS,bin1 Efficiency, sys- weighted rapidity_%.1f_%.1f,pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
							 nm_bins,fMassLow, fMassHigh,ndm_bins,0,dm_max);  
				((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effwEm_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetXaxis()->SetTitle("Mass,GeV");
				((TH1D*)fpHistFile->FindObjectAny(Form("m_dm_%dS_effwEm_bin1_rapidity_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetYaxis()->SetTitle("dm,MeV");
			}//pt loop
		}//y loop
	}//mode 5
	
	


	for(int ipt =0; ipt<fNpt1; ++ipt){
		
		k = new TH2D(Form("dm_y_ptbin1_%.1f_%.1f",fPTbin1[ipt],fPTbin1[ipt+1]),Form("dm vs y in P_{T} bins %.1f-%.1f;y;dm(MeV)", fPTbin1[ipt],fPTbin1[ipt+1]),100,0,2.5,100,0,dm_max);
		
	}
	
	fpHistFile->cd("dm_hist");
	

      for(int iy =0; iy<fNy2; ++iy){

		  k = new TH2D(Form("dm_pt_ybin%.1f_%.1f", fYbin2[iy],fYbin2[iy+1]), 
		      Form("dm vs pt in rapidity bins %.1f_%.1f", fYbin2[iy], fYbin2[iy+1]),fNpt1,fPTbin1, ndm_bins, 0,dm_max);
		  ((TH2D*)fpHistFile->FindObjectAny(Form("dm_pt_ybin%.1f_%.1f", fYbin2[iy],fYbin2[iy+1])))->GetXaxis()->SetTitle("P_{T},GeV");
		  ((TH2D*)fpHistFile->FindObjectAny(Form("dm_pt_ybin%.1f_%.1f", fYbin2[iy],fYbin2[iy+1])))->GetYaxis()->SetTitle("dm,MeV");
      }

	for(int ipt =0; ipt<fNpt2; ++ipt){
		
		k = new TH2D(Form("dm_y_ptbin_%.1f_%.1f",fPTbin2[ipt],fPTbin2[ipt+1]),Form("dm vs y in P_{T} bins %.1f-%.1f;y;dm(MeV)", fPTbin2[ipt],fPTbin2[ipt+1]),100,0,2.5,100,0,dm_max);
		
	}
	

	
	
	k = new TH2D("pt_dm","Pt and dm",NPTbins,PT0,PTmax,100,0,dm_max);
	((TH2D*)fpHistFile->FindObjectAny("pt_dm"))->GetXaxis()->SetTitle("Candidate P_{T}, GeV"); 
	((TH2D*)fpHistFile->FindObjectAny("pt_dm"))->GetYaxis()->SetTitle("dm,MeV"); 
	k = new TH2D("y_dm","y and dm",100,0,2.5,100,0,dm_max);
	((TH2D*)fpHistFile->FindObjectAny("y_dm"))->GetXaxis()->SetTitle("y(#mu#mu)"); 
	((TH2D*)fpHistFile->FindObjectAny("y_dm"))->GetYaxis()->SetTitle("dm,MeV");
	
	
	prof = new TProfile2D("pt_y_dm","Profile of dm as a function of P_{T} & y",100,0,100,25,0,1.4,0,dm_max);
	((TProfile2D*)fpHistFile->FindObjectAny("pt_y_dm"))->GetXaxis()->SetTitle("P_{T},GeV");
	((TProfile2D*)fpHistFile->FindObjectAny("pt_y_dm"))->GetYaxis()->SetTitle("y");
	
	k = new TH2D("deta_dm","Muon Delta Eta and dm", 200,-2.5,2.5,200,0,dm_max);
	k = new TH2D("phi_dm","Phi and dm",100,-3.2,3.2,200,0,dm_max);
	k = new TH2D("DR","DeltaR dm", 100,0,5,200,0,dm_max);
	

	fpHistFile->cd();

	int NPTEbins=100; 
	double PTEmin=0; 
	double PTEmax=200;

	int NEtaEbins=100;
	double NEtaEmin=0.;
	double NEtaEmax=5.;

	fpHistFile->cd("errors");

	k = new TH2D("muonptE_pt", "Muon PtE*1000 vs. Muon Pt", NPTbins, PT0,PTmax,NPTEbins,PTEmin,PTEmax);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_pt"))->GetXaxis()->SetTitle("Muon Pt, GeV");
	((TH2D*)fpHistFile->FindObjectAny("muonptE_pt"))->GetYaxis()->SetTitle("Muon PtE, MeV");
	k = new TH2D("muonptE_eta", " Muon PtE*1000 vs Muon Eta", 100,0,2.5,NPTEbins, PTEmin,PTEmax);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_eta"))->GetXaxis()->SetTitle("Muon Eta");
	((TH2D*)fpHistFile->FindObjectAny("muonptE_eta"))->GetYaxis()->SetTitle("Muon PtE, MeV");
	k = new TH2D("muonptE_dphi", "Muon dphi vs. MuonPtE*1000", 100,0,3.5, NPTEbins,PTEmin,PTEmax);
	((TH2D*)fpHistFile->FindObjectAny("muonptE_dphi"))->GetXaxis()->SetTitle("Muon dphi");
	((TH2D*)fpHistFile->FindObjectAny("muonptE_dphi"))->GetYaxis()->SetTitle("Muon PtE, MeV");


	k = new TH2D("muonetaE_pt", "Track #sigma_{#eta}*1000 vs Muon P_{T} GeV", NPTbins, PT0,PTmax,NEtaEbins,NEtaEmin,NEtaEmax);
	k = new TH2D("muonetaE_eta", "Track #sigma_{#eta}*1000 vs Muon #eta", 100,0,2.5,NEtaEbins,NEtaEmin,NEtaEmax);
	k = new TH2D("muonetaE_dphi", "Track #sigma_{#eta}*1000 vs Muon #Delta#phi", 100,0,2.5, NEtaEbins,NEtaEmin,NEtaEmax);

	k = new TH2D("muondphiE_pt", "Track #sigma_{#phi}*1000 vs Muon P_{T} GeV", NPTbins, PT0,PTmax,100,0,5.);
	k = new TH2D("muondphiE_eta", "Track #sigma_{#phi}*1000 vs Muon #eta", 100,0,2.5,100,0,5.);
	k = new TH2D("muondphiE_dphi", "Track #sigma_{#phi}*1000 vs Muon #Delta#phi", 100,0,3.5, 100,0,5.);

	h = new TH1D("ptE_pt","Muon #sigma_{P_{T}}/P_{T} [%]", 100,0,5);
	((TH1D*)fpHistFile->FindObjectAny("ptE_pt"))->GetXaxis()->SetTitle("Percent Error");

	h = new TH1D("ptE_pt_etabin0_0.8","#sigma_{P_{T}}/P_{T} [%] |#eta|<0.8",100,0,5);
	h = new TH1D("ptE_pt_etabin0.8_2.5","#sigma_{P_{T}}/P_{T} [%] 0.8<|#eta|<2.5",100,0,5);

        for(int iy =0; iy<fNy; ++iy){

         h = new TH1D(Form("ptE_pt_etabin%.1f_%.1f", fYbin[iy],fYbin[iy+1]), 
		      Form("#sigma_{P_{T}}/P_{T} [%] %.1f<|#eta|<%.1f", fYbin[iy], fYbin[iy+1]), 100,0,5);

                  }

	k = new TH2D("ptE_pt_eta", "|#eta| vs. #sigma_{P_{T}}/P_{T}", 100,0,2.5,100,0,5);
	((TH2D*)fpHistFile->FindObjectAny("ptE_pt_eta"))->GetYaxis()->SetTitle("Percent Error");
	((TH2D*)fpHistFile->FindObjectAny("ptE_pt_eta"))->GetXaxis()->SetTitle("|#eta|");
	h = new TH1D("etaE_eta","Muon EtaE*100/Eta", 100,0,5);
	((TH1D*)fpHistFile->FindObjectAny("etaE_eta"))->GetXaxis()->SetTitle("Percent Error");
	h = new TH1D("dphiE_phi", "Muon dphiE*100/dphi",100,0,5);
	((TH1D*)fpHistFile->FindObjectAny("dphiE_phi"))->GetXaxis()->SetTitle("Percent Error");
	k = new TH2D("sigma_pt","#sigma_{P_{T}(#mu#mu)} vs P_{T}^{gen} GeV; P_{T}^{gen} GeV; #sigma_{P_{T}(#mu#mu) GeV}",90,10,100, 1000,0,10); 


	fpHistFile->cd("kinematics");
	k = new TH2D("genMass_RecoMass","GenMass vs. Reco Mass",nm_bins,fMassLow,fMassHigh,nm_bins,fMassLow,fMassHigh);


	for ( int iy = 0; iy < fNy; ++iy ){
	  h = new TH1D(Form("muonPt_rapidity_%dS%.1f_%.1f", UPSTYPE, fYbin[iy], fYbin[iy+1]),Form("muonPt_rapidity_%dS%.1f_%.1f, 9.3<M_{#mu#mu}<10.3", UPSTYPE, fYbin[iy], fYbin[iy+1]),250\
		       ,0,100);

	}//y loop

	k = new TH2D("MuonPt_Candpt","Muon P_{T} vs. Cand P_{T}, 9.3<M_{#mu#mu}<10.3", 100, 0, 100, 100,0,100); 
	((TH2D*)fpHistFile->FindObjectAny("MuonPt_Candpt"))->GetXaxis()->SetTitle("Cand P_{T}"); 
	((TH2D*)fpHistFile->FindObjectAny("MuonPt_Candpt"))->GetYaxis()->SetTitle("Muon P_{T}"); 
	k = new TH2D("EventsY1_y_pt", "Number of events for Y1, 9.26<M_{#mu#mu}<9.66",100,0,100, 60, 0, 1.2); 
	k = new TH2D("EventsY2_y_pt", "Number of events for Y1, 9.92<M_{#mu#mu}<10.12",100,0,100, 60, 0, 1.2);
	k = new TH2D("EventsY3_y_pt", "Number of events for Y1, 10.255<M_{#mu#mu}<10.455",100,0,100, 60, 0, 1.2);

	k = new TH2D("muon_kinematics","Muon #eta vs P_{T}",100,0,100,250,-2.5,2.5);  
	k = new TH2D("muon_kinematics_y0.0_0.6", "Muon #eta vs P_{T} for y<0.6", 100,0,100, 250,-2.5,2.5); 

	h = new TH1D("opening_angle", "Opening Angle #theta_{12}",100,-3.14,3.14 );

	k = new TH2D("muon_pt_upsilon1S_pt","Muon Pt vs Upsilon 1S Pt, 9.26<M<9.66", 100,0,100, 50,0,50); 
	
	for (int iy=0; iy<fNy2; iy++) {
		k = new TH2D(Form("eta1_eta2_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]), Form("#eta_{1}(#mu) vs. #eta_{2}(#mu) %.1f<|y(#mu#mu)|<%.1f, All P_{T}(#mu#mu) bins;#eta_{1}(#mu);#eta_{2}(#mu)",fYbin2[iy],fYbin2[iy+1]), 100,0,2.5,100,0,2.5); 
		k = new TH2D(Form("eta1_eta2_aftercuts_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]), Form("#eta_{1}(#mu) vs. #eta_{2}(#mu) %.1f<|y(#mu#mu)|<%.1f After cuts, All P_{T}(#mu#mu) bins;#eta_{1}(#mu);#eta_{2}(#mu)",fYbin2[iy],fYbin2[iy+1]), 100,0,2.5,100,0,2.5); 

	}
	
	for (int iy=0; iy<fNy2; iy++) {
		for (int ipt=0; ipt<fNpt1; ipt++) {
			k = new TH2D(Form("eta1_eta2_bin1_y%.1f_%.1f_pt%.1f_%.1f",fYbin2[iy],fYbin2[iy+1],fPTbin1[ipt],fPTbin1[ipt+1]), Form("#eta_{1}(#mu) vs. #eta_{2}(#mu) %.1f<|y(#mu#mu)|<%.1f, %.1f<P_{T}(#mu#mu)<%.1f;#eta_{1}(#mu);#eta_{2}(#mu)",fYbin2[iy],fYbin2[iy+1],fPTbin1[ipt],fPTbin1[ipt+1]), 100,0,2.5,100,0,2.5);
		}
		
	}
	
	
	
	
	fpHistFile->cd();
	fpHistFile->cd("dm_hist");	
	
	for ( int ipt = 0; ipt < fNpt2; ++ipt ){
		k = new TH2D(Form("m_dm_pt%.1f_%.1f", fPTbin2[ipt], fPTbin2[ipt+1]), 
					 Form("m_dm_pt%.1f_%.1f", fPTbin2[ipt], fPTbin2[ipt+1]), 100, fMassLow, fMassHigh,100,0,250);
		//((TH1D*)fpHistFile->FindObjectAny(Form("Rapidity_IntegratedMass,pt%.1f_%.1f",  fPTbin[ipt], fPTbin[ipt+1])))->Sumw2(); 
		}
	fpHistFile->cd();

	//Histograms for muon chi2 and vertex selection	
	
	fpHistFile->cd("vertex_optimization");
	k = new TH2D("vertex_prob_m","Vertex Probability vs M_{#mu#mu}", nm_bins, fMassLow, fMassHigh, 1000,0,1); 
	k = new TH2D("vertex_prob_m_trig","Vertex Probability vs M_{#mu#mu}", nm_bins, fMassLow, fMassHigh, 1000,0,1); 
	
	k = new TH2D("vertex_prob_pt","Vertex Probability vs p_{T}", fNpt1,fPTbin1, 1000,0,1); 
	k = new TH2D("vertex_prob_pt_trig","Vertex Probability vs p_{T}", fNpt1, fPTbin1,1000,0,1); 
	
	
	
 fpHistFile->cd();
  

// space bins for vertex optimization every 0.05

  fpHistFile->cd("acceptance"); 	
	
	if(MODE==5){	
		//Only create and fill these histograms when running of MC
		// Acceptance Histograms
		k = new TH2D(Form("AllGenRes_%dS",  UPSTYPE), Form("AllGenRes_%dS", UPSTYPE), fNy, fYbin, fNpt, fPTbin);
		k = new TH2D(Form("RecoGenRes_%dS", UPSTYPE), Form("RecoGenRes_%dS", UPSTYPE), fNy, fYbin, fNpt, fPTbin); 
		((TH2D*)fpHistFile->FindObjectAny(Form("AllGenRes_%dS", UPSTYPE)))->Sumw2();
		((TH2D*)fpHistFile->FindObjectAny(Form("RecoGenRes_%dS", UPSTYPE)))->Sumw2();
		
		// Upsilon Gun Acceptance Histograms
		k = new TH2D(Form("UG_pt1_pt2_%dS",UPSTYPE), Form("UG_pt1_pt2_%dS, |y|<%.1f, %.1f<P_{T}<%.1f",UPSTYPE, fYbin2[1], fPTbin2[0],fPTbin2[1]), 100,0,100, 100,0,100); 
		k = new TH2D(Form("UG_eta1_eta2_%dS",UPSTYPE), Form("UG_eta1_eta2_%dS, |y|<%.1f, %.1f<P_{T}<%.1f",UPSTYPE, fYbin2[1], fPTbin2[0],fPTbin2[1]), 50,0,RAPCAND, 50,0,RAPCAND); 
		
		double AccpTlo=10; 
		double AccpThi=100.;
		double AccBW=1.0;
		int NBAcc=static_cast<int> ((AccpThi-AccpTlo)/AccBW);
		
		k = new TH2D(Form("UG_AllGenRes_y_pt_%dS",  UPSTYPE), Form("UG_AllGenRes_y_pt_%dS", UPSTYPE), fNy2, fYbin2, NBAcc, AccpTlo, AccpThi); 
		h = new TH1D(Form("UG_AllGenRes_%dS",UPSTYPE), Form("UG_AllGenRes_%dS, |y|<%.1f",UPSTYPE,fYbin2[1]),NBAcc, AccpTlo, AccpThi);
		h = new TH1D(Form("UG_AllGenRes_Ep_%dS",UPSTYPE), Form("UG_AllGenRes_%dS, |y|<%.1f, +1#sigma",UPSTYPE,fYbin2[1]),NBAcc, AccpTlo, AccpThi);
		h = new TH1D(Form("UG_AllGenRes_Em_%dS",UPSTYPE), Form("UG_AllGenRes_%dS, |y|<%.1f, -1#sigma",UPSTYPE,fYbin2[1]),NBAcc, AccpTlo, AccpThi);
		k = new TH2D(Form("UG_RecoGenRes_y_pt_%dS", UPSTYPE), Form("UG_RecoGenRes_%dS", UPSTYPE), fNy2, fYbin2, NBAcc, AccpTlo, AccpThi); 
		h = new TH1D(Form("UG_RecoGenRes_Ep_%dS",UPSTYPE),Form("UG_RecoGenRes_Ep_%dS,  |y|<%.1f, +1#sigma", UPSTYPE, fYbin2[1]),NBAcc, AccpTlo, AccpThi);
		h = new TH1D(Form("UG_RecoGenRes_Em_%dS",UPSTYPE),Form("UG_RecoGenRes_Em_%dS,  |y|<%.1f, -1#sigma", UPSTYPE, fYbin2[1]),NBAcc, AccpTlo, AccpThi);
		h = new TH1D(Form("UG_RecoGenRes_%dS",UPSTYPE),Form("UG_RecoGenRes_%dS,  |y|<%.1f", UPSTYPE, fYbin2[1]),NBAcc, AccpTlo, AccpThi);
		
		((TH2D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_y_pt_%dS", UPSTYPE)))->Sumw2();
		((TH1D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_%dS", UPSTYPE)))->Sumw2();
		((TH1D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_Ep_%dS", UPSTYPE)))->Sumw2();
		((TH1D*)fpHistFile->FindObjectAny(Form("UG_AllGenRes_Em_%dS", UPSTYPE)))->Sumw2();	
		((TH2D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_y_pt_%dS", UPSTYPE)))->Sumw2();
		((TH1D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_%dS", UPSTYPE)))->Sumw2();
		((TH1D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_Ep_%dS", UPSTYPE)))->Sumw2();
		((TH1D*)fpHistFile->FindObjectAny(Form("UG_RecoGenRes_Em_%dS", UPSTYPE)))->Sumw2();
		
		// Polarization Histo
		h = new TH1D("Ups_phi_weighted", "Ups #phi weighted", 100, -4, 4.);
		h = new TH1D("Ups_y_weighted", "Ups y weighted", 100, 0, 1.2);
		h = new TH1D("Ups_pt_weighted","Ups P_{T} weighted", 90,10,100); 	
		h = new TH1D("PolarizationWeight", "Polarization Weight", 30, 0, 1);
	}
	
  fpHistFile->cd();
  // Trig Check
  k = new TH2D(Form("TrigCheck_after_%dS",UPSTYPE), Form("TrigCheck_after_%1.dS", UPSTYPE), fNy, fYbin, fNpt, fPTbin);
  k = new TH2D(Form("TrigCheck_before_%dS",UPSTYPE), Form("TrigCheck_before_%1.dS", UPSTYPE), fNy, fYbin, fNpt, fPTbin);
  ((TH2D*)fpHistFile->FindObjectAny(Form("TrigCheck_after_%dS", UPSTYPE)))->Sumw2();
  ((TH2D*)fpHistFile->FindObjectAny(Form("TrigCheck_before_%dS", UPSTYPE)))->Sumw2();     

 
  fpHistFile->cd("efficiency");

  for (int iy=0; iy<fNy2; iy++) {
	  h = new TH1D(Form("seagull_num_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]),Form("Seagull Numerator, |y(#mu#mu)|<%.1f",fYbin2[iy+1]), 90,10,100); 
	  h = new TH1D(Form("seagull_den_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]),Form("Seagull Denominator, |y(#mu#mu)|<%.1f",fYbin2[iy+1]),90,10,100);
	  k = new TH2D(Form("vertex_prob_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]),Form("Vertex Probability Numerator, |y(#mu#mu)|<%.1f",fYbin2[iy+1]),fNpt1,fPTbin1,1000,0,1); 
	  k= new TH2D(Form("dxy_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]),Form("dxy_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]), fNpt1, fPTbin1, 100,0,.5);  
	  k= new TH2D(Form("dz_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]),Form("dz_eff_y%.1f_%.1f",fYbin2[iy],fYbin2[iy+1]), fNpt1, fPTbin1, 200,0,35);  


  }		

	
  //Total Dimuon Efficiency

  
  h = new TH1D("Eff_numerator", Form("Eff Numerator, |y(#mu#mu)|<%.1f",fYbin2[1]), fNpt2, fPTbin2); 
  h = new TH1D("Eff_denominator",Form("Eff Denominator, |y(#mu#mu)|<%.1f",fYbin2[1]), fNpt2,fPTbin2); 	
	
  h = new TH1D("Eff_numerator_Ep", Form("Eff Numerator, |y(#mu#mu)|<%.1f",fYbin2[1]), fNpt2, fPTbin2); 
 	
  h = new TH1D("Eff_numerator_Em", Form("Eff Numerator, |y(#mu#mu)|<%.1f",fYbin2[1]), fNpt2, fPTbin2); 
 	
  ((TH1D*)fpHistFile->FindObjectAny("Eff_numerator"))->Sumw2();  
  ((TH1D*)fpHistFile->FindObjectAny("Eff_denominator"))->Sumw2();	
	
  ((TH1D*)fpHistFile->FindObjectAny("Eff_numerator_Ep"))->Sumw2();  
 	
  ((TH1D*)fpHistFile->FindObjectAny("Eff_numerator_Em"))->Sumw2();  
 	
  h = new TH1D("Eff_vertexsig_num","Numerator for vertex f3D/#sigma_{f3D} efficiency",fNpt2,fPTbin2); 	
  h = new TH1D("Eff_vertexsig_den","Denominator for vertex f3D/#sigma_{f3D} efficiency",fNpt2,fPTbin2); 
	

	
	
  fpHistFile->cd();
  
  fpHistFile->cd("vertex_position"); 
  k = new TH2D("CandPxPy","Candidate Px-Py", 100,-100,100,100,-100,100);
  k = new TH2D("CandVxy", "Candidate Vertex x-y position", 100,-1,1,100,-1,1);
  h = new TH1D("CandVz", "Candidate Vertex z position", 100,-15,15);
  k = new TH2D("CandVPhidxy", "Candidate vertex phi vs dxy", 100,-3.4,3.4,100,-1,1);

  fpHistFile->cd("gen_info"); 	
	h = new TH1D("genUpsilonMass", "Gen Level Upsilon Mass", 250,fMassLow,fMassHigh);
	h = new TH1D("gen_reco_ptdiff", "#DeltaP_{T}=P_{T}(gen)=P_{T}(reco) MeV; #Delta P_{T} MeV",200,-100,100);  
	h = new TH1D("Egamma","Photon Energy",10000,0,10);
	k = new TH2D("GenVPhidxy", "Gen vertex phi vs. dxy", 100,-3.4,3.4, 100,-1,1);
	k = new TH2D("GenVxy", "Gen Vertex x-y position", 100,-1,1,100,-1,1);
	h = new TH1D("GenVz", "Gen Vertex z position", 100,-15,15);

  fpHistFile->cd("summary"); 
	h = new TH1D("UpsilonMass", "UpsilonMass", BIN, fMassLow, fMassHigh); 
	h = new TH1D("UpsilonMasscut","UpsilonMass, Pt>10GeV, |dpz|<12",BIN,fMassLow,fMassHigh);
	h = new TH1D("CandMass", "M_{#mu#mu}", 44, 1, 12.);
	h = new TH1D("CandPhi", "Dimuon Phi", 100,-4,4);
	h = new TH1D("CandPt", "P_{T}(#mu#mu) GeV", 50, 0, 50.);
	h = new TH1D("CandRapidity", "y(#mu#mu)", 60, 0, 3.);
	h = new TH1D("CandEta", "#eta(#mu#mu)", 100, -5, 5.);
	((TH1D*)fpHistFile->FindObjectAny("UpsilonMass"))->Sumw2(); 
	h = new TH1D("SigMuEta", "Single Muon #eta", 60, -3, 3.);
	h = new TH1D("SigMuPt", "Single Muon P_{T} GeV", 100, 0, 50.);
	k = new TH2D("SigMuEtaPt", "Single Muon #Eta P_{T}", 48, -2.4, 2.4, 50, 0, 50);
	k = new TH2D("CandMuPt", "CandMuPt", 50, 0, 50., 50, 0., 50.);
	
	h = new TH1D("Y1Pt", "unbinned P_{T} distribution for #Upsilon(1S), |y|<0.6", 100,0,100); // Y1s histogram
	h = new TH1D("Y2Pt", "unbinned P_{T} distribution for #Upsilon(2S), |y|<0.6", 100,0,100); 
	h = new TH1D("Y3Pt", "unbinned P_{T} distribution for #Upsilon(3S), |y|<0.6", 100,0,100); 
	
	h = new TH1D("SB_L_Pt", "unbinned P_{T} distribution for sideband low, |y|<0.6", 100,0,100); 
	h = new TH1D("SB_H_Pt", "unbinned P_{T} distribution for sideband high, |y|<0.6", 100,0,100); 

	
	h = new TH1D("Y1Pt_pre", "unbinned P_{T} distribution for #Upsilon(1S)pre-cut", 100,0,100); // Y1S before cuts
	h = new TH1D("Y2Pt_pre", "unbinned P_{T} distribution for #Upsilon(2S)pre-cut", 100,0,100); 
	h = new TH1D("Y3Pt_pre", "unbinned P_{T} distribution for #Upsilon(3S)pre-cut", 100,0,100); 
	
	
	h = new TH1D("CandRapidity", "y(#mu#mu)", 60, 0, 3.);
	h = new TH1D("CandRapidity_pre","y(#mu#mu) pre-cut",60,0,3.);
	h = new TH1D("CandEta", "#eta(#mu#mu)", 100, -5, 5.);
	h = new TH1D("UpsilonMass", "M_{#mu#mu} summary histogram", BIN, fMassLow, fMassHigh);
	h = new TH1D("UpsilonMass_pre", "M_{#mu#mu} summary histogram pre-cut", BIN, fMassLow, fMassHigh); // pre cuts  
	
	h = new TH1D("dm","Error on Mass, no binning in y or P_{T}",ndm_bins,0,dm_max);
	for (int iy=0; iy<fNy2; iy++) {
		k = new TH2D(Form("m_dm_y%d",iy),Form("dm(MeV) vs m(GeV) no P_{T} binning, %.1f<|y|<%.1f",fYbin2[iy],fYbin2[iy+1]), nm_bins,fMassLow,fMassHigh,ndm_bins,0,dm_max);
		k = new TH2D(Form("m_dmw_y%d",iy),Form("dm(MeV) vs m(GeV) no P_{T} binning, %.1f<|y|<%.1f weighted",fYbin2[iy],fYbin2[iy+1]), nm_bins,fMassLow,fMassHigh,ndm_bins,0,dm_max);
		k = new TH2D(Form("m_dmw_Ep_y%d",iy),Form("dm(MeV) vs m(GeV) no P_{T} binning, %.1f<|y|<%.1f weighted +1#sigma",fYbin2[iy],fYbin2[iy+1]), nm_bins,fMassLow,fMassHigh,ndm_bins,0,dm_max);
		k = new TH2D(Form("m_dmw_Em_y%d",iy),Form("dm(MeV) vs m(GeV) no P_{T} binning, %.1f<|y|<%.1f weighted -1#sigma",fYbin2[iy],fYbin2[iy+1]), nm_bins,fMassLow,fMassHigh,ndm_bins,0,dm_max);
	}


	

	k = new TH2D("dm_nPv", "dm(MeV) vs. Number of PV", 30,0,30,ndm_bins,0,dm_max); 
	
    fpHistFile->cd("kinematics");
	
	//Histograms for making PDFs 

	h = new TH1D("dphi","Difference in dimuon phi.", 100,-4,4);  

	k = new TH2D("muonPxPy", "Muon Px-Py", 100,-100,100,100,-100,100);
	h = new TH1D("muonphi","Single Muon Phi", 100,-4,4);

	k = new TH2D("phi1_phi2", "phi1 vs phi2", 100,-4,4,100,-4,4);
	k = new TH2D("pt_phi", " Pt vs phi", 100,-4,4, 100, 0, 100); 

	// Histograms for various quantities which we will cut on
	
	fpHistFile->cd("various_histograms"); 
	
	// Single muon Pt
	h = new TH1D("single_muon_pass","Single Muon Pt, Passed Cuts.", 100, 0, 200);
	h = new TH1D("single_muon_failed","Single Muon Pt, Failed Cuts.", 100, 0, 200); 
	
	k = new TH2D("doca_prob","Doca vs Prob", 100,0,1, 100,0,0.25); 
	((TH2D*)fpHistFile->FindObjectAny("doca_prob"))->GetXaxis()->SetTitle("Vertex Probability"); 
	((TH2D*)fpHistFile->FindObjectAny("doca_prob"))->GetYaxis()->SetTitle("Distance of Closest Approach"); 

	h = new TH1D("max_doca", "Distance of Closest Approach", 100,0,0.5); 
	h = new TH1D("dxy", "Track dxy", 50,-0.5,0.5);
	h = new TH1D("fTip", "Track dxy wrt primary vertex", 100,-5,5);
	h = new TH1D("fBsTip", "Track dxy wrt primary vertex", 100,-5,5);

	
	k = new TH2D("dxy_candphi","dxy vs cand-phi", 50,-1,1,100,-4,4);
	h = new TH1D("dz","dz Track",50,-20,20);
	h = new TH1D("fLip","dz Track wrt primary vertex",50,-30,30);
	h = new TH1D("fBsLip","dz Track wrt primary vertex",50,-30,30);


	h = new TH1D("trackchi2_dof","Track #chi^{2}/NDOF",50,0,5);
    h = new TH1D("muonchi2_dof", "Muon #chi^{2}/NDOF",50,0,5);
	h = new TH1D("vertexchi2_dof","Vertex #chi^{2}/NDOF",50,0,5);
	h = new TH1D("numberOfValidHits","Number of Valid Hits", 50,0,50); 
	h = new TH1I("bit4", "MuonID Bit 4:TrackerMuonArbitrated", 2,0,2); 
	h = new TH1I("numberOfPixelLayersWithMeasurement", "Numer Of Pixel Layers With Measurement", 24,0,6); 
	h = new TH1I("bit12","MuonID Bit 12:TMOneStationTight",2,0,2);     
	h = new TH1D("vertex_prob", "Vertex Probability ", 50,0,1.1);
	h = new TH1D("vertex_significance", "Vertex f3D/f3DE", 50,0,20);

	fpHistFile->cd(); 

	fpHistFile->cd("cuts"); 
	
	//Mass histograms for each cut
	//1D 
	h = new TH1D("mProb","Probability cut",BIN, fMassLow, fMassHigh);
	h = new TH1D("mVertexSig","Vertex Distance Selection",BIN, fMassLow, fMassHigh);
	h = new TH1D("mIndex","cut on index numbers",BIN, fMassLow, fMassHigh);
	h = new TH1D("mdz","cut on dz", BIN, fMassLow, fMassHigh);
	h = new TH1D("mdxy","cut on dxy", BIN, fMassLow, fMassHigh);
	h = new TH1D("mHits","cuts on hits",BIN, fMassLow, fMassHigh);
	h = new TH1D("mChi","cut on Chi2",BIN, fMassLow, fMassHigh);
	h = new TH1D("mbits", "cut on bits", BIN, fMassLow, fMassHigh);   
	h = new TH1D("mPt_region1", "cut on Pt Region1", BIN, fMassLow, fMassHigh);
	h = new TH1D("mPt_region2", "cut on Pt Region2", BIN, fMassLow, fMassHigh);
	h = new TH1D("meta", "cut on eta", BIN, fMassLow, fMassHigh); 
	h = new TH1D("mfQ", "cut on fQ", BIN, fMassLow, fMassHigh);    
	
	h = new TH1D("muonPtProb","Probability cut",100, 0, 200.);
	h = new TH1D("muonPtVertexSig","Decay Length cut",100, 0, 200.);
	h = new TH1D("muonPtIndex","cut on index numbers",100, 0, 200.);
	h = new TH1D("muonPtdz","cut on dz", 100, 0, 200.);
	h = new TH1D("muonPtHits","cuts on hits",100, 0, 200.);
	h = new TH1D("muonPtChi","cut on Chi2",100, 0, 200.);
	h = new TH1D("muonPtbits", "cut on bits", 100, 0, 200.);
	h = new TH1D("muonPt_region1", "cut on Pt Region1", 100, 0, 200.);
	h = new TH1D("muonPt_region2", "cut on Pt Region2", 100, 0, 200.);
	h = new TH1D("muonPteta", "cut on eta", 100, 0, 200.);
	h = new TH1D("muonPtfQ", "cut on fQ", 100, 0, 200.);
	
	h = new TH1D("cuts_sequence","Histogram of Cuts",50,0,20); 

	fpHistFile->cd();
	
  fpHistFile->cd("binned_ups");	

	//Upsilon Mass in new binning
	for ( int iy = 0; iy < fNy2; ++iy ){
		for ( int ipt = 0; ipt < fNpt1; ++ipt ){
			h = new TH1D(Form("UpsilonMass_y_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
						 Form("UpsilonMass_y_%.1f_%.1f,pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1]),
						 250,fMassLow, fMassHigh);  
			((TH1D*)fpHistFile->FindObjectAny(Form("UpsilonMass_y_%.1f_%.1f_pt%.1f_%.1f", UPSTYPE, fYbin2[iy], fYbin2[iy+1], fPTbin1[ipt], fPTbin1[ipt+1])))->GetXaxis()->SetTitle("M_{#mu#mu},GeV");
	    }//pt loop
	}//y loop
	
	//Upsilon Mass spectrum in bins of muon Pt 
	
	for ( int ipt = 1; ipt < fNpt; ++ipt ){
		h = new TH1D(Form("UpsilonMass_ptmu%.1f_%.1f", fPTbin[ipt], fPTbin[ipt+1]),
					 Form("UpsilonMass_ptmu%.1f_%.1f", fPTbin[ipt], fPTbin[ipt+1]),
					 250,fMassLow, fMassHigh);  
		((TH1D*)fpHistFile->FindObjectAny(Form("UpsilonMass_ptmu%.1f_%.1f", fPTbin[ipt], fPTbin[ipt+1])))->GetXaxis()->SetTitle("Mass,GeV");
	}//pt loop
	
	
}

void xsReader::readCuts(TString filename, int dump) {
  BIN=100;
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "Reading " << fCutFile.Data() << " for xsReader cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  char SetName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  
  int ibin; 

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f %s", CutName, &CutValue, SetName);
    
    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }
    
    if (!strcmp(CutName, "MODE")) {
      MODE = int(CutValue); ok = 1;
      if (dump) cout << "MODE:           " << MODE << endl;
    }
        
    if (!strcmp(CutName, "MUTYPE1")) {
      MUTYPE1 = int(CutValue); ok = 1;
      if (dump) cout << "MUTYPE1:           " << MUTYPE1 << endl;
    }    
    
    if (!strcmp(CutName, "MUTYPE2")) {
      MUTYPE2 = int(CutValue); ok = 1;
      if (dump) cout << "MUTYPE2:           " << MUTYPE2 << endl;
    }       
    
	  if (!strcmp(CutName, "PT_CORR")) {
		  PT_CORR = int(CutValue); ok = 1;
		  if (dump) cout << "PT_CORR:         " << PT_CORR << endl;
	  }
	  
	  if (!strcmp(CutName, "JPSI_Z_CORR")) {
		  JPSI_Z_CORR = CutValue; ok = 1;
		  if (dump) cout << "JPSI_Z_CORR:         " << JPSI_Z_CORR << endl;
	  }
	  if (!strcmp(CutName, "JPSI_Z_WEIGHT")) {
		  JPSI_Z_WEIGHT = CutValue; ok = 1;
		  if (dump) cout << "JPSI_Z_WEIGHT:         " << JPSI_Z_WEIGHT << endl;
	  }
	  
	  if (!strcmp(CutName, "PTLO")) {
		  PTLO = CutValue; ok = 1;
		  if (dump) cout << "PTLO:           " << PTLO << " GeV" << endl;
		  ibin = 11;
		  }
	  
	  if (!strcmp(CutName, "PTHI")) {
		  PTHI = CutValue; ok = 1;
		  if (dump) cout << "PTHI:           " << PTHI << " GeV" << endl;
		  ibin = 11;
		 
	  }    
	  
	  if (!strcmp(CutName, "PTMI")) {
		  PTMI = CutValue; ok = 1;
		  if (dump) cout << "PTMI:           " << PTMI << " GeV" << endl;
	  } 
	  
	  
	  if (!strcmp(CutName, "PTBARREL")) {
		  PTBARREL = CutValue; ok = 1;
		  if (dump) cout << "PTBARREL:        " << PTBARREL << " GeV" << endl;
	  }  
	  
	  if (!strcmp(CutName, "ETABARREL")) {
		  ETABARREL = CutValue; ok = 1;
		  if (dump) cout << "ETABARREL:       " << ETABARREL  << endl;
	  }     
	  
	  if (!strcmp(CutName, "ETAMI")) {
		  ETAMI = CutValue; ok = 1;
		  if (dump) cout << "ETAMI:           " << ETAMI << endl;
	  }
      if (!strcmp(CutName, "ETAHI")) {
		  ETAHI = CutValue; ok = 1;
		  if (dump) cout << "ETAHI:           " << ETAHI << endl;
		  ibin = 14;
		}
	  
	  if (!strcmp(CutName, "PTCAND_MIN")) {
		  PTCAND_MIN = CutValue; ok = 1;
		  if (dump) cout << "PTCAND_MIN:        " << PTCAND_MIN << " GeV" << endl;
	  }  
	  
	  if (!strcmp(CutName, "PTCAND")) {
		  PTCAND = CutValue; ok = 1;
		  if (dump) cout << "PTCAND:        " << PTCAND << " GeV" << endl;
	  }  
	  
	  if (!strcmp(CutName, "RAPCAND")) {
		  RAPCAND = CutValue; ok = 1;
		  if (dump) cout << "RAPCAND:        " << RAPCAND << " GeV" << endl;
	  }     
	  
	  
	  
    if (!strcmp(CutName, "RESTYPE")) {
      RESTYPE = int(CutValue); ok = 1;
      if (dump) cout << "RESTYPE:         " << RESTYPE << endl;
    }
    
    if (!strcmp(CutName, "UPSTYPE")) {
      UPSTYPE = int(CutValue); ok = 1;
      if (dump) cout << "UPSTYPE:         " << UPSTYPE << endl;
    }
    
    if (!strcmp(CutName, "MASSLO")) {
      MASSLO = CutValue; ok = 1;
      if (dump) cout << "MASSLO:          " << MASSLO << endl;
    }           
    
    if (!strcmp(CutName, "MASSHI")) {
      MASSHI = CutValue; ok = 1;
      if (dump) cout << "MASSHI:         " << MASSHI << endl;
    }   
    
    if (!strcmp(CutName, "DETA")) {
      DETA = CutValue; ok = 1;
      if (dump) cout << "DETA:           " << DETA << endl;
    } 
    
    if (!strcmp(CutName, "DPHI")) {
      DPHI = CutValue; ok = 1;
      if (dump) cout << "DPHI:           " << DPHI << endl;
    }
	  
	  if (!strcmp(CutName, "VPROB_MIN")) {
		  VPROB_MIN = CutValue; ok = 1;
		  if (dump) cout << "VPROB_MIN:           " << VPROB_MIN << endl;
	  }
	  if (!strcmp(CutName, "VERTEX_SIG")) {
		  VERTEX_SIG = CutValue; ok = 1;
		  if (dump) cout << "VERTEX_SIG:           " << VERTEX_SIG << endl;
	  }  
	  
	  if (!strcmp(CutName, "DZ_MAX")) {
		  DZ_MAX = CutValue; ok = 1;
		  if (dump) cout << "DZ_MAX:           " << DZ_MAX << endl;
	  }
	  
	  if (!strcmp(CutName, "DXY_MAX")) {
		  DXY_MAX = CutValue; ok = 1;
		  if (dump) cout << "DXY_MAX:           " << DXY_MAX << endl;
	  }
	  if (!strcmp(CutName, "NV_HITS_MIN")) {
		  NV_HITS_MIN = CutValue; ok = 1;
		  if (dump) cout << "NV_HITS_MIN:           " << NV_HITS_MIN << endl;
	  }
	  
	  if (!strcmp(CutName, "MUON_CHI2_DOF_MAX")) {
		  MUON_CHI2_DOF_MAX = CutValue; ok = 1;
		  if (dump) cout << "MUON_CHI2_DOF_MAX:           " << MUON_CHI2_DOF_MAX << endl;
	  }

    if (!strcmp(CutName, "BARREL")) {
      BARREL = int(CutValue); ok = 1;
      if (dump) cout << "BARREL:         " << BARREL << endl;
    }
    
		  
    if (!strcmp(CutName, "HLTPATH")) {
      HLTPATH = SetName; ok = 1;
      if (dump) cout << "HLTPATH:    " << HLTPATH  << endl;
    } 
    
    if (!strcmp(CutName, "HLTLABEL")) {
      HLTLABEL = SetName; ok = 1;
      if (dump) cout << "HLTLABEL:   " << HLTLABEL  << endl;
    }     
    
    if (!strcmp(CutName, "HLTPATH1")) {
      HLTPATH1 = SetName; ok = 1;
      if (dump) cout << "HLTPATH1:   " << HLTPATH1  << endl;
    }   
    
    if (!strcmp(CutName, "HLTPATH2")) {
      HLTPATH2 = SetName; ok = 1;
      if (dump) cout << "HLTPATH2:   " << HLTPATH2  << endl;
    }
     
    if (!strcmp(CutName, "HLTPATH3")) {
      HLTPATH3 = SetName; ok = 1;
      if (dump) cout << "HLTPATH3:   " << HLTPATH3  << endl;
    }   

    if (!strcmp(CutName, "HLTPATH4")) {
      HLTPATH4 = SetName; ok = 1;
       if (dump) cout << "HLTPATH4:   " << HLTPATH4  << endl;
    }      
    
    if (!strcmp(CutName, "HLTLABEL1")) {
      HLTLABEL1 = SetName; ok = 1;
      if (dump) cout << "HLTLABEL1:  " << HLTLABEL1  << endl;
    }     
     
    if (!strcmp(CutName, "HLTPATH5")) {
      HLTPATH5 = SetName; ok = 1;
      if (dump) cout << "HLTPATH5:   " << HLTPATH5  << endl;
    }   

    if (!strcmp(CutName, "HLTPATH6")) {
      HLTPATH6 = SetName; ok = 1;
       if (dump) cout << "HLTPATH6:   " << HLTPATH6  << endl;
    }     
    
    if (!strcmp(CutName, "HLTLABEL2")) {
      HLTLABEL2 = SetName; ok = 1;
      if (dump) cout << "HLTLABEL2:  " << HLTLABEL2  << endl;
    }     
     
    if (!strcmp(CutName, "HLTPATH7")) {
      HLTPATH7 = SetName; ok = 1;
      if (dump) cout << "HLTPATH7:   " << HLTPATH7  << endl;
    }   

    if (!strcmp(CutName, "HLTPATH8")) {
      HLTPATH8 = SetName; ok = 1;
       if (dump) cout << "HLTPATH8:   " << HLTPATH8  << endl;
    }     
    
    if (!strcmp(CutName, "HLTPATH9")) {
      HLTPATH9 = SetName; ok = 1;
      if (dump) cout << "HLTPATH9:   " << HLTPATH9  << endl;
    }   

    if (!strcmp(CutName, "HLTPATH10")) {
      HLTPATH10 = SetName; ok = 1;
       if (dump) cout << "HLTPATH10:   " << HLTPATH10  << endl;
    }       
    
    if (!strcmp(CutName, "HLTPATH11")) {
      HLTPATH11 = SetName; ok = 1;
       if (dump) cout << "HLTPATH11:   " << HLTPATH11  << endl;
    }  
    
    if (!strcmp(CutName, "HLTPATH12")) {
      HLTPATH12 = SetName; ok = 1;
       if (dump) cout << "HLTPATH12:   " << HLTPATH12  << endl;
    } 
    
    if (!strcmp(CutName, "HLTPATH13")) {
      HLTPATH13 = SetName; ok = 1;
       if (dump) cout << "HLTPATH13:   " << HLTPATH13  << endl;
    }  
        
    if (!strcmp(CutName, "HLTLABEL3")) {
      HLTLABEL3 = SetName; ok = 1;
      if (dump) cout << "HLTLABEL3:  " << HLTLABEL3  << endl;
    }
         
  
    
    if (!ok) cout << "==> ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;
}

