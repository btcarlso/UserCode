/*
 *  rho_script.c
 *  
 *
 *  Created by Benjamin Carlson on 5/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "rho_script_MC.h"

void rho_script_MC(){
	gROOT->SetBatch();
	set_output_tree();
	loop();
	fTree->Write();
	fOut->Close();
	
}

void fill_tree(){

	m=Ups_P->M();
	y=Ups_P->Rapidity();
	pt=Ups_P->Perp();
	etaMuP=Mu1_P->Eta();
	etaMuM=Mu2_P->Eta();
	ptMuP=Mu1_P->Perp();
	ptMuM=Mu2_P->Perp();
	if(Mu1_P->DeltaPhi(*Mu2_P)<0) SeaGull=1;
	else SeaGull=0; 
	Trig5=0; Trig7=0; Trig9=0; 
	if(D5_v1==1 || D5_v2==1 || D5_v3==1 || D5_v5==1) Trig5=1; 
	if(D7_v1==1 || D7_v4==1) Trig7=1; 
	if(D9_v1==1 || D9_v4==1) Trig9=1; 

	

}

void set_output_tree(){
	fTree = new TTree("data","data");
	fTree->Branch("JpsiP",&m,"JpsiP/F");
	fTree->Branch("JpsiPt",&pt,"JpsiPt/F");
	fTree->Branch("JpsiRapidity",&y,"JpsiRapidity/F");
	fTree->Branch("MuP_Pt",&ptMuP,"MuP_Pt/F");
	fTree->Branch("MuM_Pt",&ptMuM,"MuM_Pt/F");
	fTree->Branch("MuP_eta",&etaMuP,"MuP_eta/F");
	fTree->Branch("MuM_eta",&etaMuM,"MuM_eta/F");

	fTree->Branch("SeaGull",&SeaGull,"SeaGull/I");
	fTree->Branch("Trig5",&Trig5,"Trig5/I");
	fTree->Branch("Trig7",&Trig7,"Trig5/I");
	fTree->Branch("Trig9",&Trig9,"Trig5/I");
	
	fTree->Branch("vProb",&vProb,"vProb/F");
	fTree->Branch("distM1",&distM1,"distM1/F");


}
void set_branches(TTree *tree){

	
	tree->SetBranchAddress("runNb",&runNb); 
	tree->SetBranchAddress("JpsiP_Gen", &Ups_P);
	tree->SetBranchAddress("muPosP_Gen", &Mu1_P);
	tree->SetBranchAddress("muNegP_Gen", &Mu2_P);
	
	tree->SetBranchAddress("JpsiP", &Ups_P_reco);
	tree->SetBranchAddress("muPosP", &Mu1_P_reco);
	tree->SetBranchAddress("muNegP", &Mu2_P_reco);

	tree->SetBranchAddress("muPos_qual", &muPos_qual);
	tree->SetBranchAddress("muNeg_qual", &muNeg_qual);

	tree->SetBranchAddress("NCand", &NCand);

	
	
//	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v1",&D5_v1);
//	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v2",&D5_v2);
	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3",&D5_v3);
	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v5",&D5_v5);

	tree->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v1",&D7_v1);
//	tree->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v4",&D7_v4);
	
	tree->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v1",&D9_v1);
//	tree->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v4",&D9_v4);



	tree->SetBranchAddress("JpsiVprob", &vProb); 
	tree->SetBranchAddress("JpsiDistM1",&distM1);

}

bool select_jpsi(){
	bool pt_=false; 
	bool y_=false;
	bool m_=false;
	
	bool mu1_=false; 
	bool mu2_=false; 
	
	bool SG=false;
	
	if(Ups_P->Perp()>10 && Ups_P->Perp()<100) pt_=true;
	if(fabs(Ups_P->Rapidity())<0.6) y_=true;
	if(Ups_P->M()>2.8 && Ups_P->M()<3.35) m_=true;
	if(Mu1_P->X()>900 || Mu2_P->X()>900) return false;
	if(Mu1_P->Perp()<4.5 || Mu2_P->Perp()<4.5) return false;
	
	if(Mu1_P->DeltaPhi(*Mu2_P)<0) SG=true;
	
	if((fabs(Mu1_P->Eta())<1.2 && Mu1_P->Perp()>4.5) || ((fabs(Mu1_P->Eta())>1.2 && fabs(Mu1_P->Eta())<1.4) && Mu1_P->Perp()>4.0) || ((fabs(Mu1_P->Eta())>1.2 && fabs(Mu1_P->Eta())<1.4) && Mu1_P->Perp()>3.5))
		mu1_=true;
	if((fabs(Mu2_P->Eta())<1.2 && Mu2_P->Perp()>4.5) || ((fabs(Mu2_P->Eta())>1.2 && fabs(Mu2_P->Eta())<1.4) && Mu2_P->Perp()>4.0) || ((fabs(Mu2_P->Eta())>1.2 && fabs(Mu2_P->Eta())<1.4) && Mu2_P->Perp()>3.5))
		mu2_=true;
	
	bool good=false;
	if(pt_ && y_ && m_ && mu1_ && mu2_ && SG) good=true; 
	return good;
}

bool select_upsilons(){
	bool pt_=false; 
	bool y_=false;
	bool m_=false;
	
	bool mu1_=false; 
	bool mu2_=false; 
	
	bool rho1=false;
	bool rho2=false; 
	
	bool leadingmu_pt=false;
	
	bool SG=true;
	
	if(Ups_P->Perp()>10 && Ups_P->Perp()<100) pt_=true;
	if(fabs(Ups_P->Rapidity())<0.6) y_=true;
	if(Ups_P->M()>8.5 && Ups_P->M()<11.5) m_=true;
	if(Mu1_P->X()>900 || Mu2_P->X()>900) return false;
	if(Mu1_P->Perp()<3.5 || Mu2_P->Perp()<3.5) return false;
	
	double ptmax=TMath::Max(Mu1_P->Perp(),Mu2_P->Perp());
	if(ptmax>20)leadingmu_pt=true; 
	
	//if(Mu1_P->DeltaPhi(*Mu2_P)<0) SG=true;
	
	if((fabs(Mu1_P->Eta())<1.2 && Mu1_P->Perp()>4.5) || ((fabs(Mu1_P->Eta())>1.2 && fabs(Mu1_P->Eta())<1.4) && Mu1_P->Perp()>4.0) || ((fabs(Mu1_P->Eta())>1.2 && fabs(Mu1_P->Eta())<1.4) && Mu1_P->Perp()>3.5))
		mu1_=true;
	if((fabs(Mu2_P->Eta())<1.2 && Mu2_P->Perp()>4.5) || ((fabs(Mu2_P->Eta())>1.2 && fabs(Mu2_P->Eta())<1.4) && Mu2_P->Perp()>4.0) || ((fabs(Mu2_P->Eta())>1.2 && fabs(Mu2_P->Eta())<1.4) && Mu2_P->Perp()>3.5))
		mu2_=true;
	
	if(Mu1_P->Rho()>2.5) rho1=true;
	if(Mu2_P->Rho()>2.5) rho2=true;
	
	
	bool good=false;
	if(pt_ && y_ && m_ && mu1_ && mu2_ && rho1 && rho2) good=true; 
	return good;
}
int eta_to_bin(double eta){
	
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
double muonEff(double pt,double eta,string sys){

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
	

	delete gr;
	if(RV<0) cout << "return value <0" << endl; 
	return RV; 
}


void chain(){
	//onia2MuMu_tree_13_1_G4M.root
	tree = new TChain("data"); 
	tree->Add("dcap:///pnfs/cms/WAX/11/store/user/btcarlso/UpsilonAnalysis/Ups1S_dz30/onia2MuMu_tree*.root");
}

void loop(){

//	TFile *F=new TFile("TTree_SingleMu_Onia2MuMu_v20_PromptReco_AB.root","READ"); 
//	TTree *tree = (TTree*)F->FindObjectAny("data"); 
	chain();
	Long64_t entries = tree->GetEntries(); 
	set_branches(tree);
	
	float pt_bins[]={10,20,30,40,50,60,70,100};
	
	
	TEfficiency *rho = new TEfficiency("TEff_rho_pt","#rho; p_{T}(#mu#mu) [GeV]; #rho",25,10,100);

	
	TEfficiency *seagull = new TEfficiency("TEff_sg_pt","#epsilon_{sg} ; p_{T}(#mu#mu) [GeV]; #epsilon_{sg}",25,10,100); 
	TEfficiency *VP_eff = new TEfficiency("TEff_VP","#epsilon_{vp} ; p_{T}(#mu#mu) [GeV]; #epsilon_{vp}",25,10,100); 
	
	TEfficiency *MuQualP = new TEfficiency("MuQualP","#epsilon_{#mu+} ; p_{T}(#mu) [GeV]; #epsilon_{#mu+}",100,0,100); 
	TEfficiency *MuQualN = new TEfficiency("MuQualN","#epsilon_{#mu- ; p_{T}(#mu) [GeV]; #epsilon_{#mu-}",100,0,100); 

	
	TH1D *hpt_lead = new TH1D("hpt_lead","p_{T} of leading Muon;p_{T} [GeV]",100,0,100); 
	TH1D *hpt_trail = new TH1D("hpt_trail","p_{T} of trailing Muon;p_{T} [GeV]",100,0,100); 

	TH1D *hpt_den = new TH1D("pt_den","p_{T}",50,10,100); 
	TH1D *hpt_num_uW = new TH1D("hpt_num_uW",";p_{T}",50,10,100);
	TH1D *hpt_num = new TH1D("pt_num",";p_{T}",50,10,100); 
	
	TH1D *hdR_num = new TH1D("dR_num",";p_{T}",50,0,3); 
	TH1D *hdR_den = new TH1D("dR_den",";p_{T}",50,0,3); 

	
	
	TH1D *mass = new TH1D("mass","M_{#mu#mu}",100, 8.5,11.5); 
	TH2D *mass_pt = new TH2D("mass_pt",";M_{#mu#mu};p_{T}(#mu#mu)", 100,8.5,11.5,9,10,100);
	TH2D *mass_pt_num = new TH2D("mass_pt_num",";M_{#mu#mu};p_{T}(#mu#mu)", 100,8.5,11.5,9,10,100);
	
	TH1F *h_vp=new TH1F("h_vp","Vertex Probability", 1000,0,100); 

	cout << "Entries: " << entries << endl;
	entries=1000000;
	for(Long64_t i=0; i<entries;i++){
		if(i%100000==0) {
			cout << "Event: " << i << endl;
			cout << "Run Nb: " << runNb << endl; 
		}
		tree->GetEntry(i); 

		if(D0_v1==1 || D0_v2==1 || D0_v3==1 || D0_v5==1 || D0_v6==1 || D0_v9==1 || DMu3_v1==1 || DUM0_v1==1 || DUM0_v2==1 || DUM0_v3==1 || DUM0_v4==1 || DUM0_v6==1) continue;
		//if(runNb<165088 || runNb>179889) continue;
		if(Mu1_P->X()>900 || Mu2_P->X()>900) continue;

		
		fill_tree();
		if(!select_upsilons()) continue; 
		
		MuQualP->Fill((muPos_qual==-9),Mu1_P->Perp()); 
		MuQualN->Fill((muNeg_qual==-9),Mu2_P->Perp()); 
		
		if(muPos_qual==-9 || muNeg_qual==-9) continue;
		
		bool Trig=false; 
		if(D5_v1==1 || D5_v2==1 || D5_v3==1 || D5_v5==1 || D7_v1==1 || D7_v4==1 || D9_v1==1 || D9_v4==1) Trig=true; 
		bool SG=false; 
		if(Mu1_P->DeltaPhi(*Mu2_P)<0) SG=true; 
		
		if(Trig)seagull->Fill(SG,Ups_P->Perp());
		
		if(!SG) continue;
		bool good_vertex=false; 
		if(vProb>0.005) good_vertex=true; 
		VP_eff->Fill(good_vertex,Ups_P->Perp());
		h_vp->Fill(vProb*100);
		
		if(vProb<0.005) continue;
		
		mass->Fill(Ups_P->M()); 
		mass_pt->Fill(Ups_P->M(),Ups_P->Perp()); 
		if(Ups_P->M()>9.26 && Ups_P->M()<9.66){ 	
			hpt_den->Fill(Ups_P->Perp());
			hdR_den->Fill(Mu1_P->DeltaR(*Mu2_P));
		}

		double ptmin=TMath::Min(Mu1_P->Perp(),Mu2_P->Perp());
		double ptmax=TMath::Max(Mu1_P->Perp(),Mu2_P->Perp()); 
		
		hpt_lead->Fill(ptmax);
		hpt_trail->Fill(ptmin); 

		
		double dimuon_eff=muonEff(Mu1_P->Perp(),Mu1_P->Eta(),"nom")*muonEff(Mu2_P->Perp(),Mu2_P->Eta(),"nom");
		bool Mu_qual = false;
		if(muPos_qual && muNeg_qual) Mu_qual=true; 
		
		if(Trig && Mu_qual && Ups_P->M()>9.26 && Ups_P->M()<9.66) {
			hpt_num_uW->Fill(Ups_P->Perp()); 
			hdR_num->Fill(Mu1_P->DeltaR(*Mu2_P),1.0/dimuon_eff); 
			hpt_num->Fill(Ups_P->Perp(),1.0/dimuon_eff);
		}
		
		
		if(Trig) mass_pt_num->Fill(Ups_P->M(),Ups_P->Perp(),1./dimuon_eff);
		
		rho->Fill((Trig && Mu_qual),Ups_P->Perp());
		
		
	//	rho->SetWeight(0.97/dimuon_eff);
		
		//if(HLT_Dimuon7_Upsilon_Barrel_v1) hpt_num->Fill(Ups_P->Perp(),1./dimuon_eff); 
		
	}
	
	cout << "Vertex Prob: " << Form("%.2f",h_vp->Integral(h_vp->FindBin(1),h_vp->GetNbinsX())/h_vp->Integral())<< endl;
	cout << "Rho stat uncert: " << rho->GetEfficiencyErrorLow(25) << endl;
	
	
	TCanvas *tmp = new TCanvas();
	hpt_num->SetStats(kFALSE); 
	fOut->cd();
	hpt_den->Write();
	hpt_num->Write();
	hpt_lead->Write();
	hpt_trail->Write();
	hpt_num_uW->Write();
	hdR_num->Write();
	hdR_den->Write();
	seagull->Write();
	VP_eff->Write();
	
	hpt_num->Divide(hpt_den);
	hdR_num->Divide(hdR_den);
	hpt_num->SetLineColor(kRed);
	
	TCanvas *Efficiency_plots = new TCanvas("Efficiency_plots");
	Efficiency_plots->Divide(2,2);
	Efficiency_plots->cd(1);
	seagull->Draw();
	Efficiency_plots->cd(2);
	VP_eff->Draw();
	Efficiency_plots->cd(3);
	//rho->Draw();
	hpt_num->Draw();
	
	Efficiency_plots->cd(4);
	hdR_num->Draw();
	
	
	
	//MuQualP->Draw();
	//MuQualN->Draw("same");
	
	
	
	Efficiency_plots->Write();
	rho->Write();
	mass->Write();
	
	MuQualN->Write();
	MuQualP->Write();
	mass_pt->Write();
	h_vp->Write();
	mass_pt_num->Write();
	hpt_num->SetName("rho_pt"); 
	hpt_num->Divide(hpt_den); 
	hpt_num->Draw(); 
	tmp->Print("rho_factor.png");
	
	TCanvas *M = new TCanvas();
	mass->Draw();
	M->Print("dimuon_mass.png");
	
}