/*
 *  rho_script.c
 *  
 *
 *  Created by Benjamin Carlson on 5/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "rho_script.h"

void rho_script(){
	gROOT->SetBatch();
	set_output_tree();
	loop();
	writeHisto();
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
	fTree->Fill();
	

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
	
	float VP=vProb;
	float DM1=distM1;
	
	fTree->Branch("vProb",&VP,"vProb/F");
	fTree->Branch("distM1",&DM1,"distM1/F");


}
void set_branches(TTree *tree){

	
	tree->SetBranchAddress("runNb",&runNb); 
	tree->SetBranchAddress("JpsiP", &Ups_P);
	tree->SetBranchAddress("muPosP", &Mu1_P);
	tree->SetBranchAddress("muNegP", &Mu2_P);
	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v1",&D5_v1);
	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v2",&D5_v2);
	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3",&D5_v3);
	tree->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v5",&D5_v5);

	tree->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v1",&D7_v1);
	tree->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v4",&D7_v4);
	
	tree->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v1",&D9_v1);
	tree->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v4",&D9_v4);


	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_v1", &D0_v1); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_v2", &D0_v2); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_v3", &D0_v3); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_v5", &D0_v5);
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_v6", &D0_v6); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_v9", &D0_v9); 
	tree->SetBranchAddress("HLT_DoubleMu3_Upsilon_v1", &DMu3_v1); 

	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v1", &DUM0_v1);
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v2", &DUM0_v2); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v3", &DUM0_v3); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v4", &DUM0_v4); 
	tree->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v6", &DUM0_v6); 

	
	tree->SetBranchAddress("HLT_Mu5_L2Mu2_v1", &HLT_Mu5_L2Mu2_v1);
	tree->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6", &HLT_Dimuon10_Jpsi_Barrel_v6); 

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
	if(fabs(Ups_P->Rapidity())<1.2) y_=true;
	if(Ups_P->M()>2.8 && Ups_P->M()<3.35) m_=true;
	if(Mu1_P->X()>900 || Mu2_P->X()>900) return false;
	if(Mu1_P->Perp()<3.5 || Mu2_P->Perp()<3.5) return false;
	
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
	
	bool leadingmu_pt=false;
	
	bool SG=false;
	
	if(Ups_P->Perp()>10 && Ups_P->Perp()<100) pt_=true;
	if(fabs(Ups_P->Rapidity())<1.2) y_=true;
	if(Ups_P->M()>8.5 && Ups_P->M()<11.5) m_=true;
	if(Mu1_P->X()>900 || Mu2_P->X()>900) return false;
	if(Mu1_P->Perp()<3. || Mu2_P->Perp()<3.) return false;
	
	double ptmax=TMath::Max(Mu1_P->Perp(),Mu2_P->Perp());
	if(ptmax>20)leadingmu_pt=true; 
	
	if(Mu1_P->DeltaPhi(*Mu2_P)<0) SG=true;
	
	if((fabs(Mu1_P->Eta())<1.2 && Mu1_P->Perp()>4.5) || ((fabs(Mu1_P->Eta())>1.2 && fabs(Mu1_P->Eta())<1.4) && Mu1_P->Perp()>3.5) || ((fabs(Mu1_P->Eta())>1.2 && fabs(Mu1_P->Eta())<1.4) && Mu1_P->Perp()>3.0))
		mu1_=true;
	if((fabs(Mu2_P->Eta())<1.2 && Mu2_P->Perp()>4.5) || ((fabs(Mu2_P->Eta())>1.2 && fabs(Mu2_P->Eta())<1.4) && Mu2_P->Perp()>3.5) || ((fabs(Mu2_P->Eta())>1.2 && fabs(Mu2_P->Eta())<1.4) && Mu2_P->Perp()>3.0))
		mu2_=true;
	
	bool good=false;
	if(pt_ && y_ && m_ && mu1_ && mu2_ && SG && leadingmu_pt) good=true; 
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
	
	//string name =Form("gEff_DATA_PT_AETA%d",bin); 
	string name= Form("gEff_nominal_AETA%d",bin);//nominal

	if(print) cout << name << endl; 
	fEffL1L2->cd();
	if(fEffL1L2->IsOpen()!=1 || fEffL3->IsOpen()!=1 || fEff->IsOpen()!=1) {
		cout << "Efficiency File not open. " << endl; 
		return 0; 
	}
		
	TGraph *grL1L2=(TGraph*)fEffL1L2->Get(name.c_str()); 
	TGraph *grL3=(TGraph*)fEffL3->Get(name.c_str()); 
	TGraph *gr=(TGraph*)fEff->Get(name.c_str()); 
	
	//double emu=grL1L2->Eval(pt)*grL3->Eval(pt); 
	double emu=gr->Eval(pt)*0.99/0.97;
	double RV=-9; 
	
	if(sys=="nom") RV=emu; 
	
	if(print) cout << "Pt: " << pt << "effmu: " << emu << endl; 
	

	delete grL1L2;
	delete grL3; 
	delete gr; 
	if(RV<0) cout << "return value <0" << endl; 
	return RV; 
}


void loop(){

	TFile *F=new TFile("TTree_SingleMu_Onia2MuMu_v20_PromptReco_AB.root","READ"); 
	TTree *tree = (TTree*)F->FindObjectAny("data"); 
	int entries = tree->GetEntries(); 
	set_branches(tree);
	
	double pt_bins[]={10,20,30,40,50,60,70,100};
	double dR_bins[]={0,0.25,0.3,0.35,0.4,0.8,1,1.5};
	int ndR=sizeof(dR_bins)/sizeof(double)-1; 
	int npT=sizeof(pt_bins)/sizeof(double)-1;
	//TEfficiency *rho = new TEfficiency("TEff_rho_pt","#rho p_{T}; p_{T}(#mu#mu) [GeV]; #rho",7,pt_bins); 
	
	CreateHistogram("hpt_lead","p_{T}(#mu) of higher p_{T} muon","p_{T}(#mu) [GeV]","#rho",100,0,100); 
	CreateHistogram("hpt_trail","p_{T}(#mu) of lower p_{T} muon","p_{T}(#mu) [GeV]","#rho",100,0,100); 

	CreateHistogram("h_mass_pt0","M_{#mu#mu}, p_{T}<35","M_{#mu#mu}","Events",100,8.5,11.5); 
	CreateHistogram("h_mass_pt1","M_{#mu#mu}, p_{T}>35","M_{#mu#mu}","Events",100,8.5,11.5); 

	CreateHistogram("h_massTrig_pt0","M_{#mu#mu}, p_{T}<35","M_{#mu#mu}","Events",100,8.5,11.5); 
	CreateHistogram("h_massTrig_pt1","M_{#mu#mu}, p_{T}>35","M_{#mu#mu}","Events",100,8.5,11.5); 
	
	CreateHistogram("DeltaRPtE_num","#DeltaR_p_{T}^{ellpictic}","#DeltaR_p_{T}^{ellpictic}","",100,0,3); 
	CreateHistogram("DeltaRPtE_num_uW","#DeltaR_p_{T}^{ellpictic}","#DeltaR_p_{T}^{ellpictic}","",100,0,3); 
	CreateHistogram("DeltaRPtE_den","#DeltaR_p_{T}^{ellpictic}","#DeltaR_p_{T}^{ellpictic}","",100,0,3); 

	
	CreateHistogram("hpt_den","#rho denominator","p_{T} [GeV]","",9,10,100); 
	CreateHistogram("hpt_num_uW","#rho numerator","p_{T} [GeV]","",9,10,100); 
	CreateHistogram("hpt_num","#rho numerator","p_{T} [GeV]","",9,10,100); 

	CreateHistogram("mass","mass","M_{#mu#mu}","",100,8.5,11.5); 
	
	CreateHistogram("mass_DeltaRPtE_num","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,100,0,3);
	CreateHistogram("mass_DeltaRPtE_den","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,100,0,3);

	
	CreateHistogram("mass_y0_DeltaRPtE_numBin","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,ndR,dR_bins);
	CreateHistogram("mass_y0_DeltaRPtE_numBin_uW","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,ndR,dR_bins);
	CreateHistogram("mass_y0_DeltaRPtE_denBin","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,ndR,dR_bins);

	CreateHistogram("mass_y1_DeltaRPtE_numBin","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,ndR,dR_bins);
	CreateHistogram("mass_y1_DeltaRPtE_numBin_uW","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,ndR,dR_bins);
	CreateHistogram("mass_y1_DeltaRPtE_denBin","M_{#mu#mu}","M_{#mu#mu}","#DeltaR_p_{T}^{ellpictic}",100,8.5,11.5,ndR,dR_bins);


	CreateHistogram("mass_pt","M_{#mu#mu}","p_{T} [GeV]","",100,8.5,11.5,npT,pt_bins);
	CreateHistogram("mass_pt_num","M_{#mu#mu}","p_{T} [GeV]","",100,8.5,11.5,npT,pt_bins);
	CreateHistogram("mass_pt_num_uW","M_{#mu#mu}","p_{T} [GeV]","",100,8.5,11.5,npT,pt_bins);

	CreateHistogram("DeltaRPtE_pt","#DeltaR_p_{T}^{ellpictic} vs p_{T}","p_{T} [GeV]", "#DeltaR_p_{T}^{ellpictic}",9,10,100,20,0,3); 
	CreateHistogram("DeltaRPtE_pt_Trig","#DeltaR_p_{T}^{ellpictic} vs p_{T}","p_{T} [GeV]", "#DeltaR_p_{T}^{ellpictic}",9,10,100,20,0,3); 


	cout << "Entries: " << entries << endl;
	//entries=1000;
	for(int i=0; i<entries;i++){
		if(i%100000==0) {
			cout << "Event: " << i << endl;
			cout << "Run Nb: " << runNb << endl; 
		}
		tree->GetEntry(i); 

		if(D0_v1==1 || D0_v2==1 || D0_v3==1 || D0_v5==1 || D0_v6==1 || D0_v9==1 || DMu3_v1==1 || DUM0_v1==1 || DUM0_v2==1 || DUM0_v3==1 || DUM0_v4==1 || DUM0_v6==1) continue;
		if(runNb<165088 || runNb>179889) continue;
		if(Mu1_P->X()>900 || Mu2_P->X()>900) continue;
		fill_tree();
		if(!select_upsilons()) continue; 
//		if(!select_jpsi()) continue; 

		double DeltaPhi=TMath::Abs(Mu1_P->DeltaPhi(*Mu2_P)); 
		double DeltaEta=TMath::Abs(Mu1_P->Eta()-Mu2_P->Eta());
		double DeltaR=Mu1_P->DeltaR(*Mu2_P);
		double DeltaPt=TMath::Abs(Mu1_P->Perp()-Mu2_P->Perp());
		double DeltaRPtE=TMath::Sqrt(TMath::Power(0.00157*DeltaPt,2)+TMath::Power(1.2*DeltaPhi,2)+TMath::Power(DeltaEta,2));//deltaR elliptic_deltaPt

		if(Ups_P->Perp()<35)hName["h_mass_pt0"]->Fill(Ups_P->M()); 
		else hName["h_mass_pt1"]->Fill(Ups_P->M()); 
		
		hName["mass"]->Fill(Ups_P->M()); 
		hName2D["DeltaRPtE_pt"]->Fill(Ups_P->Perp(),DeltaRPtE); 
		
		hName2D["mass_pt"]->Fill(Ups_P->M(),Ups_P->Perp()); 
		int iy=-1;
		if(fabs(Ups_P->Rapidity()<0.6)) iy=0; 
		else iy=1; 
		hName2D["mass_DeltaRPtE_den"]->Fill(Ups_P->M(),DeltaRPtE); 
		hName2D[Form("mass_y%d_DeltaRPtE_denBin",iy)]->Fill(Ups_P->M(),DeltaRPtE); 
		if(Ups_P->M()>9.26 && Ups_P->M()<9.66){
			hName["hpt_den"]->Fill(Ups_P->Perp());
			hName["DeltaRPtE_den"]->Fill(DeltaRPtE);
		}


		double ptmin=TMath::Min(Mu1_P->Perp(),Mu2_P->Perp());
		double ptmax=TMath::Max(Mu1_P->Perp(),Mu2_P->Perp()); 
		
		hName["hpt_lead"]->Fill(ptmax);
		hName["hpt_trail"]->Fill(ptmin); 

		double etamin=0; 
		
		if(ptmin<Mu1_P->Perp()) etamin=Mu2_P->Eta();
		else etamin=Mu1_P->Eta();
		
		double dimuon_eff=muonEff(ptmin,etamin,"nom");
		bool Trig=false; 
		if(D5_v1==1 || D5_v2==1 || D5_v3==1 || D5_v5==1 || D7_v1==1 || D7_v4==1 || D9_v1==1 || D9_v4==1) Trig=true; 
		if(Trig && Ups_P->M()>9.26 && Ups_P->M()<9.66){
			hName["hpt_num_uW"]->Fill(Ups_P->Perp()); 
			hName["hpt_num"]->Fill(Ups_P->Perp(),1./dimuon_eff);
			hName["DeltaRPtE_num"]->Fill(DeltaRPtE,1./dimuon_eff);
			hName["DeltaRPtE_num_uW"]->Fill(DeltaRPtE); 
		}
		if(Trig){
			if(Ups_P->Perp()<35)hName["h_massTrig_pt0"]->Fill(Ups_P->M()); 
			else hName["h_massTrig_pt1"]->Fill(Ups_P->M()); 

			hName2D["mass_pt_num"]->Fill(Ups_P->M(),Ups_P->Perp(),1.0/dimuon_eff);
			hName2D["mass_pt_num_uW"]->Fill(Ups_P->M(),Ups_P->Perp());
			hName2D["mass_DeltaRPtE_num"]->Fill(Ups_P->M(),DeltaRPtE,1.0/dimuon_eff); 
			hName2D[Form("mass_y%d_DeltaRPtE_numBin",iy)]->Fill(Ups_P->M(), DeltaRPtE,1.0/dimuon_eff);
			hName2D[Form("mass_y%d_DeltaRPtE_numBin_uW",iy)]->Fill(Ups_P->M(), DeltaRPtE);
			hName2D["DeltaRPtE_pt_Trig"]->Fill(Ups_P->Perp(),DeltaRPtE);

		}

		//if(Ups_P->M()>9.26 && Ups_P->M()<9.66)rho->Fill(Trig,Ups_P->Perp());
	//	rho->SetWeight(0.97/dimuon_eff);
		
		//if(HLT_Dimuon7_Upsilon_Barrel_v1) hpt_num->Fill(Ups_P->Perp(),1./dimuon_eff); 
		
	}
		
}


void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Double_t* xBins)
{
	TH1F* h = new TH1F(name, title, nBinsX, xBins);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}

void CreateHistogram(const char* name,   const char* title,
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

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, Double_t xLow, Double_t xUp,
					Int_t nBinsY, const Double_t* yBins)
{
	TH2F* h = new TH2F(name, title, nBinsX, xLow,xUp,nBinsY,yBins);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName2D[name] = h;
}


void writeHisto(){
	fOut->cd();
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		hName[it->first]->Write();
	}
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++) {
		hName2D[it->first]->Write();
	}
	
}
