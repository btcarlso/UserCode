/*
 *  analyze_trees
 *  
 *
 *  Created by Benjamin Carlson on 6/28/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "analyze_trees.h"


void analyze_trees(int jobnumber){
	setup_files(jobnumber);
	
	bookHisto();
	open_files();
	
	divide_BW();
	writeHisto();
}

void initialize_tree(TChain *tree){
	entries=tree->GetEntries();
	
	tree->SetBranchAddress("jets_n",&nJets); 
	
	tree->SetBranchAddress("st",&st); 
	tree->SetBranchAddress("vertices_n",&nPv); 
	tree->SetBranchAddress("muons_n",&nMuons); 
	tree->SetBranchAddress("photons_n",&nPhotons); 
	tree->SetBranchAddress("electrons_n",&nElectrons); 
	tree->SetBranchAddress("met_et",&met);
	tree->SetBranchAddress("met_phi",&met_phi);

	
	tree->SetBranchAddress("electron_px",&electron_px, &b_electron_py); 
	tree->SetBranchAddress("electron_py",&electron_py, &b_electron_py);   
	tree->SetBranchAddress("electron_pz",&electron_pz, &b_electron_pz);
	tree->SetBranchAddress("electron_e",&electron_e, &b_electron_e); 
	
	tree->SetBranchAddress("muon_px",&muon_px, &b_muon_py); 
	tree->SetBranchAddress("muon_py",&muon_py, &b_muon_py);   
	tree->SetBranchAddress("muon_pz",&muon_pz, &b_muon_pz);
	tree->SetBranchAddress("muon_e",&muon_e, &b_muon_e); 
	tree->SetBranchAddress("muon_charge",&muon_charge,&b_muon_charge); 
	
	
	tree->SetBranchAddress("photon_px",&photon_px, &b_photon_py); 
	tree->SetBranchAddress("photon_py",&photon_py, &b_photon_py);   
	tree->SetBranchAddress("photon_pz",&photon_pz, &b_photon_pz);
	tree->SetBranchAddress("photon_e",&photon_e, &b_photon_e); 
	
	tree->SetBranchAddress("jet_px",&jet_px, &b_jet_py); 
	tree->SetBranchAddress("jet_py",&jet_py, &b_jet_py);   
	tree->SetBranchAddress("jet_pz",&jet_pz, &b_jet_pz);
	tree->SetBranchAddress("jet_e",&jet_e, &b_jet_e); 
	tree->SetBranchAddress("jet_bTag",&jet_btag, &b_jet_btag); 

	
	
	
}

void bookHisto(){
	cout << "Book Histograms: " << endl; 
	//float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1600,1800,2000,2200,2500,3000};
	float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	
	CreateHistogram("rapidity","y","y","Events",50,-4,4); 
	CreateHistogram("rapidity_st","y vs S_{T}","S_{T} [GeV]","y",100,0,3000, 50,-4,4); 

	CreateHistogram("st_total","st","S_{T} [GeV]","Events",N,st_bins); 
	CreateHistogram("ht_total","ht","H_{T} [GeV]","Events",N,st_bins); 
	CreateHistogram("nJets","nJets","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu","nJets","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu","nJets","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_W","nJets, #mu=1, 70<M_{T}(#mu,M_{ET})<250 GeV","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_Z","nJets, #mu=2, 75<M_{#mu#mu}<105 GeV","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El","nJets, #mu=1, 1 El","n-jets","Events",10,0.5,10.5); 

	
	CreateHistogram("nJets_stl1000","n-jets S_{T}<1000 GeV","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_stg1000","n-jets S_{T}>1000 GeV","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_st500_1000","n-jets 500<S_{T}<1000 GeV","n-jets","Events",10,0.5,10.5); 

	CreateHistogram("h_met","M_{ET} [GeV]","M_{ET} [GeV]","Events",50,0,800); 
	CreateHistogram("h_npv","pile-up","N_{pv}","Events",50,0,50); 

	CreateProfile("R_dataMC_eff","Ratio #epsilon(data)/#epsilon(MC)","p_{T}(#mu) [GeV]","R(data/MC)",50,0,250); 
	CreateProfile("pr_nJets_npv", "pile-up", "N_{pv}", "n-jets",50,0,50); 
	
	for(int ist=0; ist<N; ist++){
		TString name=Form("nJets_st%.0f-%.0f",st_bins[ist],st_bins[ist+1]);
		TString title=Form("N_{jets}, %.0f<S_{T}<%.0f",st_bins[ist],st_bins[ist+1]);
		CreateHistogram(name,title,"N_{jets}","Events",10,0.5,10.5); 
	}
	
	for(int mu=1; mu<=2; mu++){
			CreateHistogram(Form("nJets_stl1000_%dMu",mu),Form("n-jets S_{T}<1000 GeV, #mu=%d",mu),"n-jets","Events",10,0.5,10.5); 
		CreateHistogram(Form("nJets_stg1000_%dMu",mu),Form("n-jets S_{T}>1000 GeV, #mu=%d",mu),"n-jets","Events",10,0.5,10.5); 
		CreateHistogram(Form("nJets_st500_1000_%dMu",mu),Form("n-jets 500<S_{T}<1000 GeV, #mu=%d",mu),"n-jets","Events",10,0.5,10.5); 
		
		CreateHistogram(Form("nBJets_stl1000_%dMu",mu),Form("n-bjets S_{T}<1000 GeV, #mu=%d",mu),"n-bjets","Events",11,-0.5,10.5); 
		CreateHistogram(Form("nBJets_stg1000_%dMu",mu),Form("n-bjets S_{T}>1000 GeV, #mu=%d",mu),"n-bjets","Events",11,-0.5,10.5); 
		CreateHistogram(Form("nBJets_st500_1000_%dMu",mu),Form("n-bjets 500<S_{T}<1000 GeV, #mu=%d",mu),"n-bjets","Events",11,-0.5,10.5); 
		
	}
	

	for(int nJ=0; nJ<=nJetmax; nJ++){
		
		CreateHistogram(Form("jet%d_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);
		
		CreateHistogram(Form("jet%d_muMass_pt1",nJ),"M(#mu,j), p_{T}(#mu,j)<45","M(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_muMass_pt2",nJ),"M(#mu,j)","M(#mu,j) [GeV], 45<p_{T}(#mu,j)<100 ","Events",100,0,500);
		CreateHistogram(Form("jet%d_muMass_pt3",nJ),"M(#mu,j)","M(#mu,j) [GeV], p_{T}(#mu,j)>100","Events",100,0,500);

		
		CreateHistogram(Form("jet%d_munuMass",nJ),"M(#mu,,#nu,j)","M(#mu,#nu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_nuMass",nJ),"M(#nu,j)","M(#nu,j) [GeV]","Events",100,0,500);
		
		CreateHistogram(Form("jet%d_mujetptM1",nJ),"p_{T}(#mu,j), M(#mu,j)<90","p_{T}(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_mujetptM2",nJ),"p_{T}(#mu,j), M(#mu,j):90-120","p_{T}(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_mujetptM3",nJ),"p_{T}(#mu,j), M(#mu,j)>120","p_{T}(#mu,j) [GeV]","Events",100,0,500);

		
		CreateHistogram(Form("jet%d_ptmu35_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_ptmu40_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_ptmu45_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);

		CreateHistogram(Form("jet%d_ptjet35_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_ptjet40_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);
		CreateHistogram(Form("jet%d_ptjet45_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);
		
		CreateHistogram(Form("jet%d_deltaeta_RF",nJ),"#Delta #eta jet-#mu Rest Frame","#Delta #eta","Events",50,-5,5);
		CreateHistogram(Form("jet%d_deltaphi_RF",nJ),"#Delta #phi jet-#mu Rest Frame","#Delta #phi","Events",50,-3.15,3.15);

		
		CreateHistogram(Form("jet%d_munuMass_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","p_{T}(#mu) [GeV]",100,0,500,100,0,500);

		CreateHistogram(Form("jet%d_mupt_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","p_{T}(#mu) [GeV]",100,0,500,100,0,500);
		CreateHistogram(Form("jet%d_jetpt_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","p_{T}(jet) [GeV]",100,0,500,200,0,1000);

		CreateHistogram(Form("jet%d_met_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","M_{ET} [GeV]",100,0,500,100,0,500);

		CreateHistogram(Form("jet%d_jetmupt_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","p_{T}(#mu,j) [GeV]",100,0,500,100,0,500);
		CreateHistogram(Form("jet%d_jetmuy_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","|y(jet,#mu)|",100,0,500,50,0,2.5);
		CreateHistogram(Form("jet%d_etajet_etamu_muMass",nJ),"#eta(jet),#eta(#mu)","#eta(#mu)","#eta(jet)",50,-2.5,2.5,50,-2.5,2.5);
	
		CreateHistogram(Form("jetrandomy%d_muMass",nJ),"M(#mu,j)","M(#mu,j) [GeV]","Events",100,0,500);

		CreateHistogram(Form("bjet%d_muMass",nJ),"M(#mu,bj)","M(#mu,bj) [GeV]","Events",100,0,500);

		CreateHistogram(Form("rapidity_nJets%d",nJ),"y, n-jets","y","Events",50,-4,4); 
		CreateHistogram(Form("st_nJets%d_uW",nJ),Form("%d jets",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d",nJ),Form("%d jets",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_50GeV",nJ),Form("%d jets",nJ),"S_{T} [GeV]","Events",120,0,6000);
		CreateHistogram(Form("st_variable_pt_nJet%d",nJ),"S_{T}", "S_{T} [GeV]", "Events",N,st_bins); 
		
		CreateHistogram(Form("h_mumu_nJet%d",nJ),"M_{#mu#mu}","Mass [GeV]","Events",100,0,500); 
		
		CreateHistogram(Form("ht_nJets%d",nJ),Form("%d jets",nJ),"H_{T} [GeV]", "Events",N,st_bins);
		for(int pT=30; pT<=100; pT+=10){
			CreateHistogram(Form("st_nJets%d_pTmin%d",nJ, pT),Form("%d jets, p_{T}^{min}=%d",nJ,pT), "S_{T} [GeV]", "Events",N,st_bins); 
			CreateHistogram(Form("st_nJets%d_pTmin%d_uW",nJ, pT),Form("%d jets, p_{T}^{min}=%d",nJ,pT), "S_{T} [GeV]", "Events",N,st_bins); 
			
		}//pT loop 
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_uW",nJ),Form("%d jets, #mu,e",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El",nJ),Form("%d jets, #mu,e",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_W",nJ),Form("%d jets, %d#mu, 70<M_{T}(#mu,M_{ET})<250 GeV",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_W_uW",nJ),Form("%d jets, %d#mu, 70<M_{T}(#mu,M_{ET})<250 GeV",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z",nJ),Form("%d jets, %d#mu, 75<M_{#mu#mu}<105 GeV",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_uW",nJ),Form("%d jets, %d#mu, 75<M_{#mu#mu}<105 GeV",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		
		
		
	}//first nJ loop for inclusive bins (useless!!!) 
	
	for(int mu=1; mu<=2; mu++){		
		for(int nJ=0; nJ<=nJetmax; nJ++){
			
			CreateHistogram(Form("nJets%d_DeltaPhiMuonMet_%dMu",nJ,mu),"#Delta #phi","#Delta #phi(#mu,M_{ET})","Events",20,0,3.15); 
			CreateHistogram(Form("nJets%d_muonpT%d",nJ,mu),"Muon pT", "p_{T}(#mu) [GeV]","Events",300,0,1500); 
			CreateHistogram(Form("jetpT%d_%dMu",nJ,mu),Form("Jet p_{T} spectrum for %dth jet",nJ),"Jet p_{T} [GeV]","Events",300,0,1500);			

			CreateHistogram(Form("st_nJets%d_%dMu_uW",nJ,mu),Form("%d jets, %d#mu",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			CreateHistogram(Form("st_nJets%d_%dMu",nJ,mu),Form("%d jets, %d#mu",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);

			
			CreateHistogram(Form("mt_nJets%d_%dMu",nJ,mu),Form("M_{T} %d-Jets, #mu=%d",nJ,mu),"M_{T} [GeV]","Events",100,0,500);  
			CreateHistogram(Form("met_nJets%d_%dMu",nJ,mu),Form("M_{ET} %d-Jets, #mu=%d",nJ,mu),"M_{ET} [GeV]","Events",200,0,1000);  

			CreateHistogram(Form("st_nJets%d_%dMu_uW_50GeV",nJ,mu),Form("%d jets, %d#mu",nJ,mu),"S_{T} [GeV]","Events",120,0,6000);
			CreateHistogram(Form("st_nJets%d_%dMu_50GeV",nJ,mu),Form("%d jets, %d#mu",nJ,mu),"S_{T} [GeV]","Events",120,0,6000);
			
			
			for(int gamma=0; gamma<=2; gamma++){
				CreateHistogram(Form("st_nJets%d_%dMu_%dg_uW",nJ,mu,gamma),Form("%d jets, %d#mu, %d#gamma",nJ,mu,gamma),"S_{T} [GeV]","Events",N,st_bins);
				CreateHistogram(Form("st_nJets%d_%dMu_%dg",nJ,mu,gamma),Form("%d jets, %d#mu, %d#gamma",nJ,mu,gamma),"S_{T} [GeV]","Events",N,st_bins);
				
			}//gamma loop 
			
			CreateHistogram(Form("DeltaPhi_Jet%d_%dMu",nJ,mu),"","#DeltaR","Events",20,0,6); 
			CreateHistogram(Form("DeltaR_Jet%d_%dMu",nJ,mu),"","#DeltaR","Events",20,0,6); 
			CreateHistogram(Form("DeltaPhiJetMet_Jet%d_%dMu",nJ,mu),"","#Delta #Phi(jet,M_{ET})","Events",20,0,3.15); 
			
		}//nJ loop 
	}//mu loop 
		
	
	
	CreateHistogram("bin_width","bin_width","S_{T} GeV","Events",N,st_bins);
	for (int i=1; i<=hName["bin_width"]->GetNbinsX(); i++) {
		//	cout << hName["bin_width"]->GetBinWidth(i) << " "; 
		hName["bin_width"]->SetBinContent(i,hName["bin_width"]->GetBinWidth(i)); 
	}
	//cout << endl; 
	
	CreateHistogram("Muon_pT","Muon p_{T}","p_{T} GeV","Events",100,0,st_max); 
	CreateHistogram("nJ_sT","nJ_sT","S_{T} GeV", "number of jets", 100, 0, 3000, 6, 0.5, 6.5); 
	
	CreateHistogram("data","data, nJ>=4","S_{T} GeV", "Events", 5700,300,6000);  
	CreateHistogram("eta_sT","eta vs sT", "S_{T} GeV", "|#eta|", 10, 0, 5, 100, 0, 3000); 
	

	CreateHistogram("h_mumu","M_{#mu#mu}","Mass [GeV]","Events",250,0,500); 

	
}

float event_rapidity(){
	TLorentzVector p;
	TLorentzVector Ptot;
	for(int j=0; j<jet_px->size(); j++){
		p.SetPxPyPzE(jet_px->at(j), jet_py->at(j), jet_pz->at(j),jet_e->at(j)); 
		Ptot+=p; 
	}
	for(int j=0; j<muon_px->size(); j++){
		p.SetPxPyPzE(muon_px->at(j), muon_py->at(j), muon_pz->at(j),muon_e->at(j)); 
		Ptot+=p; 
	}
	
	for(int j=0; j<photon_px->size(); j++){
		p.SetPxPyPzE(photon_px->at(j), photon_py->at(j), photon_pz->at(j),photon_e->at(j)); 
		Ptot+=p; 
	}
	for(int j=0; j<electron_px->size(); j++){
		p.SetPxPyPzE(electron_px->at(j), electron_py->at(j), electron_pz->at(j),electron_e->at(j)); 
		Ptot+=p; 
	}
	return Ptot.Rapidity();
	
}



float calcMt(){
	//compute the transverse mass between the leading pt lepton and the met
	
	TLorentzVector pNu;
	pNu.SetPtEtaPhiM(met,0,met_phi,0);
	TLorentzVector pM;
	if(muon_px->size()>=1) {
		pM.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0));
		pM+=pNu;
		mt=pM.Mt(); 
	}

}


float calcSt(){
	float hT=calcHt(); 
	float leptpT=0; 
	float gammapT=0;
	TLorentzVector tmp;
	for(int i=0; i<muon_px->size(); i++){
		tmp.SetPxPyPzE(muon_px->at(i),muon_py->at(i),muon_pz->at(i),muon_e->at(i));
		leptpT+=tmp.Perp();
	}
	
	for(int i=0; i<electron_px->size(); i++){
		tmp.SetPxPyPzE(electron_px->at(i),electron_py->at(i),electron_pz->at(i),electron_e->at(i));
		leptpT+=tmp.Et();
	}
	
	for(int i=0; i<photon_px->size(); i++){
		tmp.SetPxPyPzE(photon_px->at(i),photon_py->at(i),photon_pz->at(i),photon_e->at(i));
		gammapT+=tmp.Et();
	}
	
	float sT=hT+leptpT+gammapT;
	if(met>15)sT+=met;
	
	
	return sT; 
}

float calcHt(){
	float hT=0; 
	for(int j=0; j<jet_px->size(); j++){
		TLorentzVector pJ;
		pJ.SetPxPyPzE(jet_px->at(j), jet_py->at(j), jet_pz->at(j),jet_e->at(j)); 
		hT+=pJ.Perp();
	}//loop over jets
	return hT; 
}



void get_muTrigEff(){

	open_graph("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");
	
	open_graph("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");
}

void open_graph(TString name){
	TGraphAsymmErrors *gr = (TGraphAsymmErrors*)fEff->FindObjectAny(name); 
	grName[name]=gr; 
}

double get_eff(double pt, double eta, TString mode){
	
	double Eff=1; 
	
	if (fabs(eta)<0.9 && mode=="DATA"){
		Eff=grName["IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD"]->Eval(pt); 
	}
	
	if (fabs(eta)>0.9 && fabs(eta)<1.2 && mode=="DATA"){
		Eff=grName["IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD"]->Eval(pt); 
	}
	
	if (fabs(eta)>1.2 && fabs(eta)<2.1 && mode=="DATA"){
		Eff=grName["IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD"]->Eval(pt); 
	}
	
	if (fabs(eta)<0.9 && mode=="MC"){
		Eff=grName["IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD"]->Eval(pt); 
	}
	
	if (fabs(eta)>0.9 && fabs(eta)<1.2 && mode=="MC"){
		Eff=grName["IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD"]->Eval(pt); 
	}
	
	if (fabs(eta)>1.2 && fabs(eta)<2.1 && mode=="MC"){
		Eff=grName["IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD"]->Eval(pt); 
	}
	
	
	return Eff; 
	
}

void efficiency(){
	TLorentzVector pM;
	
	double EffData1=0; 
	double EffMC1=0;
	double EffData2=0;
	double EffMC2=0; 
	double pt=0; 
	
	for(int mu=0; mu<muon_px->size(); mu++){
		pM.SetPxPyPzE(muon_px->at(mu),muon_py->at(mu), muon_pz->at(mu),muon_e->at(mu));
		if(mu==0){
			//eff first muon
			pt=pM.Perp();
		 EffData1=get_eff(pM.Perp(),pM.Eta(),"DATA"); 
		 EffMC1=get_eff(pM.Perp(),pM.Eta(),"MC"); 
		}
		if(mu==1){
			//eff second muon 
			EffData2=get_eff(pM.Perp(),pM.Eta(),"DATA"); 
			EffMC2=get_eff(pM.Perp(),pM.Eta(),"MC"); 
		}
		
	}	//for loop
	//cout << "eff1Data: " << EffData1 << " " << EffData2 << endl; 
	//cout << "eff1MC: " << EffMC1 << " " << EffMC2 << endl; 
	//cout << " data: " << 1-(1-EffData1)*(1-EffData2) << " MC: " << 1-(1-EffMC1)*(1-EffMC2) << endl; 
	
	double R_dataMC=(1-(1-EffData1)*(1-EffData2))/(1-(1-EffMC1)*(1-EffMC2)); 
	//cout << "R_data/MC: " << R_dataMC << endl; 
	if(_MC)prName["R_dataMC_eff"]->Fill(pt,R_dataMC); 
	//cout << "weight: " << weight << endl; 
	if(_MC)weight=weight*R_dataMC; 
	//cout << "weight: " << weight << endl; 
}

void quantile_ratio(TH1F *h_num, TH1F *h_den, TH1F *h_num_uW, TH1F *h_den_uW, TGraphAsymmErrors *g_rat){
	//pass 3 histogrmas, numerator, denominator, and ratio
	
	for (size_t i = 1; i <= h_num->GetNbinsX(); i++) {
		Float_t nNum = h_num_uW->GetBinContent(i);
		Float_t nDen = h_den_uW->GetBinContent(i);
		Float_t nRat = nNum / nDen;
		
		// tail = (1 - CL) / 2; For 95% CL, tail = (1 - 0.95) / 2 = 0.025
		Float_t tail = 0.16;
		Float_t qHi  = 0.0;
		
		if (nDen > 0.0)
			qHi = ROOT::Math::fdistribution_quantile_c(tail, (nNum + 1.0) * 2.0, nDen * 2.0);
		
		Float_t qLo = ROOT::Math::fdistribution_quantile_c(1.0 - tail, nNum * 2.0, (nDen + 1.0) * 2.0);
		
		Float_t rHi = qHi * (nNum + 1.0) / nDen;
		Float_t rLo = 0.0;
		
		if (nDen > 0.0)
			rLo = qLo * nNum / (nDen + 1.0);
		
		if (nDen > 0.0) {
			g_rat->SetPoint     (i - 1, h_num->GetBinCenter(i), nRat);
			g_rat->SetPointError(i - 1, h_num->GetBinWidth (i) / 2.0, h_num->GetBinWidth  (i) / 2.0, nRat - rLo, rHi - nRat);
		}
		//cout << "pT: " << h_rat->GetBinCenter(i) << " N: "  << nNum << " D: " << nDen << endl;
		//if(nDen>0) cout << "pT: " << h_rat->GetBinCenter(i) << " R: " << nRat << " +/-: " << rHi-nRat << " " << nRat-rLo << endl;
	}	
	h_num->Divide(h_den); 
	for(int i=0; i<g_rat->GetN(); i++){
		double x=g_rat->GetX()[i];
		double y=h_num->GetBinContent(h_num->FindBin(x)); 
		
		double exh=g_rat->GetErrorXhigh(i); 
		double exl=g_rat->GetErrorXlow(i); 
		
		double eyh=g_rat->GetErrorYhigh(i); 
		double eyl=g_rat->GetErrorYlow(i); 
		
		double yR=g_rat->GetY()[i]; 
		
		g_rat->SetPoint(i,x,y);
		g_rat->SetPointError(i,exl,exh,(eyl/yR)*y,(eyh/yR)*y); 
		
	}
}

void divide_histo(TString name){

	//Divide histogram by BW
	TH1F *n=(TH1F*)hName[name]->Clone(name+"_GeV");
	
	n->Divide(hName["bin_width"]); 
	n->GetYaxis()->SetTitle("Events/GeV"); 
	hName[n->GetName()]=n; 

}

void divide_BW(){
	
	
	//divide by bin width
	output_file->cd();
	for(int nJ=0; nJ<=nJetmax; nJ++){
		divide_histo(Form("st_nJets%d",nJ));
		divide_histo(Form("ht_nJets%d",nJ)); 
		divide_histo(Form("st_nJets%d_1Mu_1El",nJ)); 
		
		divide_histo(Form("st_nJets%d_W",nJ)); 
		divide_histo(Form("st_nJets%d_Z",nJ)); 

		for(int pT=30; pT<=100; pT+=10){
			divide_histo(Form("st_nJets%d_pTmin%d",nJ, pT));
		}//loop over pT bins
		
		for(int mu=1; mu<=2; mu++){
			divide_histo(Form("st_nJets%d_%dMu_uW",nJ,mu)); 
			divide_histo(Form("st_nJets%d_%dMu",nJ,mu)); 
					
			for(int gamma=0; gamma<=2; gamma++){
				divide_histo(Form("st_nJets%d_%dMu_%dg_uW",nJ,mu,gamma)); 
				divide_histo(Form("st_nJets%d_%dMu_%dg",nJ,mu,gamma)); 
				
			}
		}
		
		}//loop over jets



}


void fill_st(){
	bool print=false;
	
	if(!QCD && muon_px->size()<1) return; 
	//if(met<50) return; 
	//if(photon_px->size()!=2) return; 
	int Nmu=muon_px->size();
	int NEl=electron_px->size();
	int nGamma=photon_px->size();
	
	if(!QCD && Nmu>2) return; 
	//if(nBJets<1) return; 
	if(!QCD && Nmu==2){
		if(muon_charge->at(0)*muon_charge->at(1)==1) return; 
	}
		
	double M_mumu=0; 

	if(!QCD && Nmu==2) {
		TLorentzVector p1; 
		TLorentzVector p2; 
		p1.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0)); 
		p2.SetPxPyPzE(muon_px->at(1),muon_py->at(1),muon_pz->at(1),muon_e->at(1));
		TLorentzVector pmumu;
		pmumu=p1+p2; 
		M_mumu=pmumu.M(); 
		if(pmumu.M()<50) return; 
	}
	if(print)cout << "Fill some general histograms: " << endl; 
	float HT=calcHt();
	//if(!QCD)fill_jets(); 
	hName["h_npv"]->Fill(nPv); 
	prName["pr_nJets_npv"]->Fill(nPv, nJets,weight); 
	if(nJets>=5)hName["h_met"]->Fill(met,weight);
	hName2D["nJ_sT"]->Fill(st,nJets,weight);
	hName["st_total"]->Fill(st,weight);
	hName["ht_total"]->Fill(HT,weight); 
	hName["nJets"]->Fill(nJets,weight);
	
	if(print)cout << "Fill nJets with 1<mu<2 "<< endl; 
	if(st<1000)hName["nJets_stl1000"]->Fill(nJets,weight);
	else hName["nJets_stg1000"]->Fill(nJets,weight); 
	if(st>500 && st<1000)hName["nJets_st500_1000"]->Fill(nJets,weight);

	
	if(print)cout << "nJ loop: " << endl; 
	if(nJets>nJetmax)nJets=nJetmax; 
	if(print) cout << "nJ: " << nJets << " nMu: " << Nmu << endl; 
	if(print) cout << "st: " << st << " mt: " << mt << " ht " << HT << " met " << met << endl; 
	int nJ=nJets;

	if(print) cout << "fill mumu hist: " << endl; 
	if(!QCD) hName[Form("h_mumu_nJet%d",nJ)]->Fill(M_mumu,weight); 
	if(!QCD && Nmu==2 && st>500 && st<1000)hName["h_mumu"]->Fill(M_mumu,weight); 

	//fill 3 control regions for normalization 
	
	//fill Z
	if(Nmu==2 && M_mumu>75 && M_mumu<=105){
		hName["nJets_Z"]->Fill(nJ,weight); 
		hName[Form("st_nJets%d_Z_uW",nJ)]->Fill(st);
		hName[Form("st_nJets%d_Z",nJ)]->Fill(st,weight);
	}
	
	//Fill W
	if(Nmu==1 && mt>70 && mt<250){
		hName["nJets_W"]->Fill(nJ,weight); 
		hName[Form("st_nJets%d_W_uW",nJ)]->Fill(st);
		hName[Form("st_nJets%d_W",nJ)]->Fill(st,weight);
	}
	//fill ttbar
	if(Nmu==1 && NEl==1){
		hName["nJets_1Mu_1El"]->Fill(nJ,weight); 
		hName[Form("st_nJets%d_1Mu_1El_uW",nJ)]->Fill(st); 
		hName[Form("st_nJets%d_1Mu_1El",nJ)]->Fill(st,weight); 
	}
	
	if(!QCD && Nmu==2 && (M_mumu>75 && M_mumu<=105)) return; 
	if(print)cout << "fill main njets histograms: " << endl; 
	hName[Form("nJets_%dMu",Nmu)]->Fill(nJ,weight);
	if(print)cout << "fill main st histograms: " << endl; 

	hName[Form("st_nJets%d_uW",nJ)]->Fill(st); 
	hName[Form("st_nJets%d",nJ)]->Fill(st,weight); 
	hName[Form("ht_nJets%d",nJ)]->Fill(HT,weight);
	hName[Form("st_nJets%d_50GeV",nJ)]->Fill(st,weight); 
	if(QCD)return; 
	if(print)cout << "muon binned: " << endl; 
	
	hName[Form("st_nJets%d_%dMu_uW",nJ,Nmu)]->Fill(st);
	hName[Form("st_nJets%d_%dMu",nJ,Nmu)]->Fill(st,weight);
	
	hName[Form("st_nJets%d_%dMu_uW_50GeV",nJ,Nmu)]->Fill(st);
	hName[Form("st_nJets%d_%dMu_50GeV",nJ,Nmu)]->Fill(st,weight);
	
	
	
	//hName[Form("st_nJets%d_%dMu_%dg_uW",nJ,Nmu,nGamma)]->Fill(st); 
	//hName[Form("st_nJets%d_%dMu_%dg",nJ,Nmu,nGamma)]->Fill(st,weight); 
	if(print) cout << "fill met " << endl; 
	hName[Form("met_nJets%d_%dMu",nJ,Nmu)]->Fill(met,weight); 
	hName[Form("mt_nJets%d_%dMu",nJ,Nmu)]->Fill(mt,weight); 
		
	
	/*
	for(int nJ=nJets; nJ>=1; nJ--){		
		//cout << "nJets: " << nJets << " nJ : " << nJ << endl; 
		hName[Form("st_nJets%d_uW",nJ)]->Fill(st); 
		hName[Form("st_nJets%d",nJ)]->Fill(st,weight); 
		hName[Form("ht_nJets%d",nJ)]->Fill(HT,weight);
		hName[Form("st_nJets%d_50GeV",nJ)]->Fill(st,weight); 
		hName[Form("st_nJets%d_%dMu_uW",nJ,Nmu)]->Fill(st);
		hName[Form("st_nJets%d_%dMu",nJ,Nmu)]->Fill(st,weight);
		hName[Form("st_nJets%d_%dMu_%dg_uW",nJ,Nmu,nGamma)]->Fill(st); 
		hName[Form("st_nJets%d_%dMu_%dg",nJ,Nmu,nGamma)]->Fill(st,weight); 
		
		if(st>500 && st<1000)	hName[Form("h_mumu_nJet%d",nJ)]->Fill(M_mumu,weight); 
		for(int mu=1; mu<Nmu; mu++){
			hName[Form("mt_nJets%d_%dMu",nJ,mu)]->Fill(mt->at(mu),weight); 
		}	
		if(nJ<=nJetGT) break; 

	}
	 
	if(print)cout << "dimuon mass in slice of st: " << endl; 

	
	
	if(nJets>=4){
		hName["st_nJets4p_uW"]->Fill(st);
		hName["st_nJets4p"]->Fill(st,weight);
	}
	else {
		hName["st_nJets23_uW"]->Fill(st);
		hName["st_nJets23"]->Fill(st,weight);
	}

	if(print)cout << "bjets and jets in st bins: " << endl; 

	if(st<1000){
		hName[Form("nJets_stl1000_%dMu",Nmu)]->Fill(nJets,weight);
		hName[Form("nBJets_stl1000_%dMu",Nmu)]->Fill(nBJets,weight);
		
	}
	else {
		hName[Form("nJets_stg1000_%dMu",Nmu)]->Fill(nJets,weight); 
		hName[Form("nBJets_stg1000_%dMu",Nmu)]->Fill(nBJets,weight); 
		
	}
	if(st>500 && st<1000){
		hName[Form("nJets_st500_1000_%dMu",Nmu)]->Fill(nJets,weight);
		hName[Form("nBJets_st500_1000_%dMu",Nmu)]->Fill(nBJets,weight);
		
	}
*/
	
}

void fill_bjet_mu_mass(){
	//Form("jet%d_muMass",nJ)
	TLorentzVector pJ;
	TLorentzVector pM;
	
	TLorentzVector pMRF;
	TLorentzVector pJRF;
	
	if(muon_px->size()!=1) return; 

	int imu=0; 
	pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
	pMRF=pM;
	
	TLorentzVector pNu;
	pNu.SetPtEtaPhiM(met,0,met_phi,0);
	
	for(int ijet=0; ijet<jet_px->size(); ijet++){
		int jet=ijet+1;
		if(jet>nJetmax)jet=nJetmax; 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		if(jet_btag->at(ijet)){
			//fill histogram for b,mu mass 
			pJ+=pM; 
			double BMuMass=pJ.M();
			hName[Form("bjet%d_muMass",jet)]->Fill(BMuMass,weight);
		}
		//now fill histograms for any flavor jet,mu mass 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		pJRF=pJ; 
		double ptjet=pJ.Perp();
		double jetEta=pJ.Eta();
		pJ+=pM; 
		pMRF.Boost(pJ.BoostVector()); 
		pJRF.Boost(pJ.BoostVector()); 
		
		double jetMuMass=pJ.M();
		double jetMuPt=pJ.Perp();
		double jetMuy=fabs(pJ.Rapidity()); 
		
		pJ+=pNu;
		double	 jetMuNuMass=pJ.M(); 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		pJ+=pNu;
		double jetNuMass=pJ.M();
		
		
		hName[Form("jet%d_muMass",jet)]->Fill(jetMuMass,weight);

		if(jetMuPt<45) hName[Form("jet%d_muMass_pt1",jet)]->Fill(jetMuMass,weight);
		if(jetMuPt>45 && jetMuPt<100) hName[Form("jet%d_muMass_pt2",jet)]->Fill(jetMuMass,weight);
		if(jetMuPt>100)hName[Form("jet%d_muMass_pt3",jet)]->Fill(jetMuMass,weight);

		hName[Form("jet%d_munuMass",jet)]->Fill(jetMuNuMass,weight);
		hName[Form("jet%d_nuMass",jet)]->Fill(jetNuMass,weight);


		if(pM.Perp()>35)	hName[Form("jet%d_ptmu35_muMass",jet)]->Fill(jetMuMass,weight);
		if(pM.Perp()>40)	hName[Form("jet%d_ptmu40_muMass",jet)]->Fill(jetMuMass,weight);
		if(pM.Perp()>45)	hName[Form("jet%d_ptmu45_muMass",jet)]->Fill(jetMuMass,weight);

		
		if(ptjet>35)	hName[Form("jet%d_ptjet35_muMass",jet)]->Fill(jetMuMass,weight);
		if(ptjet>40)	hName[Form("jet%d_ptjet40_muMass",jet)]->Fill(jetMuMass,weight);
		if(ptjet>45)	hName[Form("jet%d_ptjet45_muMass",jet)]->Fill(jetMuMass,weight);
		
		if(jetMuMass<=90) hName[Form("jet%d_mujetptM1",jet)]->Fill(jetMuPt,weight);
		if(jetMuMass>90 && jetMuMass<120) hName[Form("jet%d_mujetptM2",jet)]->Fill(jetMuPt,weight);
		if(jetMuMass>=120) hName[Form("jet%d_mujetptM3",jet)]->Fill(jetMuPt,weight);

		hName2D[Form("jet%d_munuMass_muMass",jet)]->Fill(jetMuMass,jetMuNuMass,weight);
		hName2D[Form("jet%d_mupt_muMass",jet)]->Fill(jetMuMass,pM.Perp(),weight);
		hName2D[Form("jet%d_jetpt_muMass",jet)]->Fill(jetMuMass,ptjet,weight);
		hName2D[Form("jet%d_met_muMass",jet)]->Fill(jetMuMass,met,weight);
		hName2D[Form("jet%d_jetmupt_muMass",jet)]->Fill(jetMuMass,jetMuPt,weight);
		hName2D[Form("jet%d_jetmuy_muMass",jet)]->Fill(jetMuMass,jetMuy,weight);
		hName2D[Form("jet%d_etajet_etamu_muMass",jet)]->Fill(pM.Eta(),jetEta,weight);

		hName[Form("jet%d_deltaeta_RF",jet)]->Fill(pJRF.Eta()-pMRF.Eta(),weight); 
		hName[Form("jet%d_deltaphi_RF",jet)]->Fill(pJRF.DeltaPhi(pMRF),weight); 
		
		pJ.SetPxPyPzE(jet_px->at(ijet),-jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		pJ+=pM; 
		jetMuMass=pJ.M();
		hName[Form("jetrandomy%d_muMass",jet)]->Fill(jetMuMass,weight);

		
	}	
}

void nBtag(){
	int nB=0; 
	for(int ib=0; ib<jet_btag->size();ib++){
		bool bT=jet_btag->at(ib);
		if(bT==1) nB++;
	}
	nBJets=nB; 
}

void fill_rapidity(){
	float y=event_rapidity();
	hName["rapidity"]->Fill(y,weight); 
	hName2D["rapidity_st"]->Fill(st,y,weight);
	if(nJets<nJetmax ) hName[Form("rapidity_nJets%d",nJets)]->Fill(y,weight); 	
}

void fill_muon(){
	TLorentzVector pM; 
	for(int imu=0; imu<muon_px->size(); imu++){
		int mu=imu+1;
		if(mu>2) continue; 
		pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
		hName[Form("nJets%d_muonpT%d",nJets,mu)]->Fill(pM.Perp(),weight); 
	}
}

void fill_jets(){
	TLorentzVector pJ;
	for(int ijet=0; ijet<jet_px->size(); ijet++){
		int jet=ijet+1;
		if(jet>nJetmax)jet=nJetmax; 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		hName[Form("jetpT%d_%dMu",nJets,muon_px->size())]->Fill(pJ.Perp(),weight); 
	}
}

void fill_deltaR(){
	TLorentzVector pJ; 
	TLorentzVector pM; 
	for(int ijet=0; ijet<jet_px->size(); ijet++){
		int jet=ijet+1;
		if(jet>nJetmax)jet=nJetmax; 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		for(int imu=0; imu<muon_px->size(); imu++){
			pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
			double dR=pJ.DeltaR(pM); 
			double dPhi=pJ.DeltaPhi(pM);
			int mu=imu+1;
			if(mu>2) continue; 
			hName[Form("DeltaR_Jet%d_%dMu",jet,mu)]->Fill(dR,weight); 
			hName[Form("DeltaPhi_Jet%d_%dMu",jet,mu)]->Fill(dPhi,weight); 

		}
	}
}

void fill_deltaPhi(){
	TLorentzVector pJ; 
	TLorentzVector pMet;
	TLorentzVector pM;
	pMet.SetPtEtaPhiM(met,0,met_phi,0);

	int Nmu=muon_px->size();
	
	for(int ijet=0; ijet<jet_px->size(); ijet++){
		int jet=ijet+1;
		if(jet>nJetmax)jet=nJetmax; 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		double dPhi=TMath::Abs(pJ.DeltaPhi(pMet));
		hName[Form("DeltaPhiJetMet_Jet%d_%dMu",jet,Nmu)]->Fill(dPhi,weight); 

	}
	
	for(int imu=0; imu<muon_px->size(); imu++){
		pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
		double dPhi=TMath::Abs(pM.DeltaPhi(pMet)); 
		int mu=imu+1;
		if(mu>2)continue;
		hName[Form("nJets%d_DeltaPhiMuonMet_%dMu",nJets,mu)]->Fill(dPhi,weight); 

	}
	
}

void remove_muon(int i){
	muon_px->erase(muon_px->begin()+i);
	muon_py->erase(muon_py->begin()+i);
	muon_pz->erase(muon_pz->begin()+i);
	muon_e->erase(muon_e->begin()+i);
}

void remove_jet(int ijet){
	jet_px->erase( jet_px->begin() + ijet );
	jet_py->erase( jet_py->begin() + ijet );
	jet_pz->erase( jet_pz->begin() + ijet );
	jet_e->erase( jet_e->begin() + ijet );
}

void print_jets(){
	TLorentzVector pJ;
	for(int ijet=0; ijet<jet_px->size(); ijet++){
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
		cout << "pT: " << pJ.Perp() << " eta: " << pJ.Eta() << endl; 
	}
}

void print_muons(){
	TLorentzVector pM;
	for(int imu=0; imu<muon_px->size(); imu++){
		pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
		cout << "pT: " << pM.Perp() << " eta: " << pM.Eta() << endl; 
	}
}

bool cuts(){
	bool print=0; 
	
	int NEl=electron_px->size(); 
	int Ngamma=photon_px->size(); 
	
	double jetCut_pt1=125; 
	double jetCut_pt2=45;
	double MetCut=55;
	double DeltaPhiCut_jetMet=0.5; 
	double DeltaPhiCut_MetMuon=0.8;
	double muonPtCut=45; 
	TLorentzVector pMet;
	pMet.SetPtEtaPhiM(met,0,met_phi,0);
	
	
	//if(met<MetCut) return false; 
	

	TLorentzVector pJ;
	TLorentzVector pM;
	
	
	if(print)cout << "before muon cuts: "<< endl; 
	if(print)print_muons();
	
	std::vector<float>::size_type imu = 0;
	
	while ( imu < muon_px->size() ) {
		pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
		if(pM.Perp()<muonPtCut) remove_muon(imu);
		else {
			++imu;
		}
	}
	imu=0; 
	
	while ( imu < muon_px->size() ) {
		pM.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu)); 
		if(TMath::Abs(pM.DeltaPhi(pMet))<DeltaPhiCut_MetMuon) remove_muon(imu); 
		else {
			++imu;
		}
	}
	
	if(print)cout << "after muon cuts: " << endl; 
/*
	if(print)cout << "jets before cuts: " << endl; 
	if(print)print_jets();
	
	if(nJets==1) {
		pJ.SetPxPyPzE(jet_px->at(0),jet_py->at(0),jet_pz->at(0),jet_e->at(0)); 

		if(pJ.Perp()<jetCut_pt1){
			remove_jet(0);
		}
	}
	
	if(nJets>=2) {
		pJ.SetPxPyPzE(jet_px->at(0),jet_py->at(0),jet_pz->at(0),jet_e->at(0)); 
		
		if(pJ.Perp()<jetCut_pt1){
			remove_jet(0);
		}
		
		pJ.SetPxPyPzE(jet_px->at(0),jet_py->at(0),jet_pz->at(0),jet_e->at(0)); 
		
		if(pJ.Perp()<jetCut_pt2){
			remove_jet(0);
		}
	}
	*/
	std::vector<float>::size_type ijet = 0;
	
	while ( ijet < jet_px->size() ) {
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
		if(TMath::Abs(pJ.DeltaPhi(pMet))<DeltaPhiCut_jetMet) remove_jet(ijet);		
		else {
			++ijet;
		}
	}
	if(print)cout << "after delta Phi cut: " << endl; 
	if(print)print_jets();
	
	nJets=jet_px->size();
	
	int Nmu=muon_px->size(); 

	if(Nmu<1 || Nmu>2) return false; 

	st=calcSt();
	if(st<300) return false; 
	
}

void event_loop(TChain *tree){
	int maxEvt=10000;
	//entries=maxEvt;
	
	TStopwatch t; 
	t.Start(kFALSE); 
	
	double time_evt=0; 
	double W=weight; 
	cout << "Entries: " << entries << endl; 
	for (int iEvt=0; iEvt<entries; iEvt++){
		if(iEvt%1000000==0)cout << iEvt << endl;
		if(iEvt%1000000==0){
			time_evt=t.RealTime()/static_cast<double> (iEvt);
			t.Continue();
			cout << "Real Time: " << t.RealTime() << endl; 
			cout << "Time/Evt: " << time_evt << endl; 
			cout << "Predicted completion time: " << time_evt*static_cast<double> ((entries-iEvt))/60 << " minutes" << endl;
		}
		tree->GetEntry(iEvt);
		if(muon_px->size()<1 || muon_px->size()>2) continue; 

		nBtag(); 
		weight=W; 
		efficiency(); 
		calcMt();
		if(nJets>nJetmax)nJets=nJetmax;
		fill_bjet_mu_mass();
	
		fill_deltaR(); 
		fill_deltaPhi();
		fill_rapidity();
		if(st<300) continue; 
		fill_jets();
		fill_muon(); 
		//if(!cuts()) continue; 

		fill_st();
//		cout << "object loops: " << endl; 
//		object_loops(15, 30);
		
	}
	
	
}

void analyze_file(TString folder){
	get_muTrigEff();
	//TTree *tree = (TTree*)infile->Get("tree");
	TChain *tree = new TChain("tree"); 
	TString FN=folder+"/hist_analysis*.root";
	cout << "FN: " << FN << endl; 
	tree->Add(FN);
	cout << "Intitialize Tree: " << endl; 
	initialize_tree(tree);
	event_loop(tree);
	delete tree; 
}

void initialize(){
	cout <<"initialize: " << endl; 
	for(int nJ=1; nJ<=4; nJ++){
		fileName[Form("w%dJets",nJ)]=Form("w%dJetsToLNu",nJ);
		fileName[Form("dy%dJets",nJ)]=Form("dy%dJetsToLL",nJ); 
	}

	
	xs["w1Jets"]=6662;
	xs["w2Jets"]=2159;
	xs["w3Jets"]=640;
	xs["w4Jets"]=264;
	
	if(LO){
		xs["w1Jets"]=5400;
		xs["w2Jets"]=1750;
		xs["w3Jets"]=519;
		xs["w4Jets"]=214;
	}
	
	Ngen["w1Jets"]=23141598;
	Ngen["w2Jets"]=34044921;
	Ngen["w3Jets"]=15539503;
	Ngen["w4Jets"]=13382803;
	
	xs["wJetsHT150"]=290.7;
	xs["wJetsHT200"]=111.37;
	xs["wJetsHT250"]=59.24;
	xs["wJetsHT300"]=47.25;
	xs["wJetsHT400"]=31.12;
	
	
	if(LO){
		xs["wJetsHT150"]=235.6;
		xs["wJetsHT200"]=90.27;
		xs["wJetsHT250"]=48.01;
		xs["wJetsHT300"]=38.3;
		xs["wJetsHT400"]=25.22;
	}
	
	Ngen["wJetsHT150"]=21686209;
	Ngen["wJetsHT200"]=10039771;
	Ngen["wJetsHT250"]=6575572;
	Ngen["wJetsHT300"]=6840509;
	Ngen["wJetsHT400"]=6619654;
	
	fileName["wJetsHT150"]="wJets_HT_150To200";
	fileName["wJetsHT200"]="wJets_HT_200To250";
	fileName["wJetsHT250"]="wJets_HT_250To300";
	fileName["wJetsHT300"]="wJets_HT_300To400";
	fileName["wJetsHT400"]="wJets_HT_400Toinf";
	
	xs["dy1Jets"]=666;
	xs["dy2Jets"]=215;
	xs["dy3Jets"]=61;
	xs["dy4Jets"]=27;
	
	Ngen["dy1Jets"]=24045248;
	Ngen["dy2Jets"]=21852156;
	Ngen["dy3Jets"]=11015445;
	Ngen["dy4Jets"]=6402827;
	
	xs["ttSemiLept"]=107.6;
	xs["ttFullLept"]=25.6;
	
	fileName["ttSemiLept"]="ttJetsSemiLept";
	fileName["ttFullLept"]="ttJetsFullLept";

	Ngen["ttSemiLept"]=86711159;
	Ngen["ttFullLept"]=12119013;
	
	
	xs["ttG"]=25.24;
	Ngen["ttG"]=1719954;
	fileName["ttG"]="ttGJets";
	
	
	//single top xs
	
	fileName["TBar_t"]="TBar_t";
	fileName["TBar_s"]="TBar_s";
	fileName["TBar_tW"]="TBar_tW";
	
	fileName["T_t"]="T_t";
	fileName["T_s"]="T_t";
	fileName["T_tW"]="T_t"; 
	
	xs["TBar_t"]=30.7;
	xs["TBar_s"]=1.76; 
	xs["TBar_tW"]=11.1;
	
	xs["T_t"]=56.4;
	xs["T_s"]=3.79; 
	xs["T_tW"]=11.1; 
	
	Ngen["TBar_t"]=1935072;
	Ngen["TBar_s"]=139974; 
	Ngen["TBar_tW"]=493460;
	
	Ngen["T_t"]=3758227;
	Ngen["T_s"]=259961; 
	Ngen["T_tW"]=497658;
	
	
	xs["wJets_inclusive"]=37509; 
	Ngen["wJets_inclusive"]=18393090; 
	fileName["wJets_inclusive"]="wJetsToLNu";
	
	fileName["singleMu"]="singleMu";
	
	fileName["RPV300"]="RPV_300";
	xs["RPV300"]=1.99608; 
	Ngen["RPV300"]=39198;
	
	fileName["RPV400"]="RPV_400";
	xs["RPV400"]=0.35683; 
	Ngen["RPV400"]=39999;
	
	fileName["RPV500"]="RPV_500";
	xs["RPV500"]=0.0855847; 
	Ngen["RPV500"]=34397;
	
	fileName["RPV600"]="RPV_600";
	xs["RPV600"]=0.0248009; 
	Ngen["RPV600"]=39193;
	
	fileName["RPV700"]="RPV_700";
	xs["RPV700"]=0.0081141; 
	Ngen["RPV700"]=39998;

	fileName["RPV800"]="RPV_800";
	xs["RPV800"]=0.00289588; 
	Ngen["RPV800"]=39998;
	
	fileName["RPV900"]="RPV_900";
	xs["RPV900"]=0.00109501; 
	Ngen["RPV900"]=39997;
	
	fileName["RPV1000"]="RPV_1000";
	xs["RPV1000"]=0.000435488; 
	Ngen["RPV1000"]=39996;
	
	fileName["stealth_ww1400_300"]="stealth_ww1400_300"; 
	xs["stealth_ww1400_300"]=0.00014128; 
	Ngen["stealth_ww1400_300"]=743/0.36; 
	
	fileName["stealth_ww500_300"]="stealth_ww500_300"; 
	xs["stealth_ww500_300"]=0.847051; 
	Ngen["stealth_ww500_300"]=777/0.36; 
	
	fileName["QCD_30-50"]="QCD_30-50"; 
	fileName["QCD_50-80"]="QCD_50-80"; 
	fileName["QCD_80-120"]="QCD_80-120"; 
	fileName["QCD_120-170"]="QCD_120-170"; 	
	fileName["QCD_170-300"]="QCD_170-300"; 
	fileName["QCD_300-470"]="QCD_300-470"; 
	fileName["QCD_470-600"]="QCD_470-600"; 
	fileName["QCD_600-800"]="QCD_600-800"; 
	fileName["QCD_800-1000"]="QCD_800-1000"; 

	
	xs["QCD_30-50"]=806298; //xs*filtereff
	xs["QCD_50-80"]=176187.6;  
	xs["QCD_80-120"]=40448; 
	xs["QCD_120-170"]=7463.9;
	xs["QCD_170-300"]=2299.75;
	xs["QCD_300-470"]=151;
	xs["QCD_470-600"]=11.796;
	xs["QCD_600-800"]=0.267;
	xs["QCD_800-1000"]=0.36;
	
	Ngen["QCD_30-50"]=9560265;
	Ngen["QCD_50-80"]=10365230;
	Ngen["QCD_80-120"]=9238642;
	Ngen["QCD_120-170"]=8501935;
	Ngen["QCD_170-300"]=7669947;
	Ngen["QCD_300-470"]=7832261;
	Ngen["QCD_470-600"]=3783069;
	Ngen["QCD_600-800"]=4119000;
	Ngen["QCD_800-1000"]=4107853;
	
	fileName["QCD4Jet_Pt100-180"]="QCD_4Jets_Pt100-180"; 
	fileName["QCD4Jet_Pt250-400"]="QCD_4Jets_Pt250-400"; 
	fileName["QCD4Jet_Pt400-5600"]="QCD_4Jets_Pt400-5600"; 
	
	fileName["QCD6Jet_Pt100-180"]="QCD_6Jets_Pt100-180"; 
	fileName["QCD6Jet_Pt180-250"]="QCD_6Jets_Pt180-250"; 
	fileName["QCD6Jet_Pt250-400"]="QCD_6Jets_Pt250-400"; 
	fileName["QCD6Jet_Pt400-5600"]="QCD_6Jets_Pt400-5600"; 
	
	
	xs["QCD4Jet_Pt100-180"]=141163.0;
	xs["QCD4Jet_Pt250-400"]=4480.4;
	xs["QCD4Jet_Pt400-5600"]=469.15;
	
	xs["QCD6Jet_Pt100-180"]=9535.76;
	xs["QCD6Jet_Pt180-250"]=1915.08;
	xs["QCD6Jet_Pt250-400"]=758.0;
	xs["QCD6Jet_Pt400-5600"]=102.467;
	
	Ngen["QCD4Jet_Pt100-180"]=141163.0;
	Ngen["QCD4Jet_Pt250-400"]=4480.4;
	Ngen["QCD4Jet_Pt400-5600"]=469.15;
	
	Ngen["QCD6Jet_Pt100-180"]=9535.76;
	Ngen["QCD6Jet_Pt180-250"]=1915.08;
	Ngen["QCD6Jet_Pt250-400"]=758.0;
	Ngen["QCD6Jet_Pt400-5600"]=102.467;
	
	
	fileName["QCD_inc"]="QCD_Pt15-3000"; 
	xs["QCD_inc"]=29981599700;
	Ngen["QCD_inc"]=9620046; 
	
	fileName["WWJetsTo2L2Nu"]="WWJetsTo2L2Nu";
	fileName["WZJetsTo2L2Q"]="WZJetsTo2L2Q";
	fileName["WZJetsTo3LNu"]="WZJetsTo3LNu";
	fileName["WZJetsTo2QLNu"]="WZJetsTo2QLNu";
	
	xs["WWJetsTo2L2Nu"]=4.7;
	xs["WZJetsTo2L2Q"]=1.755;
	xs["WZJetsTo3LNu"]=0.8674; 
	xs["WZJetsTo2QLNu"]=3.1;
	
	Ngen["WWJetsTo2L2Nu"]=1933235;
	Ngen["WZJetsTo2L2Q"]=3215990;
	Ngen["WZJetsTo3LNu"]=1947979;
	Ngen["WZJetsTo2QLNu"]=2908657;
	
	
	fileName["ZZJetsTo2L2Q"]="ZZJetsTo2L2Q";
	fileName["ZZJetsTo2L2Nu"]="ZZJetsTo2L2Nu";
	fileName["ZZJetsTo4L"]="ZZJetsTo4L";

	xs["ZZJetsTo2L2Q"]=0.91;
	xs["ZZJetsTo2L2Nu"]=0.28;
	xs["ZZJetsTo4L"]=0.1296;
	
	Ngen["ZZJetsTo2L2Q"]=1936727;
	Ngen["ZZJetsTo2L2Nu"]=954911;
	Ngen["ZZJetsTo4L"]=4807893;
	
	fileName["ggjets"]="gg4jets";
	xs["ggjets"]=0.83; 
	Ngen["ggjets"]=6000; 
	
	xs["singleMu"]=1.0;
	Ngen["singleMu"]=1.0; 
}

void setup_files(int jobnumber){

	
	if(jobnumber==0) {
		sample_list.push_back("singleMu"); 
		_MC=false; 
		output_file_name="output_file_singleMu.root"; 
		selections="_trigger";
	}
	if(jobnumber==1) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("w%dJets",nJ);
			sample_list.push_back(tag); 
		}
		_MC=true; 
		output_file_name="output_file_wJets.root"; 
		selections="_trigger";
	}
	if(jobnumber==2) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag); 
		}
		_MC=true; 
		output_file_name="output_file_dy.root"; 
		selections="_trigger";

	}
	
	if(jobnumber==3) {
		sample_list.push_back("ttSemiLept"); 
		sample_list.push_back("ttFullLept"); 
		//sample_list.push_back("ttG"); 
		
		_MC=true; 
		output_file_name="output_file_ttbar.root"; 
		selections="_trigger";

	}
	
	if(jobnumber==4) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("w%dJets",nJ);
			sample_list.push_back(tag); 
		}
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag); 
		}
		sample_list.push_back("ttSemiLept"); 
		sample_list.push_back("ttFullLept"); 
		
		sample_list.push_back("QCD_30-50"); 
		sample_list.push_back("QCD_50-80");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_120-170");
		sample_list.push_back("QCD_170-300");
		sample_list.push_back("QCD_300-470");
		sample_list.push_back("QCD_470-600");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_600-800");
		
		sample_list.push_back("TBar_t");
		sample_list.push_back("TBar_s");
		sample_list.push_back("TBar_tW");
		
		sample_list.push_back("T_t");
		sample_list.push_back("T_s");
		sample_list.push_back("T_tW");
		
		sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		sample_list.push_back("WZJetsTo2QLNu");
		
		sample_list.push_back("ZZJetsTo2L2Q");
		sample_list.push_back("ZZJetsTo2L2Nu");
		sample_list.push_back("ZZJetsTo4L");
		
		//sample_list.push_back("ttG"); 
		
		_MC=true; 
		output_file_name="output_file_allMC.root"; 
		selections="_trigger";

	}
	
	if(jobnumber==5) {
		sample_list.push_back("ttSemiLept"); 
			
		_MC=true; 
		output_file_name="output_file_ttSemiLept.root"; 
		selections="_trigger";

	}
	if(jobnumber==6) {
		sample_list.push_back("ttFullLept"); 
		
		_MC=true; 
		output_file_name="output_file_ttFullLept.root"; 
		selections="_trigger";

	}
	if(jobnumber==7) {
		sample_list.push_back("QCD_30-50"); 
		sample_list.push_back("QCD_50-80");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_120-170");
		sample_list.push_back("QCD_170-300");
		sample_list.push_back("QCD_300-470");
		sample_list.push_back("QCD_470-600");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_600-800");
		sample_list.push_back("QCD_800-1000");

		_MC=true; 
		output_file_name="output_file_QCD.root"; 
		selections="_trigger";
		
	}
	if(jobnumber==8) {
		sample_list.push_back("RPV300"); 
		
		_MC=true; 
		output_file_name="output_file_RPV300.root"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==9) {
		sample_list.push_back("RPV400"); 
		
		_MC=true; 
		output_file_name="output_file_RPV400.root"; 
		selections="_trigger";
		
	}
	
	
	if(jobnumber==10) {
		sample_list.push_back("RPV500"); 
		
		_MC=true; 
		output_file_name="output_file_RPV500.root"; 
		selections="_trigger";

	}
	if(jobnumber==11) {
		sample_list.push_back("RPV600"); 
		
		_MC=true; 
		output_file_name="output_file_RPV600.root"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==12) {
		sample_list.push_back("RPV700"); 
		
		_MC=true; 
		output_file_name="output_file_RPV700.root"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==13) {
		sample_list.push_back("RPV800"); 
		
		_MC=true; 
		output_file_name="output_file_RPV800.root"; 
		selections="_trigger";
		
	}
	
	
	if(jobnumber==14) {
		sample_list.push_back("RPV900"); 
		
		_MC=true; 
		output_file_name="output_file_RPV900.root"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==15) {
		sample_list.push_back("RPV1000"); 
		
		_MC=true; 
		output_file_name="output_file_RPV1000.root"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==16) {
		sample_list.push_back("stealth_ww1400_300"); 
		
		_MC=true; 
		output_file_name="output_file_stealth_ww1400_300.root"; 
		selections="_trigger";
		
	}
	if(jobnumber==17) {
		sample_list.push_back("stealth_ww500_300"); 
		
		_MC=true; 
		output_file_name="output_file_stealth_ww500_300.root"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==18){
		sample_list.push_back("TBar_t");
		sample_list.push_back("TBar_s");
		sample_list.push_back("TBar_tW");
		
		sample_list.push_back("T_t");
		sample_list.push_back("T_s");
		sample_list.push_back("T_tW");
		_MC=true; 
		output_file_name = "output_file_singleTop.root"; 
		selections="_trigger"; 

	}
	/*
	if(jobnumber==19){	
		
		sample_list.push_back("QCD4Jet_Pt100-180"); 
		sample_list.push_back("QCD4Jet_Pt250-400"); 
		sample_list.push_back("QCD4Jet_Pt400-5600"); 

		sample_list.push_back("QCD6Jet_Pt100-180"); 
		sample_list.push_back("QCD6Jet_Pt180-250"); 
		sample_list.push_back("QCD6Jet_Pt250-400"); 
		sample_list.push_back("QCD6Jet_Pt400-5600"); 
		
		_MC=true; 
		output_file_name = "output_file_QCD4Jets6Jets.root"; 
		selections=""; 
		
	}

	
	if(jobnumber==20){	
		
		sample_list.push_back("QCD_inc"); 
		_MC=true; 
		QCD=true;
		output_file_name = "output_file_QCDincl.root"; 
		selections=""; 
		
	}
	*/
	
	if(jobnumber==19){	
		sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		sample_list.push_back("WZJetsTo2QLNu");
		
		sample_list.push_back("ZZJetsTo2L2Q");
		sample_list.push_back("ZZJetsTo2L2Nu");
		sample_list.push_back("ZZJetsTo4L");
		
		_MC=true; 
		QCD=false;
		output_file_name = "output_file_diboson.root"; 
		selections="_trigger"; 
		
	}
	
	if(jobnumber==20){
		sample_list.push_back("wJets_inclusive");
		_MC=true;
		QCD=false; 
		output_file_name="output_file_wJets_inclusive.root";
		selections="_trigger"; 
		
	}
	
}

void open_files(){
	TString dir="/eos/uscms/store/user/btcarlso/trees/";

	initialize();
		
	//TString sample_list[]="wJetsHT150","wJetsHT200","wJetsHT250","wJetsHT300","wJetsHT400",
	//TString sample_list[]={"RPV500"};
	
	int N_files=sample_list.size(); 
	cout << "N Files: " << N_files << endl; 

	weight=1;
	cout << "output file: " << output_file_name << endl; 
	output_file=new TFile(output_file_name,"RECREATE"); 
	cout << "Analyzing samples: " << endl; 
	for (int iF=0; iF<N_files; iF++){
		cout << sample_list.at(iF) << " " << fileName[sample_list.at(iF)] << endl; 
	}
	cout << endl; 
	
	for (int iF=0; iF<N_files; iF++) {
		cout << sample_list[iF] << endl; 
		if(_MC) weight=xs[sample_list[iF]]*19600/(Ngen[sample_list[iF]]);
		cout << "weight: " << weight << endl; 
		TString FN=dir+fileName[sample_list.at(iF)]+selections;

		analyze_file(FN);
	
	}
	
}

//work functions, eg create histograms

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
					 Int_t       nBinsX, const Float_t* xBins)
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

void CreateProfile(const char* name,   const char* title,
									const char* xTitle, const char* yTitle,
									Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TProfile* pr = new TProfile(name, title, nBinsX, xLow, xUp);
	
	pr->GetXaxis()->SetTitle(xTitle);
	pr->GetYaxis()->SetTitle(yTitle);
	
	prName[name] = pr;
}

void writeHisto(){
	output_file->cd();
	int N1D=0, N2D=0, NP1D=0; 
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++)N1D++;
	//	for (std::map<TString,TProfile*>::iterator it=prName.begin(); it!=prName.end(); it++) NP1D++;
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++)N2D++; 
	
	//	CreateHistogram("TProfile_names","","","",NP1D,0.,10.); 

	CreateHistogram("TH1F_names","","","",N1D,0.,10.); 
	CreateHistogram("TH2F_names","","","",N2D,0.,10.); 
	
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		if(it->first=="TH1F_names" || it->first=="TH2F_names")continue;
		hName["TH1F_names"]->Fill(it->first,1);
	}
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++){
		if(it->first=="TH1F_names" || it->first=="TH2F_names")continue;
		hName["TH2F_names"]->Fill(it->first,1); 
	}
	
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		hName[it->first]->Write();
	}
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++) {
		hName2D[it->first]->Write();
	}

	for (std::map<TString,TProfile*>::iterator it=prName.begin(); it!=prName.end(); it++){
		prName[it->first]->Write();
	}

	output_file->Close();
	delete output_file; 
}