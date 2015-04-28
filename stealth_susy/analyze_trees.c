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
	if(MakeEfficiency) bookEff();
	open_files();
	
	divide_BW();
	writeHisto();
	gSystem->Exit(1,1); 
}

void initialize_tree(TChain *tree){
	entries=tree->GetEntries();
	
	tree->SetBranchAddress("jets_n",&nJets); 
	
	tree->SetBranchAddress("st",&st); 
	tree->SetBranchAddress("gen_st",&gen_st); 
	tree->SetBranchAddress("gen_ht",&gen_ht); 

	tree->SetBranchAddress("vertices_n",&nPv); 
	tree->SetBranchAddress("muons_n",&nMuons); 
	tree->SetBranchAddress("photons_n",&nPhotons); 
	tree->SetBranchAddress("electrons_n",&nElectrons); 
	tree->SetBranchAddress("met_et",&met);
	tree->SetBranchAddress("met_phi",&met_phi);

	tree->SetBranchAddress("puWeight_nom",&puWt_nom); 
	tree->SetBranchAddress("puWeight_up",&puWt_up); 
	tree->SetBranchAddress("puWeight_down",&puWt_down); 
	
	tree->SetBranchAddress("electron_px",&electron_px, &b_electron_px); 
	tree->SetBranchAddress("electron_py",&electron_py, &b_electron_py);   
	tree->SetBranchAddress("electron_pz",&electron_pz, &b_electron_pz);
	tree->SetBranchAddress("electron_e",&electron_e, &b_electron_e); 
	tree->SetBranchAddress("electron_charge",&electron_charge, &b_electron_charge); 
	
	//loose electrons
	tree->SetBranchAddress("loose_electron_px",&loose_electron_px, &b_loose_electron_px); 
	tree->SetBranchAddress("loose_electron_py",&loose_electron_py, &b_loose_electron_py);   
	tree->SetBranchAddress("loose_electron_pz",&loose_electron_pz, &b_loose_electron_pz);
	tree->SetBranchAddress("loose_electron_e",&loose_electron_e, &b_loose_electron_e); 
	tree->SetBranchAddress("loose_electron_charge",&loose_electron_charge, &b_loose_electron_charge); 
	

	//get tight muons
	tree->SetBranchAddress("muon_px",&muon_px, &b_muon_px); 
	tree->SetBranchAddress("muon_py",&muon_py, &b_muon_py);   
	tree->SetBranchAddress("muon_pz",&muon_pz, &b_muon_pz);
	tree->SetBranchAddress("muon_e",&muon_e, &b_muon_e); 
	tree->SetBranchAddress("muon_charge",&muon_charge,&b_muon_charge); 
	//get loose muons
	
	tree->SetBranchAddress("loose_muon_px",&loose_muon_px, &b_loose_muon_px); 
	tree->SetBranchAddress("loose_muon_py",&loose_muon_py, &b_loose_muon_py);   
	tree->SetBranchAddress("loose_muon_pz",&loose_muon_pz, &b_loose_muon_pz);
	tree->SetBranchAddress("loose_muon_e",&loose_muon_e, &b_loose_muon_e); 
	tree->SetBranchAddress("loose_muon_charge",&loose_muon_charge,&b_loose_muon_charge); 
	
	tree->SetBranchAddress("photon_px",&photon_px, &b_photon_px); 
	tree->SetBranchAddress("photon_py",&photon_py, &b_photon_py);   
	tree->SetBranchAddress("photon_pz",&photon_pz, &b_photon_pz);
	tree->SetBranchAddress("photon_e",&photon_e, &b_photon_e); 
	
	tree->SetBranchAddress("jet_px",&jet_px, &b_jet_px); 
	tree->SetBranchAddress("jet_py",&jet_py, &b_jet_py);   
	tree->SetBranchAddress("jet_pz",&jet_pz, &b_jet_pz);
	tree->SetBranchAddress("jet_e",&jet_e, &b_jet_e); 
	tree->SetBranchAddress("jet_unc",&jet_unc, &b_jet_unc); 
	tree->SetBranchAddress("jet_bTagL",&jet_bTagL, &b_jet_bTagL); 
	tree->SetBranchAddress("jet_bTagM",&jet_bTagM, &b_jet_bTagM); 
	tree->SetBranchAddress("jet_bTagT",&jet_bTagT, &b_jet_bTagT); 
	tree->SetBranchAddress("jet_algFlavor",&jet_algFlavor, &b_jet_algFlavor); 
	tree->SetBranchAddress("jet_phyFlavor",&jet_phyFlavor, &b_jet_phyFlavor); 

	
	tree->SetBranchAddress("w_px",&w_px, &b_w_px); 
	tree->SetBranchAddress("w_py",&w_py, &b_w_py);   
	tree->SetBranchAddress("w_pz",&w_pz, &b_w_pz);
	tree->SetBranchAddress("w_e",&w_e, &b_w_e); 
	
	tree->SetBranchAddress("z_px",&z_px, &b_z_px); 
	tree->SetBranchAddress("z_py",&z_py, &b_z_py);   
	tree->SetBranchAddress("z_pz",&z_pz, &b_z_pz);
	tree->SetBranchAddress("z_e",&z_e, &b_z_e); 
	
	tree->SetBranchAddress("top_px",&top_px, &b_top_px); 
	tree->SetBranchAddress("top_py",&top_py, &b_top_py);   
	tree->SetBranchAddress("top_pz",&top_pz, &b_top_pz);
	tree->SetBranchAddress("top_e",&top_e, &b_top_e); 
	
	
	
}

void bookEff(){

	double etaBin[]={0,0.8,1.6,2.4};
	
	for(int ieta=0; ieta<3; ieta++){
		for(int iflavor=5; iflavor>=3; iflavor--){
			TString title; 
			if(iflavor==5) title="b"; 
			if(iflavor==4) title="c";
			if(iflavor==3) title="uds";
			title = title+Form(", %.1f < |#eta| < %.1f", etaBin[ieta],etaBin[ieta+1]); 
			
			CreateEfficiency(Form("btagM_eff_eta%d_flavor%d",ieta,iflavor),title);
			CreateEfficiency(Form("btagT_eff_eta%d_flavor%d",ieta,iflavor),title);

			CreateEfficiency(Form("btagM_eff_eta%d_flavor%d_stmin",ieta,iflavor),title);
			CreateEfficiency(Form("btagT_eff_eta%d_flavor%d_stmin",ieta,iflavor),title);
			
		}
        int iflavor=21;
        TString title="gluon";
		title = title+Form(", %.1f < |#eta| < %.1f", etaBin[ieta],etaBin[ieta+1]);
        CreateEfficiency(Form("btagM_eff_eta%d_flavor%d",ieta,iflavor),title);
        CreateEfficiency(Form("btagT_eff_eta%d_flavor%d",ieta,iflavor),title);
        
        CreateEfficiency(Form("btagM_eff_eta%d_flavor%d_stmin",ieta,iflavor),title);
        CreateEfficiency(Form("btagT_eff_eta%d_flavor%d_stmin",ieta,iflavor),title);

    
    }
	
}

string jet_label(int flv){
    flv=TMath::Abs(flv);
    if(flv==1) return "d";
    if(flv==2) return "u";
    if(flv==3) return "s";
    if(flv==4) return "c";
    if(flv==5) return "b";
    if(flv==21) return "g";
    
    return "un-assigned";

    
}

void fill_event_category(){
  
    int ist=0;
    if(st>300) ist=1;
    if(st>500) ist=2;
    if(st>1000) ist=3;

    if(ist==0) cout << "ist is 0" << endl;
    


    
    TString scales[]={"_noSF","_SF","_SFP","_SFM","_SFPbc","_SFPl","_SFMbc","_SFMl","_TopCor","_TopCorP","_TopCorM"};
    int Nscale=sizeof(scales)/sizeof(TString); ;
    
    std::map<TString,double> event_weight;
    event_weight["_noSF"]=weight;
    
    event_weight["_SF"]=weight*B_weight;
    event_weight["_SFP"]=weight*B_weightP;
    event_weight["_SFM"]=weight*B_weightM;
    event_weight["_SFPbc"]=weight*B_weightPbc;
    event_weight["_SFMbc"]=weight*B_weightMbc;
    
    event_weight["_SFPl"]=weight*B_weightPl;
    event_weight["_SFMl"]=weight*B_weightMl;
    
    event_weight["_TopCor"]=weight*TTbar_corr*B_weight_tight;
    event_weight["_TopCorP"]=weight*TTbar_corrP*B_weight_tight;
    event_weight["_TopCorM"]=weight*TTbar_corrM*B_weight_tight;
    
	event_weight["_tight_noSF"]=weight; 
	
    event_weight["_tight_SF"]=weight*B_weight_tight;
    event_weight["_tight_SFP"]=weight*B_weight_tightP;
    event_weight["_tight_SFM"]=weight*B_weight_tightM;
    event_weight["_tight_SFPbc"]=weight*B_weight_tightPbc;
    event_weight["_tight_SFMbc"]=weight*B_weight_tightMbc;
    
    event_weight["_tight_SFPl"]=weight*B_weight_tightPl;
    event_weight["_tight_SFMl"]=weight*B_weight_tightMl;
    
    event_weight["_tight_TopCor"]=weight*TTbar_corr*B_weight_tight;
    event_weight["_tight_TopCorP"]=weight*TTbar_corrP*B_weight_tight;
    event_weight["_tight_TopCorM"]=weight*TTbar_corrM*B_weight_tight;
    
    TString tag[]={"","_tight"};
    int Ntag=sizeof(tag)/sizeof(TString);

    for(int itag=0; itag<Ntag; itag++){
        int NB=nBJets;
        if(tag[itag]=="_tight") NB=nTightBJets;
        if(NB>2)NB=2;
        
        for(int i=0; i<Nscale; i++){
            for(int j=1; j<=ist; j++){
            int iBin=0;
            if(j==1) iBin=nJets;
            if(j==2) iBin=nJets+nJetmax;
            if(j==3) iBin=nJets+2*nJetmax;
           
            if(tt_enriched) {
                sumW_EMu+=TTbar_corr;
                sum_EMu+=1;
                hName[Form("EventCategories_1Mu_1El_%dbtag",NB) + tag[itag] + scales[i]]->Fill(iBin,event_weight[tag[itag]+scales[i]]);
            }
            if(Z_enriched_loose) hName[Form("EventCategories_2Mu_onZ_%dbtag",NB)+tag[itag]+scales[i]]->Fill(iBin,event_weight[tag[itag]+scales[i]]);
            if(sig_2mu) hName[Form("EventCategories_2Mu_offZ_%dbtag",NB)+tag[itag]+scales[i]]->Fill(iBin,event_weight[tag[itag]+scales[i]]);
        }//ist loop for insclusive bin

    }
  
    }
}

TString create_category(int nJ, bool inclusive, int Nmu, int Nel, int charge, int Nb, TString tag, TString sf, TString mass){
    TString name="";
  
    if(nJ>0 && inclusive==0)name=name+Form("nJets%d_",nJ);
    if(nJ>0 && inclusive) name=name+Form("inclusive_nJets%d",nJ);
    
    for(int imu=0; imu<Nmu;imu++){
        name=name+"Mu";
    }

    for(int iel=0; iel<Nel;iel++){
        name=name+"El";
    }
    
    if(Nmu>=2 || Nmu+Nel>=2){
        if(charge==1)name=name+"_SS";
        else if(charge==-1) name=name+"_OS";
    }
    name=name+mass;
    
    if(Nb>0)name=name+Form("_%db_",Nb);
    if(Nb>0)name=name+tag+"-tag";
    if(Nb>0)name=name+sf;
    
    return name;
    
}

void bookHisto(){
    int NMU[]={1,2};
    int NEL[]={-1,1};
    int NB[]={-1,0,1,2};
    int Q[]={-1,1};
    int NJET[]={-1,1,2,3,4,5,6,7,8};
    TString mass[]={"","_on-z","_off-z"};
    TString TAG[]={"M","T"};
    TString SF[]={"","_SF","_SFP","_SFM"};
    bool INC[]={0,1};

    int totHist=0;
    
    for(int imu=0; imu<2;imu++)
        for(int iq=0; iq<2;iq++)
            for(int iel=0; iel<2; iel++)
                for(int inb=0; inb<4;inb++)
                    for(int inj=0;inj<9;inj++)
                        for(int im=0; im<3;im++)
                            for(int inc=0; inc<2; inc++)
                                for(int itag=0; itag<2; itag++)
                                    for(int isf=0; isf<3; isf++){
                                        if(imu>=1 && iel>=1) continue;
                                        if(INC[inc]==1 && NJET[inj]>6)continue;
                                       // cout << create_category(NJET[inj],INC[inc],NMU[imu],NEL[iel],Q[iq],NB[inb],TAG[itag],SF[isf],mass[im]) << endl;
                                        totHist++;
                            }
    cout << "totHist: " << totHist << endl;
    
	cout << "Book Histograms: " << endl;
	//float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1600,1800,2000,2200,2500,3000};
	float st_bins[]={200,300,450,600,750,1000,3000}; 
	float st_bins2[]={200,300,700,3000};

	double Dst_bins[]={200,300,450,600,750,1000,3000};

	int N=sizeof(st_bins)/sizeof(float)-1;
	int N2=sizeof(st_bins2)/sizeof(float)-1;
	
	float ht_bins[]={0,200,300,450,600,750,1000,3000};
	int Nht=sizeof(ht_bins)/sizeof(float)-1;
	
	float Xbins[]={0,30,50,100,150,250,350,600,800,1500};
	int NX=sizeof(Xbins)/sizeof(float)-1;
	
	float hMumu_bins[]={0,50,60,70,80,85,90,95,100,105,110,115,120,130,140,160,180,210,260,310,500,800};
	int NHmumu=sizeof(hMumu_bins)/sizeof(float)-1;
	
    string scales[]={"_noSF","_SF","_SFP","_SFM","_SFPbc","_SFPl","_SFMbc","_SFMl","_TopCor","_TopCorP","_TopCorM","_tight_noSF","_tight_SF","_tight_SFP","_tight_SFM","_tight_SFPbc","_tight_SFPl","_tight_SFMbc","_tight_SFMl","_tight_TopCor","_tight_TopCorP","_tight_TopCorM"};
    int Nscale=sizeof(scales)/sizeof(string);
    for(int nb=0; nb<=2; nb++){
        TString title=Form("#mu,e, %d-btags",nb);
        if(nb==2)title=Form("#mu,e, #geq %d-btags",nb);


        for(int is=0; is<Nscale; is++){
			TString sc=scales[is];
			if(sc.Contains("tight")) title=title+" tight";
			CreateHistogram(Form("EventCategories_1Mu_1El_%dbtag%s",nb,scales[is].c_str()),title,"","Events",(nJetmax)*3,0.5,nJetmax*3+0.5);
			CreateHistogram(Form("EventCategories_2Mu_offZ_%dbtag%s",nb,scales[is].c_str()),title,"","Events",(nJetmax)*3,0.5,nJetmax*3+0.5);
			CreateHistogram(Form("EventCategories_2Mu_onZ_%dbtag%s",nb,scales[is].c_str()),title,"","Events",(nJetmax)*3,0.5,nJetmax*3+0.5);
		}
        
        int nJ=1;
        for(int i=1; i<=hName[Form("EventCategories_1Mu_1El_%dbtag_noSF",nb)]->GetNbinsX();i++){
            TString binLabel=Form("%d",nJ);
            if(nJ==nJetmax)binLabel=Form("#geq %d",nJ);
           for(int is=0; is<Nscale; is++) hName[Form("EventCategories_1Mu_1El_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);
           for(int is=0; is<Nscale; is++) hName[Form("EventCategories_2Mu_offZ_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);
           for(int is=0; is<Nscale; is++) hName[Form("EventCategories_2Mu_onZ_%dbtag%s",nb,scales[is].c_str())]->GetXaxis()->SetBinLabel(i,binLabel);

            nJ++;
            if(nJ>nJetmax)nJ=1;
        }

    }
    

	CreateHistogram("gen_top_pt","","p_{T}(t) [GeV]", "Events", 50,0,1000);
	CreateHistogram("gen_top_pt_reweighted","","p_{T}(t) [GeV]", "Events", 50,0,1000);
    
    for(int ib=0; ib<=2; ib++)	CreateHistogram(Form("gen_top_pt_%db",ib),"","p_{T}(t) [GeV]", "Events", 50,0,1000);
    for(int ib=0; ib<=2; ib++)	CreateHistogram(Form("gen_top_pt_%db_reweighted",ib),"","p_{T}(t) [GeV]", "Events", 50,0,1000);

    
	CreateHistogram("top_ptweight_avg","","p_{T}(t) [GeV]", "Events", 100,0,1.5); 

    CreateHistogram("theta_E","","E [GeV]","#theta",100,0,1000,20,0,3.14);
    CreateHistogram("theta_M","","M_{bb} [GeV]","#theta",100,0,1000,20,0,3.14);

    CreateHistogram("deltaPhi_M","#Delta #Phi (bb RF) vs M(bb)","M_{bb} [GeV]","#Delta #Phi",100,0,1000,20,0,3.14);
    CreateHistogram("M_bb","#Delta #Phi (bb RF) vs M(bb)","M_{bb} [GeV]","Events",100,0,1000);
    
    
    CreateHistogram("Acceptance_num_0b_st_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);
    CreateHistogram("Acceptance_num_1b_st_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);
    CreateHistogram("Acceptance_num_2b_st_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);

    CreateHistogram("Acceptance_den_st_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);
    
    CreateHistogram("Acceptance_num_0b_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",10,0.5,10.5);
    CreateHistogram("Acceptance_num_1b_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",10,0.5,10.5);
    CreateHistogram("Acceptance_num_2b_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",10,0.5,10.5);

    CreateHistogram("Acceptance_den_nJets_1Mu_1El","#mu=1","S_{T} [GeV]","n-jets",10,0.5,10.5);

    
    
	CreateHistogram("st_nJets_1Mu","#mu=1","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);
	CreateHistogram("st_nJets_2Mu","#mu=2","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);

	CreateHistogram("st_nJets_W","W","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);
	CreateHistogram("st_nJets_Z","Z","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);
	CreateHistogram("st_nJets_tt","tt","S_{T} [GeV]","n-jets",N,Dst_bins,8,0.5,8.5);

    for(int ib=0; ib<=2; ib++){
        CreateHistogram(Form("st_1Mu_1El_%dbtag_tight_SF",ib),"#mu,e","S_{T} [GeV]","Events",N,st_bins);
        CreateHistogram(Form("st_1Mu_1El_%dbtag_tight_TopCor",ib),"#mu,e","S_{T} [GeV]","Events",N,st_bins);
        
        CreateHistogram(Form("st_1Mu_1El_%dbtag_SF",ib),"#mu,e","S_{T} [GeV]","Events",N,st_bins);
        CreateHistogram(Form("st_1Mu_1El_%dbtag_TopCor",ib),"#mu,e","S_{T} [GeV]","Events",N,st_bins);
        
        CreateHistogram(Form("met_1Mu_1El_%dbtag_tight_SF",ib),"#mu,e","S_{T} [GeV]","Events",NX,Xbins);
        CreateHistogram(Form("met_1Mu_1El_%dbtag_tight_TopCor",ib),"#mu,e","S_{T} [GeV]","Events",NX,Xbins);
        
        CreateHistogram(Form("bjet1pt_1Mu_1El_%dbtag_tight_SF",ib),"#mu,e","S_{T} [GeV]","Events",NX,Xbins);
        CreateHistogram(Form("bjet1pt_1Mu_1El_%dbtag_tight_TopCor",ib),"#mu,e","S_{T} [GeV]","Events",NX,Xbins);
        
        CreateHistogram(Form("bjet2pt_1Mu_1El_%dbtag_tight_SF",ib),"#mu,e","S_{T} [GeV]","Events",NX,Xbins);
        CreateHistogram(Form("bjet2pt_1Mu_1El_%dbtag_tight_TopCor",ib),"#mu,e","S_{T} [GeV]","Events",NX,Xbins);
        
        CreateHistogram(Form("ptmu_1Mu_1El_%dbtag_tight_SF",ib),"#mu,e","p_{T}(#mu) [GeV]","Events",NX,Xbins);
        CreateHistogram(Form("ptmu_1Mu_1El_%dbtag_tight_TopCor",ib),"#mu,e","p_{T}(#mu) [GeV]","Events",NX,Xbins);

        CreateHistogram(Form("pte_1Mu_1El_%dbtag_tight_SF",ib),"#mu,e","p{T}(e) [GeV]","Events",NX,Xbins);
        CreateHistogram(Form("pte_1Mu_1El_%dbtag_tight_TopCor",ib),"#mu,e","p{T}(e) [GeV]","Events",NX,Xbins);

    }
	
	
	CreateHistogram("h_wpt","w p_{T}", "p_{T} [GeV]", "Events", 50,0,400); 
	CreateHistogram("h_wrapidity","w y", "y", "Events", 50,-5,5); 
		
	CreateProfile("pr_deltapt_deltaphi","#delta p_{T} vs #delta #phi(W,jet system)","#delta #phi(W,jet system) [GeV]", "#delta p_{T}[GeV]", 50,-3.14,3.14);
	CreateProfile("pr_deltapt_wpt","#delta p_{T} vs p_{T}(W)","p_{T}(W) [GeV]", "#delta p_{T}[GeV]", 50,0,500);
	CreateProfile("pr_deltapt_ttpt","#delta p_{T} vs p_{T}(tt)","p_{T}(tt) [GeV]", "#delta p_{T}[GeV]", 100,0,1000);
	CreateProfile("pr_deltapt_zpt","#delta p_{T} vs p_{T}(Z)","p_{T}(Z) [GeV]", "#delta p_{T}[GeV]", 50,0,500);

	CreateProfile("pr_deltaphi_wpt","#delta #phi(W,jet system) vs p_{T}(W)","p_{T}(W) [GeV]", "#delta #phi(W,jet system)", 50,0,500);

	CreateHistogram("h_deltapt","#delta p_{T} ", "#delta p_{T}[GeV]", "Events", 50,-50,50);
	CreateHistogram("h_deltaphi","#delta #phi ", "#delta #phi", "Events", 50,-3.14,3.14);

    CreateHistogram("T20_nJets_num","","n-jets","T20",10,0.5,10.5);
    CreateHistogram("T20_nJets_den","","n-jets","T20",10,0.5,10.5);

    
    for(int nJ=2; nJ<=nJetmax; nJ++){
        for(int ntag=0; ntag<=3; ntag++){
            CreateHistogram(Form("deltaPt_1Mu_1El_nJets%d_%dbtag",nJ,ntag),"","#Delta p_{T}","Events",100,0,1000);

        }
        for(int ib=0; ib<=2; ib++){
            if(ib==2)CreateHistogram(Form("h_T20_nJets%d_nbtruth%d",nJ,ib),"","T_{20}","Events",200,0,2);
            else CreateHistogram(Form("h_T20_nJets%d_nbtruth%d",nJ,ib),"","T_{20}","Events",1000,0,1000);
        }
    }
    
	CreateHistogram("h_wpt_jpt","jet p_{T} vs w p_{T}", "p_{T} [GeV](W)", "p_{T}(jet systemt)", 50,0,400,50,0,400); 

    CreateHistogram("jet_pt_1Mu_1El_btag_tight_SF","","jet p_{T} [GeV]","Events", 50,30,1000);
    
    CreateHistogram("btag_weight_tight","jet configruation weight","n-jet","n_{b}",10,0.5,10.5,3,0.5,3.5);
    CreateHistogram("btag_weight","jet configruation weight","n-jet","n_{b}",10,0.5,10.5,3,0.5,3.5);

    CreateProfile("probMC_nJet_nB_tight","jet configruation probability (MC) tight","n-jet","n_{b}",10,0.5,10.5,5,-0.5,4.5,0,1);

    
	CreateHistogram("h_Loose_Medium","n_{b} loose vs medium","n_{b} medium","n_{b} loose",10,0.5,10.5,10,0.5,10.5); 
	CreateHistogram("h_Medium_Tight","n_{b} medium vs tight","n_{b} tight","n_{b} medium",10,0.5,10.5,10,0.5,10.5); 
	CreateHistogram("h_Loose_Tight","n_{b} loose vs tight","n_{b} tight","n_{b} loose",10,0.5,10.5,10,0.5,10.5); 

	
	CreateHistogram("rapidity","y","y","Events",50,-4,4); 
	CreateHistogram("rapidity_st","y vs S_{T}","S_{T} [GeV]","y",100,0,3000, 50,-4,4); 

	CreateHistogram("nJets_1Mu","#mu=1","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu","#mu=2","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_1Mu_0btag","#mu=1, 0 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1btag","#mu=1, 1 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_2btag","#mu=1, #geq 2 b-tags","n-jets","Events",10,0.5,10.5); 

	CreateHistogram("nJets_2Mu_0btag","#mu=2, 0 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu_1btag","#mu=2,  1 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu_2btag","#mu=2,  #geq 2 b-tags","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_1Mu_0btag_tight","#mu=1, 0 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1btag_tight","#mu=1, 1 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_2btag_tight","#mu=1, #geq 2 b-tags","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_2Mu_0btag_tight","#mu=2, 0 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu_1btag_tight","#mu=2,  1 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu_2btag_tight","#mu=2,  #geq 2 b-tags","n-jets","Events",10,0.5,10.5); 
	
    CreateHistogram("nJets_ss2Mu_0btag_tight","#mu=2, 0 b-tags","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ss2Mu_1btag_tight","#mu=2,  1 b-tags","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ss2Mu_2btag_tight","#mu=2,  #geq 2 b-tags","n-jets","Events",10,0.5,10.5);
    
	
	CreateHistogram("nJets_1Mu_0btag_SF","#mu=1, 0 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1btag_SF","#mu=1, 1 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_2btag_SF","#mu=1, #geq 2 b-tags","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_2Mu_0btag_SF","#mu=2, 0 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu_1btag_SF","#mu=2,  1 b-tags","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_2Mu_2btag_SF","#mu=2,  #geq 2 b-tags","n-jets","Events",10,0.5,10.5); 
	
	
	CreateHistogram("nBJetsLoose_1Mu","b-tags csv-loose, #mu=1","b-tags","Events",11,-0.5,10.5); 
	CreateHistogram("nBJetsLoose_2Mu","b-tags csv-loose, #mu=2","b-tags","Events",11,-0.5,10.5); 
	
	CreateHistogram("nBJetsMedium_1Mu","b-tags csv-medium, #mu=1","b-tags","Events",11,-0.5,10.5); 
	CreateHistogram("nBJetsMedium_2Mu","b-tags csv-medium, #mu=2","b-tags","Events",11,-0.5,10.5); 
	
	CreateHistogram("nBJetsTight_1Mu","b-tags csv-Tight, #mu=1","b-tags","Events",11,-0.5,10.5);
	CreateHistogram("nBJetsTight_2Mu","b-tags csv-Tight, #mu=2","b-tags","Events",11,-0.5,10.5);
    CreateHistogram("nBJetsTight_1Mu_1El","b-tags csv-Tight, 1El, #mu=1","b-tags","Events",11,-0.5,10.5);

    
	CreateHistogram("nJets_W","#mu=1","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_W_0btag","#mu=1, 0 b-tags ","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_W_1btag","#mu=1, 1 b-tags ","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_W_2btag","#mu=1, #geq 2 b-tags ","n-jets","Events",10,0.5,10.5); 

	CreateHistogram("nJets_W_0btag_SF","#mu=1, 0 b-tags ","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_W_1btag_SF","#mu=1, 1 b-tags ","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_W_2btag_SF","#mu=1, #geq 2 b-tags ","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_Z","nJets, #mu=2, 87<M_{#mu#mu}<95 GeV","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_0btag","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_SF","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_SF","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_SF","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_SFP","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_SFP","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_SFP","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_SFM","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_SFM","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_SFM","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_tight","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_tight","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_tight","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_tight_SF","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_tight_SF","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_tight_SF","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_tight_SFP","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_tight_SFP","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_tight_SFP","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_Z_0btag_tight_SFM","#mu=2, 87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_1btag_tight_SFM","#mu=2, 87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_Z_2btag_tight_SFM","#mu=2, 87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
    //ss
    CreateHistogram("nJets_ssZ_0btag_tight_SF","#mu=2, ss,87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ssZ_1btag_tight_SF","#mu=2, ss,87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ssZ_2btag_tight_SF","#mu=2, ss,87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_ssZ_0btag_tight_SFP","#mu=2, ss,87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ssZ_1btag_tight_SFP","#mu=2, ss,87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ssZ_2btag_tight_SFP","#mu=2, ss,87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
	CreateHistogram("nJets_ssZ_0btag_tight_SFM","#mu=2, ss,87<M_{#mu#mu}<95 GeV, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ssZ_1btag_tight_SFM","#mu=2, ss,87<M_{#mu#mu}<95 GeV, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_ssZ_2btag_tight_SFM","#mu=2, ss,87<M_{#mu#mu}<95 GeV, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
    CreateProfile("pr_probMC_nJets_1Mu_1El_0btag_tight","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
    CreateProfile("pr_probMC_nJets_1Mu_1El_1btag_tight","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
    CreateProfile("pr_probMC_nJets_1Mu_1El_2btag_tight","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
    
    CreateHistogram("h_probMC_nJets_1Mu_1El_0btag_tight","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
    CreateHistogram("h_probMC_nJets_1Mu_1El_1btag_tight","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
    CreateHistogram("h_probMC_nJets_1Mu_1El_2btag_tight","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
	CreateHistogram("nJets_1Mu_1El","#mu=1, 1 El","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_0btag","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
    
    CreateProfile("pr_weight_nJets_1Mu_1El_0btag_tight","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
    CreateProfile("pr_weight_nJets_1Mu_1El_1btag_tight","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
    CreateProfile("pr_weight_nJets_1Mu_1El_2btag_tight","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
	CreateHistogram("nJets_1Mu_1El_0btag_tight","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_tight","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_tight","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
    for(int ib=0; ib<=2; ib++){
        for(int nbt=0; nbt<=2;nbt++){
            CreateHistogram(Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",nbt,ib),"#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
        }
    }
	
	CreateHistogram("nJets_ss1Mu_1El_0btag_tight_SF","ss #mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_ss1Mu_1El_1btag_tight_SF","ss #mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_ss1Mu_1El_2btag_tight_SF","ss #mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_ss1Mu_1El_0btag_tight_SFP","ss #mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_ss1Mu_1El_1btag_tight_SFP","ss #mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_ss1Mu_1El_2btag_tight_SFP","ss #mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_ss1Mu_1El_0btag_tight_SFM","ss #mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_ss1Mu_1El_1btag_tight_SFM","ss #mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_ss1Mu_1El_2btag_tight_SFM","ss #mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
    
	CreateHistogram("nJets_1Mu_1El_0btag_tight_SF","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SF","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SF","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_1Mu_1El_0btag_tight_SFP","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SFP","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SFP","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_1Mu_1El_0btag_tight_SFM","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SFM","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SFM","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
    CreateHistogram("nJets_1Mu_1El_0btag_tight_SFP_bc","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SFP_bc","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SFP_bc","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
	
    CreateHistogram("nJets_1Mu_1El_0btag_tight_SFM_bc","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SFM_bc","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SFM_bc","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
    CreateHistogram("nJets_1Mu_1El_0btag_tight_SFP_l","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SFP_l","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SFP_l","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
	CreateHistogram("nJets_1Mu_1El_0btag_tight_SFM_l","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_1btag_tight_SFM_l","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5);
	CreateHistogram("nJets_1Mu_1El_2btag_tight_SFM_l","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5);
    
	CreateHistogram("nJets_1Mu_1El_0btag_SF","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_SF","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_SF","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_1Mu_1El_0btag_SFP","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_SFP","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_SFP","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("nJets_1Mu_1El_0btag_SFM","#mu=1, 1 El, 0 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_1btag_SFM","#mu=1, 1 El, 1 b-tag","n-jets","Events",10,0.5,10.5); 
	CreateHistogram("nJets_1Mu_1El_2btag_SFM","#mu=1, 1 El, #geq 2 b-tag","n-jets","Events",10,0.5,10.5); 
	
	CreateHistogram("jetpt_nJets","p_{T}(jet system) vs nJets","nJets","p_{T} (jet system)[GeV]",10,0.5,10.5, 20,0,1000); 

	
	CreateHistogram("h_before_npv","pile-up","N_{pv}","Events",50,0,50); 
	CreateHistogram("h_npv","pile-up","N_{pv}","Events",50,0,50); 
	CreateHistogram("h_npv_nJets4","pile-up","N_{pv}","Events",50,0,50); 

	CreateProfile("R_dataMC_eff","Ratio #epsilon(data)/#epsilon(MC)","p_{T}(#mu) [GeV]","R(data/MC)",50,0,250); 
	CreateProfile("pr_nJets_npv", "pile-up", "N_{pv}", "n-jets",50,0,50); 
	
    CreateHistogram("h_mumu","M_{#mu#mu} inclusive","Mass [GeV]","Events",200,0,1000);
    CreateHistogram("h_ssmumu","M_{#mu#mu} same sign inclusive","Mass [GeV]","Events",200,0,1000);

    CreateHistogram("jet_flavor","jet flavor", "jet flavor", "Events",7,0,7);
    for(int i=0; i<=5;i++)CreateHistogram(Form("jet_pt_flavor%d",i),Form("jet pt, flavor %d",i), "jet p_{T} [GeV]", "Events",50,0,1000);
    CreateHistogram("jet_pt_flavor21","jet pt, flavor 21", "jet p_{T} [GeV]", "Events",50,0,1000);
    
	for(int nJ=0; nJ<=nJetmax; nJ++){
		
        CreateProfile(Form("jet_flavor_nJets%d",nJ),Form("jet flavor,%d-jets",nJ), "jet flavor", "N-jets/Event",4,0,4);
        for(int ib=0; ib<=2; ib++)CreateProfile(Form("jet_flavor_nJets%d_%dbtag",nJ,ib),Form("jet flavor,%d-jets, %d-btags",nJ,ib), "jet flavor", "N-jets/Event",4,0,4);
        for(int ib=0; ib<=2; ib++){
            CreateHistogram(Form("nb_nJets%d_%dbtag",nJ,ib),Form("nb,%d-jets, %d-btags",nJ,ib), "n_{b}^{MC Truth}", "Events",nJ,0.5,nJ);
            CreateHistogram(Form("nc_nJets%d_%dbtag",nJ,ib),Form("nc,%d-jets, %d-btags",nJ,ib), "n_{c}^{MC Truth}", "Events",nJ,0.5,nJ);
            CreateHistogram(Form("ng_nJets%d_%dbtag",nJ,ib),Form("ng,%d-jets, %d-btags",nJ,ib), "n_{g}^{MC Truth}", "Events",nJ,0.5,nJ);
            CreateHistogram(Form("nl_nJets%d_%dbtag",nJ,ib),Form("nl,%d-jets, %d-btags",nJ,ib), "n_{l}^{MC Truth}", "Events",nJ,0.5,nJ);

        }

        
		CreateHistogram(Form("st_ptsys_nJets%d",nJ),Form("%d-jets",nJ),"S_{T} [GeV]","p_{T}(jet system) [GeV]",100,0,3000, 50, 0, 1000); 
		CreateHistogram(Form("st_met_nJets%d",nJ),Form("%d-jets",nJ),"S_{T} [GeV]","M_{ET} [GeV]",100,0,3000, 50, 0, 1000); 

		CreateProfile(Form("pr_corr_nJets%d_st",nJ), Form("Correction Factor, %d-jets",nJ), "S_{T} [GeV]","Corr", N,st_bins); 
		CreateProfile(Form("pr_corr_nJets%d_st_1Mu_1El_btag",nJ), Form("Correction Factor, %d-jets",nJ), "S_{T} [GeV]","Corr", N,st_bins); 

		for(int ijet=0; ijet<nJ; ijet++){
            CreateHistogram(Form("W_jetpt_jet%d_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, W",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000);
			CreateHistogram(Form("Z_jetpt_jet%d_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, Z",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000);
			CreateHistogram(Form("tt_jetpt_jet%d_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, tt",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000);
            for(int ib=0; ib<=2; ib++) CreateHistogram(Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight",ijet,nJ,ib),Form("jet p_{T}, %dth jet, %d-jets, %db-tags, tt",ijet,nJ,ib),"p_{T}(jet) [GeV]","Events", 20,0,1000);
              for(int ib=0; ib<=2; ib++) CreateHistogram(Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight_SFP",ijet,nJ,ib),Form("jet p_{T}, %dth jet, %d-jets, %db-tags, tt",ijet,nJ,ib),"p_{T}(jet) [GeV]","Events", 20,0,1000);
            for(int ib=0; ib<=2; ib++) CreateHistogram(Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight_SFM",ijet,nJ,ib),Form("jet p_{T}, %dth jet, %d-jets, %db-tags, tt",ijet,nJ,ib),"p_{T}(jet) [GeV]","Events", 20,0,1000);
              for(int ib=0; ib<=2; ib++) CreateHistogram(Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight_SF",ijet,nJ,ib),Form("jet p_{T}, %dth jet, %d-jets, %db-tags, tt",ijet,nJ,ib),"p_{T}(jet) [GeV]","Events", 20,0,1000);

			CreateHistogram(Form("W_jetpt_jet%d_st1_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, W",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000);
			CreateHistogram(Form("Z_jetpt_jet%d_st1_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, Z",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			CreateHistogram(Form("tt_jetpt_jet%d_st1_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, tt",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 

			CreateHistogram(Form("W_jetpt_jet%d_st2_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, W",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			CreateHistogram(Form("Z_jetpt_jet%d_st2_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, Z",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			CreateHistogram(Form("tt_jetpt_jet%d_st2_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, tt",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			
			CreateHistogram(Form("W_jetpt_jet%d_st3_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, W",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			CreateHistogram(Form("Z_jetpt_jet%d_st3_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, Z",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			CreateHistogram(Form("tt_jetpt_jet%d_st3_nJets%d",ijet,nJ),Form("jet p_{T}, %dth jet, %d-jets, tt",ijet,nJ),"p_{T}(jet) [GeV]","Events", 20,0,1000); 
			
		}
		
		
		CreateHistogram(Form("Z_nJets%d_dileptonpT",nJ),Form("%d-jets, Z",nJ), "p_{T}(#mu#mu) [GeV]","Events",75,0,1500); 
		CreateHistogram(Form("tt_nJets%d_dileptonpT",nJ),Form("%d-jets, tt",nJ), "p_{T}(ll) [GeV]","Events",75,0,1500); 
		CreateHistogram(Form("W_nJets%d_muonpT",nJ),Form("%d-jets, W",nJ), "p_{T}(#mu) [GeV]","Events",75,0,1500); 

		for(int mu=1; mu<=2; mu++){
			CreateHistogram(Form("Z_nJets%d_muonpT%d",nJ,mu),Form("%d-jets,#mu=%d, Z",nJ,mu), "p_{T}(#mu) [GeV]","Events",75,0,1500); 
			CreateHistogram(Form("tt_nJets%d_muonpT%d",nJ,mu),Form("%d-jets,#mu=%d, tt",nJ,mu), "p_{T}(#mu) [GeV]","Events",75,0,1500); 
		}
		
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
		
		CreateHistogram(Form("h_mumu_nJet%d",nJ),"M_{#mu#mu}","Mass [GeV]","Events",200,0,1000);
		CreateHistogram(Form("h_mumu_nJet%d_0btag",nJ),"M_{#mu#mu}, 0 b-tag","Mass [GeV]","Events",200,0,1000);
		CreateHistogram(Form("h_mumu_nJet%d_1btag",nJ),"M_{#mu#mu}, 1 b-tag","Mass [GeV]","Events",200,0,1000);
		CreateHistogram(Form("h_mumu_nJet%d_2btag",nJ),"M_{#mu#mu}, #geq 2 b-tag","Mass [GeV]","Events",200,0,1000);
		if(nJ==5){
			nJ=56;
			
			CreateHistogram(Form("h_mumu_nJet%d_0btag_st1",nJ),"M_{#mu#mu}, 0 b-tag, 300<S_{T}<700","Mass [GeV]","Events",200,0,1000);
			CreateHistogram(Form("h_mumu_nJet%d_1btag_st1",nJ),"M_{#mu#mu}, 1 b-tag, 300<S_{T}<700","Mass [GeV]","Events",200,0,1000);
			CreateHistogram(Form("h_mumu_nJet%d_2btag_st1",nJ),"M_{#mu#mu}, #geq 2 b-tag, 300<S_{T}<700","Mass [GeV]","Events",200,0,1000);
			
			CreateHistogram(Form("h_mumu_nJet%d_0btag_st2",nJ),"M_{#mu#mu}, 0 b-tag, 700<S_{T}<3000","Mass [GeV]","Events",200,0,1000);
			CreateHistogram(Form("h_mumu_nJet%d_1btag_st2",nJ),"M_{#mu#mu}, 1 b-tag, 700<S_{T}<3000","Mass [GeV]","Events",200,0,1000);
			CreateHistogram(Form("h_mumu_nJet%d_2btag_st2",nJ),"M_{#mu#mu}, #geq 2 b-tag,700<S_{T}<3000","Mass [GeV]","Events",200,0,1000);
			nJ=5;
		}
		CreateHistogram(Form("h_muel_nJet%d",nJ),"M_{#mu e}","Mass [GeV]","Events",200,0,1000);

		CreateHistogram(Form("h_mumu_nJet%d_0p5GeV",nJ),"M_{#mu#mu}","Mass [GeV]","Events",2000,0,1000);

		
		CreateHistogram(Form("ht_nJets%d",nJ),Form("%d jets",nJ),"H_{T} [GeV]", "Events",N,st_bins);
		for(int pT=30; pT<=100; pT+=10){
			CreateHistogram(Form("st_nJets%d_pTmin%d",nJ, pT),Form("%d jets, p_{T}^{min}=%d",nJ,pT), "S_{T} [GeV]", "Events",N,st_bins); 
			CreateHistogram(Form("st_nJets%d_pTmin%d_uW",nJ, pT),Form("%d jets, p_{T}^{min}=%d",nJ,pT), "S_{T} [GeV]", "Events",N,st_bins); 
			
		}//pT loop 
		
		CreateHistogram(Form("met_nJets%d_1Mu_1El",nJ),Form("%d jets, #mu,e",nJ),"M_{ET} [GeV]","Events",NX, Xbins);
		CreateHistogram(Form("met_nJets%d_W",nJ),Form("%d jets, %d#mu, 70<M_{T}(#mu,M_{ET})<250 GeV",nJ,1),"M_{ET} [GeV]","Events",NX, Xbins);
		CreateHistogram(Form("met_nJets%d_Z",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV",nJ,2),"M_{ET} [GeV]","Events",NX, Xbins);
		
        CreateHistogram(Form("met_nJets%d_ss1Mu_1El_0btag",nJ),Form("%d jets, #mu,e",nJ),"M_{ET} [GeV]","Events",NX, Xbins);
        CreateHistogram(Form("met_nJets%d_ss1Mu_1El_1btag",nJ),Form("%d jets, #mu,e",nJ),"M_{ET} [GeV]","Events",NX, Xbins);
        CreateHistogram(Form("met_nJets%d_ss1Mu_1El_2btag",nJ),Form("%d jets, #mu,e",nJ),"M_{ET} [GeV]","Events",NX, Xbins);

		CreateHistogram(Form("metsig_nJets%d_1Mu_1El",nJ),Form("%d jets, #mu,e",nJ),"M_{ET}/#sqrt{H_{T}} #sqrt{[GeV]}","Events",60,0,30);
		CreateHistogram(Form("metsig_nJets%d_W",nJ),Form("%d jets, %d#mu, 70<M_{T}(#mu,M_{ET})<250 GeV",nJ,1),"M_{ET}/#sqrt{H_{T}} #sqrt{[GeV]}","Events",60,0,30);
		CreateHistogram(Form("metsig_nJets%d_Z",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV",nJ,2),"M_{ET}/#sqrt{H_{T}} #sqrt{[GeV]}","Events",60,0,30);
		
		
		CreateHistogram(Form("ht_nJets%d_1Mu_1El",nJ),Form("%d jets, #mu,e",nJ),"H_{T} [GeV]","Events",N,ht_bins);
		CreateHistogram(Form("ht_nJets%d_W",nJ),Form("%d jets, %d#mu, 70<M_{T}(#mu,M_{ET})<250 GeV",nJ,1),"H_{T} [GeV]","Events",N,ht_bins);
		CreateHistogram(Form("ht_nJets%d_Z",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV",nJ,2),"H_{T} [GeV]","Events",N,ht_bins);
		
		CreateHistogram(Form("ht_inclusive_nJets%d_1Mu_1El",nJ),Form("%d jets, #mu,e",nJ),"H_{T} [GeV]","Events",N,ht_bins);
		CreateHistogram(Form("ht_inclusive_nJets%d_W",nJ),Form("%d jets, %d#mu, 70<M_{T}(#mu,M_{ET})<250 GeV",nJ,1),"H_{T} [GeV]","Events",N,ht_bins);
		CreateHistogram(Form("ht_inclusive_nJets%d_Z",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV",nJ,2),"H_{T} [GeV]","Events",N,ht_bins);
		
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_uW",nJ),Form("%d jets, #mu,e",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El",nJ),Form("%d jets, #mu,e",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);

		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SFP_topreweighted",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SFP_topreweighted",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SFP_topreweighted",nJ),Form("%d jets, #mu,e, #geq 2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SFM_topreweighted",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SFM_topreweighted",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SFM_topreweighted",nJ),Form("%d jets, #mu,e, #geq 2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SF_topreweighted",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SF_topreweighted",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SF_topreweighted",nJ),Form("%d jets, #mu,e, #geq 2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);

        CreateHistogram(Form("weight_nJets%d_1Mu_1El_0btag_tight",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"B_weight","Events",100,0,2);
        CreateHistogram(Form("weight_nJets%d_1Mu_1El_1btag_tight",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"B_weight","Events",100,0,2);
        CreateHistogram(Form("weight_nJets%d_1Mu_1El_2btag_tight",nJ),Form("%d jets, #mu,e, #geq 2 b-tag",nJ),"B_weight","Events",100,0,2);
        
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SF",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SF",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SF",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
        CreateHistogram(Form("stTop_nJets%d_1Mu_1El_0btag_tight_SF",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("stTop_nJets%d_1Mu_1El_1btag_tight_SF",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("stTop_nJets%d_1Mu_1El_2btag_tight_SF",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
        
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SFP",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SFP",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SFP",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SFM",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SFM",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SFM",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
	
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_0btag_tight_SF",nJ),Form("%d jets, ss #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_1btag_tight_SF",nJ),Form("%d jets, ss #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_2btag_tight_SF",nJ),Form("%d jets, ss #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_0btag_tight_SFP",nJ),Form("%d jets, ss #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_1btag_tight_SFP",nJ),Form("%d jets, ss #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_2btag_tight_SFP",nJ),Form("%d jets, ss #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_0btag_tight_SFM",nJ),Form("%d jets, ss #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_1btag_tight_SFM",nJ),Form("%d jets, ss #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_ss1Mu_1El_2btag_tight_SFM",nJ),Form("%d jets, ss #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		if(nJ==5){
			int nJ56=56;
			CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SF",nJ56),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SF",nJ56),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SF",nJ56),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
		
			CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SFP",nJ56),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SFP",nJ56),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SFP",nJ56),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
		
			CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight_SFM",nJ56),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight_SFM",nJ56),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight_SFM",nJ56),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			
			
			CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_SF",nJ56),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_SF",nJ56),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_SF",nJ56),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_SFP",nJ56),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_SFP",nJ56),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_SFP",nJ56),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_SFM",nJ56),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_SFM",nJ56),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_SFM",nJ56),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N2,st_bins2);
			
		}
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_tight",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_tight",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_tight",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_SF",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_SF",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_SF",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_SFP",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_SFP",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_SFP",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_1Mu_1El_0btag_SFM",nJ),Form("%d jets, #mu,e, 0 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_1btag_SFM",nJ),Form("%d jets, #mu,e, 1 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_1Mu_1El_2btag_SFM",nJ),Form("%d jets, #mu,e, #geq2 b-tag",nJ),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_W_0btag",nJ),Form("%d jets, %d#mu, 0 b-tags",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_W_1btag",nJ),Form("%d jets, %d#mu, 1 b-tags",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_W_2btag",nJ),Form("%d jets, %d#mu, #geq2 b-tags",nJ,1),"S_{T} [GeV]","Events",N,st_bins);

		CreateHistogram(Form("st_nJets%d_W_0btag_SF",nJ),Form("%d jets, %d#mu, 0 b-tags",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_W_1btag_SF",nJ),Form("%d jets, %d#mu, 1 b-tags",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_W_2btag_SF",nJ),Form("%d jets, %d#mu, #geq2 b-tags",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_W",nJ),Form("%d jets, %d#mu",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_W_uW",nJ),Form("%d jets, %d#mu, ",nJ,1),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);

		CreateHistogram(Form("st_nJets%d_Z_0btag_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag_tight",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_tight",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_tight",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag_tight_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_tight_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_tight_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag_tight_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_tight_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_tight_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("st_nJets%d_Z_0btag_tight_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_1btag_tight_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		CreateHistogram(Form("st_nJets%d_Z_2btag_tight_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		if(nJ==5){
			nJ=56;
			CreateHistogram(Form("st_nJets%d_Z_0btag_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_Z_0btag_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_Z_0btag_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_Z_0btag_tight",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_tight",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_tight",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_Z_0btag_tight_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_tight_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_tight_SF",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_Z_0btag_tight_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_tight_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_tight_SFP",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			
			CreateHistogram(Form("st_nJets%d_Z_0btag_tight_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 0 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_1btag_tight_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, 1 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			CreateHistogram(Form("st_nJets%d_Z_2btag_tight_SFM",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV, #geq2 b-tags",nJ,2),"S_{T} [GeV]","Events",N2,st_bins2);
			nJ=5;
		}
		
		
		CreateHistogram(Form("st_nJets%d_Z_uW",nJ),Form("%d jets, %d#mu, 87<M_{#mu#mu}<95 GeV",nJ,2),"S_{T} [GeV]","Events",N,st_bins);
		
		CreateHistogram(Form("jetpT%d_W",nJ),Form("Jet p_{T} spectrum for %dth jet,70<M_{T}(#mu,M_{ET})<250 GeV",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);			
		CreateHistogram(Form("jetpT%d_Z",nJ),Form("Jet p_{T} spectrum for %dth jet, 87<M_{#mu#mu}<95 GeV",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);
		CreateHistogram(Form("jetpT%d_ssZ",nJ),Form("Jet p_{T} spectrum for %dth jet, ss 87<M_{#mu#mu}<95 GeV",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);

        
		CreateHistogram(Form("bjetpT%d_Z",nJ),Form("Jet p_{T} spectrum for %dth jet, 87<M_{#mu#mu}<95 GeV, csv medium b-tag",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);			
		CreateHistogram(Form("bjetpT%d_Z_tight",nJ),Form("Jet p_{T} spectrum for %dth jet, 87<M_{#mu#mu}<95 GeV, csv tight b-tag",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);			

		CreateHistogram(Form("jetpT%d_1Mu_1El",nJ),Form("Jet p_{T} spectrum for %dth jet,#mu,e",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);		
		CreateHistogram(Form("bjetpT%d_1Mu_1El",nJ),Form("Jet p_{T} spectrum for %dth jet,#mu,e, csv medium b-tag",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);		
		CreateHistogram(Form("bjetpT%d_1Mu_1El_tight",nJ),Form("Jet p_{T} spectrum for %dth jet,#mu,e, csv tight b-tag",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);		
		CreateHistogram(Form("0bjetpT%d_1Mu_1El_tight",nJ),Form("Jet p_{T} spectrum for %dth jet,#mu,e, false csv tight b-tag",nJ),"Jet p_{T} [GeV]","Events",NX, Xbins);

		
		CreateHistogram(Form("jetpTsys%d_W",nJ),Form("%d-jets jet,70<M_{T}(#mu,M_{ET})<250 GeV",nJ),"Jet p_{T}(system) [GeV]","Events",NX, Xbins);		
		CreateHistogram(Form("jetpTsys%d_Z",nJ),Form("%d-jets, 87<M_{#mu#mu}<95 GeV",nJ),"Jet p_{T}(system) [GeV]","Events",NX, Xbins);			
		CreateHistogram(Form("jetpTsys%d_1Mu_1El",nJ),Form("%d-jets, #mu,e",nJ),"Jet p_{T}(system) [GeV]","Events",NX, Xbins);			
		
		
		
	}//first nJ loop for inclusive bins  
	
	for(int mu=1; mu<=2; mu++){		
		for(int nJ=0; nJ<=nJetmax; nJ++){
			
			CreateHistogram(Form("nJets%d_DeltaPhiMuonMet_%dMu",nJ,mu),"#Delta #phi","#Delta #phi(#mu,M_{ET})","Events",20,0,3.15); 
			//CreateHistogram(Form("nJets%d_muonpT%d",nJ,mu),Form("#mu=%d",mu), "p_{T}(#mu) [GeV]","Events",75,0,1500);
			CreateHistogram(Form("jetpT%d_%dMu",nJ,mu),Form("Jet p_{T} spectrum for %dth jet",nJ),"Jet p_{T} [GeV]","Events",300,0,1500);			

			CreateHistogram(Form("st_nJets%d_%dMu_uW",nJ,mu),Form("%d jets, %d#mu",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			CreateHistogram(Form("st_nJets%d_%dMu",nJ,mu),Form("%d jets, %d#mu",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
		
			CreateHistogram(Form("st_nJets%d_%dMu_0btags",nJ,mu),Form("%d jets, %d#mu, 0 b-tags",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			CreateHistogram(Form("st_nJets%d_%dMu_1btags",nJ,mu),Form("%d jets, %d#mu, 1 b-tags",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			CreateHistogram(Form("st_nJets%d_%dMu_2btags",nJ,mu),Form("%d jets, %d#mu, #geq2 b-tags",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			
			CreateHistogram(Form("st_nJets%d_%dMu_0btags_SF",nJ,mu),Form("%d jets, %d#mu, 0 b-tags",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			CreateHistogram(Form("st_nJets%d_%dMu_1btags_SF",nJ,mu),Form("%d jets, %d#mu, 1 b-tags",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
			CreateHistogram(Form("st_nJets%d_%dMu_2btags_SF",nJ,mu),Form("%d jets, %d#mu, #geq2 b-tags",nJ,mu),"S_{T} [GeV]","Events",N,st_bins);
		
			if(nJ==5){
				nJ=56;
				CreateHistogram(Form("st_nJets%d_%dMu_0btags",nJ,mu),Form("%d jets, %d#mu, 0 b-tags",nJ,mu),"S_{T} [GeV]","Events",N2,st_bins2);
				CreateHistogram(Form("st_nJets%d_%dMu_1btags",nJ,mu),Form("%d jets, %d#mu, 1 b-tags",nJ,mu),"S_{T} [GeV]","Events",N2,st_bins2);
				CreateHistogram(Form("st_nJets%d_%dMu_2btags",nJ,mu),Form("%d jets, %d#mu, #geq2 b-tags",nJ,mu),"S_{T} [GeV]","Events",N2,st_bins2);
				
				CreateHistogram(Form("st_nJets%d_%dMu_0btags_SF",nJ,mu),Form("%d jets, %d#mu, 0 b-tags",nJ,mu),"S_{T} [GeV]","Events",N2,st_bins2);
				CreateHistogram(Form("st_nJets%d_%dMu_1btags_SF",nJ,mu),Form("%d jets, %d#mu, 1 b-tags",nJ,mu),"S_{T} [GeV]","Events",N2,st_bins2);
				CreateHistogram(Form("st_nJets%d_%dMu_2btags_SF",nJ,mu),Form("%d jets, %d#mu, #geq2 b-tags",nJ,mu),"S_{T} [GeV]","Events",N2,st_bins2);
				nJ=5;
			}
			
			
			//mt
			CreateHistogram(Form("mt_nJets%d_%dMu",nJ,mu),Form("M_{T} %d-Jets, #mu=%d",nJ,mu),"M_{T} [GeV]","Events",100,0,500);  
			CreateHistogram(Form("met_nJets%d_%dMu",nJ,mu),Form("M_{ET} %d-Jets, #mu=%d",nJ,mu),"M_{ET} [GeV]","Events",50,0,1000);  
			CreateHistogram(Form("metsig_nJets%d_%dMu",nJ,mu),Form("M_{ET} %d-Jets, #mu=%d",nJ,mu),"M_{ET}/H_{T} #sqrt{[GeV]}","Events",60,0,30);  

			
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
		hName["bin_width"]->SetBinError(i,0);
	}
	//cout << endl; 
	
	CreateHistogram("Muon_pT","Muon p_{T}","p_{T} GeV","Events",100,0,st_max); 
	
	CreateHistogram("eta_sT","eta vs sT", "S_{T} GeV", "|#eta|", 10, 0, 5, 100, 0, 3000); 
	
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
		
		/*
		cout << "pT Mu: " << pM.Perp() << " met: " << met << endl; 
		cout << "phiMu: " << pM.Phi() << " met phi: " << met_phi << endl; 
		
		*/
		
		mt=TMath::Sqrt(2*met*pM.Perp()*(1-TMath::Cos(pM.DeltaPhi(pNu)) ));
				
	}
	
}

float calcSt(){
    
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
	
	float sT=ht+leptpT+gammapT;
	if(met>15)sT+=met;
	
//	cout << "st: " << st << " calcst: "<< sT << endl; 
   

	return sT; 
}

void get_muTrigEff(){

	open_graph("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");
	
	open_graph("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
	open_graph("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");
    
    TH2F *hveto=(TH2F*)fEff_elec->FindObjectAny("sfVETO");
    hName2D[hveto->GetName()]=hveto;
    
    TH2F *hloose=(TH2F*)fEff_elec->FindObjectAny("sfLOOSE");
    hName2D[hloose->GetName()]=hloose;
    
    TH2F *hmedium=(TH2F*)fEff_elec->FindObjectAny("sfMEDIUM");
    hName2D[hmedium->GetName()]=hmedium;
    
    TH2F *htight=(TH2F*)fEff_elec->FindObjectAny("sfTIGHT");
    hName2D[htight->GetName()]=htight;

    TH1F* hx=(TH1F*)htight->ProjectionX("sfTIGHT_x");
    TH1F* hy=(TH1F*)htight->ProjectionY("sfTIGHT_y");
    
    hName[hx->GetName()]=hx;
    hName[hy->GetName()]=hy;
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

void load_btagEff(){
	TString file_name="/eos/uscms/store/user/btcarlso/efficiencies/efficiency_file_"+output_file_name+".root"; 
	if(test)file_name="/eos/uscms/store/user/btcarlso/efficiencies/efficiency_file_ttbar.root";
	fBtagEff = new TFile(file_name,"READ"); 
	
	TH1F *h=(TH1F*)fBtagEff->FindObjectAny("TEfficiency_names"); 
	
	for(int i=1; i<=h->GetNbinsX();i++){
		TString name=h->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names" || name=="TEfficiency_names") continue; 
		//cout << "loading file: " << name << endl; 
		TEfficiency *pEff=(TEfficiency*)fBtagEff->FindObjectAny(name); 
		effName[pEff->GetName()]=pEff; 
	}
	
}


void efficiency(){
	TLorentzVector pM;
    TLorentzVector pE;
	
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
	
    double SFele=1;
    
    for(int el=0; el<electron_px->size();el++){
        pE.SetPxPyPzE(electron_px->at(el),electron_py->at(el),electron_pz->at(el),electron_e->at(el));
        int ieta=hName["sfTIGHT_x"]->FindBin(TMath::Abs(pE.Eta()));
        int ipt=hName["sfTIGHT_y"]->FindBin(TMath::Abs(pE.Perp()));

        int Nx_=hName["sfTIGHT_x"]->GetNbinsX();
        int Ny_=hName["sfTIGHT_x"]->GetNbinsY();

        if(ieta>Nx_) ieta=Nx_;
        if(ieta>Ny_) ipt=Ny_;
        
        SFele=hName2D["sfTIGHT"]->GetBinContent(ieta,ipt);//tight id electron sf
        
    }
    
	double R_dataMC=(1-(1-EffData1)*(1-EffData2))/(1-(1-EffMC1)*(1-EffMC2)); 
	//cout << "R_data/MC: " << R_dataMC << endl; 
	if(_MC)prName["R_dataMC_eff"]->Fill(pt,R_dataMC); 
	//cout << "weight: " << weight << endl; 
	if(_MC)weight=weight*R_dataMC*SFele;
	//cout << "weight: " << weight << endl; 
}

void get_corrections(){
	TH1F *h=(TH1F*)fCorr->FindObjectAny("data_MC"); 
	h->SetName("Correction_DataMC"); 
	hName[h->GetName()]=h; 
	
}

void reweightMC(){
	if(_MC==false || Correct_==false) return; 
	TLorentzVector pJ; 
	pJ.SetPxPyPzE(jet_px->at(0),jet_py->at(0),jet_pz->at(0),jet_e->at(0)); 
	double pt=pJ.Perp();
	double Cor=hName["Correction_DataMC"]->Interpolate(pt); 
	weight=weight*Cor; 
	
	//cout	<< "weight: " << weight << " Cor: " << Cor << endl; 
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

void fill_pileup(){
	hName["h_npv"]->Fill(nPv,weight); 
	if(nJets==4)hName["h_npv_nJets4"]->Fill(nPv,weight); 

	prName["pr_nJets_npv"]->Fill(nPv, nJets,weight); 
	
}

void fill_categories(){
    
	int Nmu=muon_px->size();
	int NEl=electron_px->size();
    int charge=1;
    for(int imu=0; imu<muon_px->size(); imu++){
        charge*=muon_charge->at(imu);
    }
    for(int iel=0; iel<electron_charge->size(); iel++){
        charge*=electron_charge->at(iel);
    }
	
	int NB=nBJets;
	if(NB>2)NB=2;
	
	int NBT=nTightBJets;
	if(NBT>2)NBT=2;
    
    bool inclusive=false;
    
   // cout << "Data category: " << create_category(nJets,inclusive,Nmu,NEl,charge,NB,"M","_SF","off-z") << endl;
    

}

void fill_st(){
  
   
	bool print=false;
	if(print) cout << "fill st: " << endl; 
	
	int Nmu=muon_px->size();
	int NEl=electron_px->size();
	int nGamma=photon_px->size();
	
	int NB=nBJets; 
	if(NB>2)NB=2; 
	
	int NBT=nTightBJets; 
//	if(NBT>2)NBT=2;
    if(tt_enriched){
        hName[Form("nBJetsTight_1Mu_1El")]->Fill(nTightBJets,weight*B_weight_tight);
        int NNBB=NBT;
        if(NNBB>=3)NNBB=3;
        if(nJets>=2)hName[Form("deltaPt_1Mu_1El_nJets%d_%dbtag",nJets,NNBB)]->Fill(deltaPt,weight*B_weight_tight);
    }
    if(NBT>=2)NBT=2;
   // if(NBT>2) return;
    
	if(nLooseBJets==0) cout << "nLooseBJets: " << nLooseBJets << endl;
	
    hName2D["btag_weight_tight"]->Fill(nJets,NBT,B_weight_tight);
    hName2D["btag_weight"]->Fill(nJets,NBT,B_weight);
    prName2D["probMC_nJet_nB_tight"]->Fill(nJets,0.,probMC0);
    prName2D["probMC_nJet_nB_tight"]->Fill(nJets,1.,probMC1);
    prName2D["probMC_nJet_nB_tight"]->Fill(nJets,2.,probMC2);

    
//	cout << "st: " << st << " Nmu " << Nmu << " Nel " << NEl << " dimuon mass " << dimuon_mass;
//	cout << " W: " << W_enriched << " dimuon pair: " << dimuon_enriched << " Z " << Z_enriched << " tt " << tt_enriched << " tt ss: " << tt_enriched_ss << " sig1 " << sig_1mu << " sig2 " << sig_2mu << endl; 
	if(print) cout << "ttbar_enriched_btag: " << tt_enriched_btag << endl; 
	
	if(print) cout << "fill hmuel and hmumu: " << endl; 
	
	if(tt_enriched && dilepton_mass>0) hName[Form("h_muel_nJet%d",nJets)]->Fill(dilepton_mass,weight); 
	if(dimuon_mass>50 && dimuon_enriched) {
        hName["h_mumu"]->Fill(dimuon_mass,weight);
		hName[Form("h_mumu_nJet%d",nJets)]->Fill(dimuon_mass,weight);
		hName[Form("h_mumu_nJet%d_%dbtag",nJets,NB)]->Fill(dimuon_mass,weight*B_weight); 
		if(nJets>=5){
			int jtmp=56;
			if(st>300 && st<=700)hName[Form("h_mumu_nJet%d_%dbtag_st1",jtmp,NB)]->Fill(dimuon_mass,weight*B_weight); 
			if(st>700 && st<=3000)hName[Form("h_mumu_nJet%d_%dbtag_st2",jtmp,NB)]->Fill(dimuon_mass,weight*B_weight); 
		}
		hName[Form("h_mumu_nJet%d_0p5GeV",nJets)]->Fill(dimuon_mass,weight); 

	}
    if(dimuon_mass>50 && dimuon_enriched_ss)  {
        hName["h_ssmumu"]->Fill(dimuon_mass,weight);
        if(dimuon_mass<87 || dimuon_mass>95)hName[Form("nJets_ss%dMu_%dbtag_tight",Nmu,NBT)]->Fill(nJets,weight*B_weight_tight);
    }

    
	if(print) cout << "fill Z enriched: " << endl; 
	//fill 3 control regions for normalization 
	
	//fill Z

    if(Z_enriched_ss){
        for(int ijet=0; ijet<jet_px->size(); ijet++){
            TLorentzVector pJ;
            int jet=ijet+1;
            if(jet>nJetmax)jet=nJetmax;
            pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
            hName[Form("jetpT%d_ssZ",jet)]->Fill(pJ.Perp(),weight);
        }
        
        hName[Form("nJets_ssZ_%dbtag_tight_SFP",NBT)]->Fill(nJets,weight*B_weight_tightP);
        hName[Form("nJets_ssZ_%dbtag_tight_SFM",NBT)]->Fill(nJets,weight*B_weight_tightM);
        hName[Form("nJets_ssZ_%dbtag_tight_SF",NBT)]->Fill(nJets,weight*B_weight_tight);
    }
    
	if(Z_enriched){
		hName[Form("nJets_Z_%dbtag",NB)]->Fill(nJets,weight); 
		hName[Form("nJets_Z_%dbtag_tight",NBT)]->Fill(nJets,weight); 
		
		hName[Form("nJets_Z_%dbtag_SF",NB)]->Fill(nJets,weight*B_weight); 
		hName[Form("nJets_Z_%dbtag_SFP",NB)]->Fill(nJets,weight*B_weightP); // + sys 
		hName[Form("nJets_Z_%dbtag_SFM",NB)]->Fill(nJets,weight*B_weightM); //- sys

		hName[Form("nJets_Z_%dbtag_tight_SF",NBT)]->Fill(nJets,weight*B_weight_tight); 
		hName[Form("nJets_Z_%dbtag_tight_SFP",NBT)]->Fill(nJets,weight*B_weight_tightP); 
		hName[Form("nJets_Z_%dbtag_tight_SFM",NBT)]->Fill(nJets,weight*B_weight_tightM); 

		
		hName[Form("st_nJets%d_Z_%dbtag",nJets,NB)]->Fill(st,weight);
		hName[Form("st_nJets%d_Z_%dbtag_tight",nJets,NBT)]->Fill(st,weight);

		hName[Form("st_nJets%d_Z_%dbtag_SF",nJets,NB)]->Fill(st,weight*B_weight);
		hName[Form("st_nJets%d_Z_%dbtag_SFP",nJets,NB)]->Fill(st,weight*B_weightP);
		hName[Form("st_nJets%d_Z_%dbtag_SFM",nJets,NB)]->Fill(st,weight*B_weightM);

		hName[Form("st_nJets%d_Z_%dbtag_tight_SF",nJets,NBT)]->Fill(st,weight*B_weight_tight);
		hName[Form("st_nJets%d_Z_%dbtag_tight_SFP",nJets,NBT)]->Fill(st,weight*B_weight_tightP);
		hName[Form("st_nJets%d_Z_%dbtag_tight_SFM",nJets,NBT)]->Fill(st,weight*B_weight_tightM);

		if(nJets>=5){
			int jtmp=nJets;
			nJets=56;
			
			hName[Form("st_nJets%d_Z_%dbtag_SF",nJets,NB)]->Fill(st,weight*B_weight);
			hName[Form("st_nJets%d_Z_%dbtag_SFP",nJets,NB)]->Fill(st,weight*B_weightP);
			hName[Form("st_nJets%d_Z_%dbtag_SFM",nJets,NB)]->Fill(st,weight*B_weightM);
			
			hName[Form("st_nJets%d_Z_%dbtag_tight_SF",nJets,NBT)]->Fill(st,weight*B_weight_tight);
			hName[Form("st_nJets%d_Z_%dbtag_tight_SFP",nJets,NBT)]->Fill(st,weight*B_weight_tightP);
			hName[Form("st_nJets%d_Z_%dbtag_tight_SFM",nJets,NBT)]->Fill(st,weight*B_weight_tightM);	
			nJets=jtmp;
		}
		
		
		hName2D["st_nJets_Z"]->Fill(st,nJets,weight);
		hName["nJets_Z"]->Fill(nJets,weight); 
		hName[Form("st_nJets%d_Z_uW",nJets)]->Fill(st);
		hName[Form("st_nJets%d_Z",nJets)]->Fill(st,weight);
		hName[Form("ht_nJets%d_Z",nJets)]->Fill(ht,weight);
		hName[Form("met_nJets%d_Z",nJets)]->Fill(met,weight);
		hName[Form("metsig_nJets%d_Z",nJets)]->Fill(met/TMath::Sqrt(ht),weight);


		hName[Form("Z_nJets%d_dileptonpT",nJets)]->Fill(dilepton_pt,weight); 
		hName[Form("Z_nJets%d_muonpT%d",nJets,1)]->Fill(pTmuon1,weight);
		hName[Form("Z_nJets%d_muonpT%d",nJets,2)]->Fill(pTmuon2,weight);
				
		if(nJets>=1) hName[Form("ht_inclusive_nJets%d_Z",1)]->Fill(ht,weight);
		if(nJets>=2) hName[Form("ht_inclusive_nJets%d_Z",2)]->Fill(ht,weight);
		if(nJets>=3) hName[Form("ht_inclusive_nJets%d_Z",3)]->Fill(ht,weight);
		if(nJets>=4) hName[Form("ht_inclusive_nJets%d_Z",4)]->Fill(ht,weight);
	}
	
	//Fill W
	if(W_enriched){
		if(nBJets==0)hName["nJets_W_0btag"]->Fill(nJets,weight);
 	    if(nBJets==0) hName[Form("st_nJets%d_W_0btag",nJets)]->Fill(st,weight);
		hName[Form("st_nJets%d_W_0btag_SF",nJets)]->Fill(st,weight*B_weight);

		if(nBJets==1)hName["nJets_W_1btag"]->Fill(nJets,weight);
		if(nBJets==1) hName[Form("st_nJets%d_W_1btag",nJets)]->Fill(st,weight);
		hName[Form("st_nJets%d_W_1btag_SF",nJets)]->Fill(st,weight*B_weight);

		
		if(nBJets>=2)hName["nJets_W_2btag"]->Fill(nJets,weight);
		if(nBJets>=2) hName[Form("st_nJets%d_W_2btag",nJets)]->Fill(st,weight);
		hName[Form("st_nJets%d_W_2btag_SF",nJets)]->Fill(st,weight*B_weight);

		hName["nJets_W_0btag_SF"]->Fill(nJets,weight*B_weight);
		hName["nJets_W_1btag_SF"]->Fill(nJets,weight*B_weight);
		hName["nJets_W_2btag_SF"]->Fill(nJets,weight*B_weight);

		
		hName2D["st_nJets_W"]->Fill(st,nJets,weight);
		hName["nJets_W"]->Fill(nJets,weight); 
		hName[Form("st_nJets%d_W_uW",nJets)]->Fill(st);
		hName[Form("st_nJets%d_W",nJets)]->Fill(st,weight);
		hName[Form("ht_nJets%d_W",nJets)]->Fill(ht,weight);
		hName[Form("met_nJets%d_W",nJets)]->Fill(met,weight);
		hName[Form("metsig_nJets%d_W",nJets)]->Fill(met/TMath::Sqrt(ht),weight);


		hName[Form("W_nJets%d_muonpT",nJets)]->Fill(pTmuon1,weight);
		
		if(nJets>=1) hName[Form("ht_inclusive_nJets%d_W",1)]->Fill(ht,weight);
		if(nJets>=2) hName[Form("ht_inclusive_nJets%d_W",2)]->Fill(ht,weight);
		if(nJets>=3) hName[Form("ht_inclusive_nJets%d_W",3)]->Fill(ht,weight);
		if(nJets>=4) hName[Form("ht_inclusive_nJets%d_W",4)]->Fill(ht,weight);
	}
	if(print) cout <<"fill b-veto: " << endl; 

	
	if(print) cout << "fill ttbar: " << endl; 
	//fill ttbar
	if(tt_enriched_ss){
        hName[Form("met_nJets%d_ss1Mu_1El_%dbtag",nJets,NB)]->Fill(met,weight*B_weight);
		hName[Form("nJets_ss1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(nJets,weight*B_weight_tight);
		hName[Form("nJets_ss1Mu_1El_%dbtag_tight_SFP",NBT)]->Fill(nJets,weight*B_weight_tightP);
		hName[Form("nJets_ss1Mu_1El_%dbtag_tight_SFM",NBT)]->Fill(nJets,weight*B_weight_tightM);
		
		hName[Form("st_nJets%d_ss1Mu_1El_%dbtag_tight_SF",nJets,NBT)]->Fill(st,weight*B_weight_tight); 
		hName[Form("st_nJets%d_ss1Mu_1El_%dbtag_tight_SFP",nJets,NBT)]->Fill(st,weight*B_weight_tightP); 
		hName[Form("st_nJets%d_ss1Mu_1El_%dbtag_tight_SFM",nJets,NBT)]->Fill(st,weight*B_weight_tightM); 
	}
	
	if(tt_enriched){
		prName[Form("pr_corr_nJets%d_st",nJets)]->Fill(st,TTbar_corr);

		if(print) cout << "2d st-njets " << endl; 
		hName2D["st_nJets_tt"]->Fill(st,nJets,weight);
		if(print) cout << "1Mu/1El" << endl;
		
		hName["nJets_1Mu_1El"]->Fill(nJets,weight);
        
        hName[Form("weight_nJets%d_1Mu_1El_%dbtag_tight",nJets,NBT)]->Fill(B_weight_tight);
		
        hName[Form("st_1Mu_1El_%dbtag_tight_TopCor",NBT)]->Fill(st,weight*B_weight_tight*TTbar_corr);
        hName[Form("st_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(st,weight*B_weight_tight);

        hName[Form("st_1Mu_1El_%dbtag_TopCor",NBT)]->Fill(st,weight*B_weight*TTbar_corr);
        hName[Form("st_1Mu_1El_%dbtag_SF",NBT)]->Fill(st,weight*B_weight);
        
        hName[Form("met_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(met,weight*B_weight_tight);
        hName[Form("met_1Mu_1El_%dbtag_tight_TopCor",NBT)]->Fill(met,weight*B_weight_tight*TTbar_corr);
        
        hName[Form("bjet1pt_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(pTbjet1,weight*B_weight_tight);
        hName[Form("bjet2pt_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(pTbjet2,weight*B_weight_tight);
        
        hName[Form("bjet1pt_1Mu_1El_%dbtag_tight_TopCor",NBT)]->Fill(pTbjet1,weight*B_weight_tight*TTbar_corr);
        hName[Form("bjet2pt_1Mu_1El_%dbtag_tight_TopCor",NBT)]->Fill(pTbjet2,weight*B_weight_tight*TTbar_corr);

        hName[Form("ptmu_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(pTmuon1,weight*B_weight_tight);
        hName[Form("ptmu_1Mu_1El_%dbtag_tight_TopCor",NBT)]->Fill(pTmuon1,weight*B_weight_tight*TTbar_corr);

        hName[Form("pte_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(pTElectron1,weight*B_weight_tight);
        hName[Form("pte_1Mu_1El_%dbtag_tight_TopCor",NBT)]->Fill(pTElectron1,weight*B_weight_tight*TTbar_corr);
        
        
		hName[Form("st_nJets%d_1Mu_1El_%dbtag",nJets,NB)]->Fill(st,weight);//first make plots without scale factors
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight",nJets,NBT)]->Fill(st,weight); 	
		
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_SF",nJets,NB)]->Fill(st,weight*B_weight); 
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_SFP",nJets,NB)]->Fill(st,weight*B_weightP); 
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_SFM",nJets,NB)]->Fill(st,weight*B_weightM); 

		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFP_topreweighted",nJets,NBT)]->Fill(st,weight*B_weight_tightP*TTbar_corr); 
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFM_topreweighted",nJets,NBT)]->Fill(st,weight*B_weight_tightM*TTbar_corr); 
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SF_topreweighted",nJets,NBT)]->Fill(st,weight*B_weight_tight*TTbar_corr); 

        hName[Form("stTop_nJets%d_1Mu_1El_%dbtag_tight_SF",nJets,NBT)]->Fill(st_top,weight*B_weight_tight);
        
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SF",nJets,NBT)]->Fill(st,weight*B_weight_tight); 
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFP",nJets,NBT)]->Fill(st,weight*B_weight_tightP); 
		hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFM",nJets,NBT)]->Fill(st,weight*B_weight_tightM); 
		
		if(nJets>=5){
			int jtmp=nJets;
			nJets=56;
			hName[Form("st_nJets%d_1Mu_1El_%dbtag_SF",nJets,NB)]->Fill(st,weight*B_weight); 
			hName[Form("st_nJets%d_1Mu_1El_%dbtag_SFP",nJets,NB)]->Fill(st,weight*B_weightP); 
			hName[Form("st_nJets%d_1Mu_1El_%dbtag_SFM",nJets,NB)]->Fill(st,weight*B_weightM); 
			
			hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SF",nJets,NBT)]->Fill(st,weight*B_weight_tight); 
			hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFP",nJets,NBT)]->Fill(st,weight*B_weight_tightP); 
			hName[Form("st_nJets%d_1Mu_1El_%dbtag_tight_SFM",nJets,NBT)]->Fill(st,weight*B_weight_tightM); 
			nJets=jtmp;
		}
		

		hName[Form("nJets_1Mu_1El_%dbtag",NB)]->Fill(nJets,weight);
        int NBTRUTH=nb_truth;
        if(NBTRUTH>=2)NBTRUTH=2;
        hName2D[Form("Acceptance_num_%db_st_nJets_1Mu_1El",NBTRUTH)]->Fill(st,nJets,weight);
        hName2D["Acceptance_den_st_nJets_1Mu_1El"]->Fill(st,nJets,weight);

        hName[Form("Acceptance_num_%db_nJets_1Mu_1El",NBTRUTH)]->Fill(nJets,weight);
        hName["Acceptance_den_nJets_1Mu_1El"]->Fill(nJets,weight);

        hName[Form("nJets_1Mu_1El_nbtruth%d_%dbtag_tight_SF",NBTRUTH,NBT)]->Fill(nJets,weight*B_weight_tight);
		hName[Form("nJets_1Mu_1El_%dbtag_tight",NBT)]->Fill(nJets,weight);
        prName[Form("pr_weight_nJets_1Mu_1El_%dbtag_tight",NBT)]->Fill(nJets,B_weight_tight);
        prName[Form("pr_probMC_nJets_1Mu_1El_%dbtag_tight",0)]->Fill(nJets, probMC0);
        prName[Form("pr_probMC_nJets_1Mu_1El_%dbtag_tight",1)]->Fill(nJets, probMC1);
        prName[Form("pr_probMC_nJets_1Mu_1El_%dbtag_tight",2)]->Fill(nJets, probMC2);

        hName[Form("h_probMC_nJets_1Mu_1El_%dbtag_tight",0)]->Fill(nJets, probMC0);
        hName[Form("h_probMC_nJets_1Mu_1El_%dbtag_tight",1)]->Fill(nJets, probMC1);
        hName[Form("h_probMC_nJets_1Mu_1El_%dbtag_tight",2)]->Fill(nJets, probMC2);

        
		hName[Form("nJets_1Mu_1El_%dbtag_SF",NB)]->Fill(nJets,weight*B_weight);	
		hName[Form("nJets_1Mu_1El_%dbtag_SFP",NB)]->Fill(nJets,weight*B_weightP);		
		hName[Form("nJets_1Mu_1El_%dbtag_SFM",NB)]->Fill(nJets,weight*B_weightM);
        
		hName[Form("nJets_1Mu_1El_%dbtag_tight_SF",NBT)]->Fill(nJets,weight*B_weight_tight);
		hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP",NBT)]->Fill(nJets,weight*B_weight_tightP);
		hName[Form("nJets_1Mu_1El_%dbtag_tight_SFM",NBT)]->Fill(nJets,weight*B_weight_tightM);
       
        hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP_bc",NBT)]->Fill(nJets,weight*B_weight_tightPbc);
        hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP_l",NBT)]->Fill(nJets,weight*B_weight_tightPl);

		hName[Form("nJets_1Mu_1El_%dbtag_tight_SFM_bc",NBT)]->Fill(nJets,weight*B_weight_tightMbc);
		hName[Form("nJets_1Mu_1El_%dbtag_tight_SFM_l",NBT)]->Fill(nJets,weight*B_weight_tightMl);

		
		if(print) cout << "st 1Mu/1El" << endl; 
		
		hName[Form("st_nJets%d_1Mu_1El_uW",nJets)]->Fill(st); 
		hName[Form("st_nJets%d_1Mu_1El",nJets)]->Fill(st,weight); 
		if(print) cout << "ht " << endl; 
		hName[Form("ht_nJets%d_1Mu_1El",nJets)]->Fill(ht,weight); 
		if(print) cout	 << "met" << endl; 
		hName[Form("met_nJets%d_1Mu_1El",nJets)]->Fill(met,weight);
		hName[Form("metsig_nJets%d_1Mu_1El",nJets)]->Fill(met/TMath::Sqrt(ht),weight);

		if(print) cout << "dilepton pt, pt1,pt2 " << endl; 
		
		hName[Form("tt_nJets%d_dileptonpT",nJets)]->Fill(dilepton_pt,weight); 
		hName[Form("tt_nJets%d_muonpT%d",nJets,1)]->Fill(pTmuon1,weight);
		hName[Form("tt_nJets%d_muonpT%d",nJets,2)]->Fill(pTmuon2,weight);
		
		if(nJets>=1) hName[Form("ht_inclusive_nJets%d_1Mu_1El",1)]->Fill(ht,weight);
		if(nJets>=2) hName[Form("ht_inclusive_nJets%d_1Mu_1El",2)]->Fill(ht,weight);
		if(nJets>=3) hName[Form("ht_inclusive_nJets%d_1Mu_1El",3)]->Fill(ht,weight);
		if(nJets>=4) hName[Form("ht_inclusive_nJets%d_1Mu_1El",4)]->Fill(ht,weight);
	}

	
	if(tt_enriched_btag){
		if(print) cout << "nJets dist." << endl; 
		if(print) cout << "st distribution" << endl; 
		prName[Form("pr_corr_nJets%d_st_1Mu_1El_btag",nJets)]->Fill(st,TTbar_corr); 

	}
	
	//fill signal region
	if(Nmu==1 && sig_1mu==false) return; 
	if(Nmu==2 && sig_2mu==false) return;
	
	if(print)cout << "fill main njets histograms: " << endl; 
	hName[Form("nJets_%dMu",Nmu)]->Fill(nJets,weight);
    hName[Form("nJets_%dMu_%dbtag",Nmu,NB)]->Fill(nJets,weight);

	hName[Form("nJets_%dMu_%dbtag_tight",Nmu,NBT)]->Fill(nJets,weight*B_weight_tight);
    
	hName[Form("nJets_%dMu_%dbtag_SF",Nmu,NB)]->Fill(nJets,weight*B_weight);

	hName[Form("nBJetsLoose_%dMu",Nmu)]->Fill(nLooseBJets,weight); 
	hName[Form("nBJetsMedium_%dMu",Nmu)]->Fill(nBJets,weight); 
	hName[Form("nBJetsTight_%dMu",Nmu)]->Fill(nTightBJets,weight); 
	hName2D["h_Loose_Medium"]->Fill(nBJets,nLooseBJets,weight);
	hName2D["h_Medium_Tight"]->Fill(nTightBJets,nBJets,weight);
	hName2D["h_Loose_Tight"]->Fill(nTightBJets,nLooseBJets,weight);


	if(print)cout << "fill main st histograms: " << endl; 

	if(QCD)return; 
	if(print)cout << "muon binned: " << endl; 
	
	hName2D[Form("st_nJets_%dMu",Nmu)]->Fill(st,nJets,weight);
	
	hName[Form("st_nJets%d_%dMu_uW",nJets,Nmu)]->Fill(st);
	hName[Form("st_nJets%d_%dMu",nJets,Nmu)]->Fill(st,weight);
	hName[Form("st_nJets%d_%dMu_%dbtags",nJets,Nmu,NB)]->Fill(st,weight);
		
	hName[Form("st_nJets%d_%dMu_%dbtags_SF",nJets,Nmu,NB)]->Fill(st,weight*B_weight);
	
	if(nJets>=5){
		int jtmp=nJets;
		nJets=56;
		if(nBJets==0)hName[Form("st_nJets%d_%dMu_0btags_SF",nJets,Nmu)]->Fill(st,weight*B_weight); 
		if(nBJets==1)hName[Form("st_nJets%d_%dMu_1btags_SF",nJets,Nmu)]->Fill(st,weight*B_weight); 
		if(nBJets>=2)hName[Form("st_nJets%d_%dMu_2btags_SF",nJets,Nmu)]->Fill(st,weight*B_weight); 
		nJets=jtmp;
	}
	
	hName[Form("st_nJets%d_%dMu_uW_50GeV",nJets,Nmu)]->Fill(st);
	hName[Form("st_nJets%d_%dMu_50GeV",nJets,Nmu)]->Fill(st,weight);
	
	if(print) cout << "fill met " << endl; 
	hName[Form("met_nJets%d_%dMu",nJets,Nmu)]->Fill(met,weight); 
	hName[Form("metsig_nJets%d_%dMu",nJets,Nmu)]->Fill(met/TMath::Sqrt(ht),weight); 

	hName[Form("mt_nJets%d_%dMu",nJets,Nmu)]->Fill(mt,weight); 
		
	if(print) cout << "end of fill_st() " << endl; 
    

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
		if(jet_bTagM->at(ijet)){
			//fill histogram for b,mu mass 
			pJ+=pM; 
			double BMuMass=pJ.M();
			hName[Form("bjet%d_muMass",jet)]->Fill(BMuMass,weight*B_weight);
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
/*
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
*/
		
	}	
}

void nBtag(){
	bool print=false;
	//if(muon_px->size()>=1 && nJets>=2) print=true;
	
	if(print) cout << "nJets: "<< nJets << endl; 
	
	int nB=0; 
	int nBL=0;
	int nBT=0;
    
    nb_truth=0;
    nc_truth=0;
    nl_truth=0;
    n_split=0;
    
    bool flv_disc=false;
    std::vector<TLorentzVector> btagged;
    TLorentzVector pJi;
	for(int ib=0; ib<jet_bTagM->size();ib++){
        pJi.SetPxPyPzE(jet_px->at(ib),jet_py->at(ib),jet_pz->at(ib),jet_e->at(ib));

		bool bT=jet_bTagM->at(ib);
		bool bTloose = jet_bTagL->at(ib); 
		bool bTtight = jet_bTagT->at(ib);
        if(bTtight) btagged.push_back(pJi);

		if(print)cout << "loose tag: " << bTloose << " tight tag: "<< bT << endl;
		
		if(bT==1) nB++;
		if(bTloose) nBL++; 
		if(bTtight) nBT++;
        int flv=TMath::Abs(jet_algFlavor->at(ib));
        int phyflv=TMath::Abs(jet_phyFlavor->at(ib));
        if(flv==5)nb_truth++;
        if(flv==4)nc_truth++;
        if(flv==21)nl_truth++;
        if(flv<=3)nl_truth++;
        if(phyflv==21 && flv!=21)n_split++;
        
	}
    
    for(int ib=0; ib<btagged.size();ib++){
        pJi.SetPxPyPzE(jet_px->at(ib),jet_py->at(ib),jet_pz->at(ib),jet_e->at(ib));
        pTbjet1=pJi.Perp();
        pTbjet2=pJi.Perp();
        
        for(int jb=ib+1; jb<btagged.size();jb++){
            if(ib==jb)continue;
            TLorentzVector pJj;
            pJj.SetPxPyPzE(jet_px->at(jb),jet_py->at(jb),jet_pz->at(jb),jet_e->at(jb));
            double theta=pJi.Angle(pJj.Vect());
            pJj+=pJi;
            TLorentzVector pJj_boost=pJj;
            pJj_boost.Boost(-pJj.BoostVector());
            pJi.Boost(-pJj.BoostVector());
            double deltaPhi=TMath::Abs(pJj_boost.DeltaPhi(pJi));
            double E=pJj.E();

            hName2D["theta_E"]->Fill(E,theta,weight);
            hName2D["theta_M"]->Fill(pJj.M(),theta,weight);

            hName2D["deltaPhi_M"]->Fill(pJj.M(),deltaPhi,weight);
            hName["M_bb"]->Fill(pJj.M(),weight);
            //cout << "deltaPhi: " << deltaPhi << " M: " << pJj.M() << endl;
            //cout << "theta: "<< ib << " " << jb << " "  << theta << " M " << pJj.M() << endl;
        }
    }
    
    nBJets=nB;
	nLooseBJets=nBL;
	nTightBJets=nBT;
    
    if(nBT>=3 && print){
    for(int ib=0; ib<jet_bTagM->size();ib++){
        pJi.SetPxPyPzE(jet_px->at(ib),jet_py->at(ib),jet_pz->at(ib),jet_e->at(ib));
       /* cout << "pt: " << pJi.Perp() << " " << pJi.Eta() << " " << pJi.Phi()<< endl;
        cout << "jet phy flavor: " << jet_phyFlavor->at(ib) << endl;
        cout << "jet alg flavor: " << jet_algFlavor->at(ib) << endl;
        cout << "tag: " << jet_bTagT->at(ib) << endl;*/
    }
        cout << "nb: " << nb_truth << " nc: " << nc_truth << " nl: "<< nl_truth << endl;
        cout << "nBTight: " << nTightBJets  << endl;
        cout << "nJets: " << nJets << endl << endl;
        cout << "n_split : "<< n_split << endl;
    }

	if(print)cout << "nBJets: " << nBJets << endl; 
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
		//hName[Form("nJets%d_muonpT%d",nJets,mu)]->Fill(pM.Perp(),weight);
	}
}

void gen_pt(){
   

	TTbar_corr=1;
    TTbar_corrP=1;
	TTbar_corrM=1;

	if(!_MC) return;
    for(int i=0; i<w_px->size(); i++){
        wP.SetPxPyPzE(w_px->at(i), w_py->at(i), w_pz->at(i), w_e->at(i));
       // cout << "wPt: " << wP.Perp() << " eta " << wP.Eta() << " phi: " << wP.Phi() << endl;
    }
	ttP.SetPxPyPzE(0,0,0,0);
	for(int i=0; i<top_px->size(); i++){
		TLorentzVector top; 
		top.SetPxPyPzE(top_px->at(i), top_py->at(i), top_pz->at(i), top_e->at(i));
      /*
        cout << "top pt: " << top.Perp() << " eta: " << top.Eta() << " phi " << top.Phi() << endl;
        for(int i=0; i<w_px->size(); i++){
            wP.SetPxPyPzE(w_px->at(i), w_py->at(i), w_pz->at(i), w_e->at(i));
            for(int ij=0; ij<jet_px->size();ij++){
                TLorentzVector pJ;
                pJ.SetPxPyPzE(jet_px->at(ij),jet_py->at(ij),jet_pz->at(ij),jet_e->at(ij));
                pJ+=wP;
                double dR=top.DeltaR(pJ);
                double dPT= top.Perp()-pJ.Perp();
                double dRPT=TMath::Sqrt(dR*dR+0.0001*dPT*dPT);
                cout << "deltaR: " << dR << " delta pT: " << dPT<< " drpt " << dRPT << " flavor " << jet_algFlavor->at(ij)<< endl;
                //cout << "wPt: " << wP.Perp() << " eta " << wP.Eta() << " phi: " << wP.Phi() << endl;
            }
        }
        cout << endl;
       */
		hName["gen_top_pt"]->Fill(top.Perp());
        int NBTRUTH=nb_truth;
        if(NBTRUTH>2)NBTRUTH=2;
        
        hName[Form("gen_top_pt_%db",NBTRUTH)]->Fill(top.Perp());

		double	 pt_top=top.Perp();
		if(pt_top>400)pt_top=400;
		TTbar_corr=TTbar_corr*expo_corr->Eval(pt_top);
        TTbar_corrP=1.2*TTbar_corrP*expo_corr->Eval(pt_top);
		TTbar_corrM=0.8*TTbar_corrM*expo_corr->Eval(pt_top);

		hName["gen_top_pt_reweighted"]->Fill(top.Perp(),TTbar_corr);
        hName[Form("gen_top_pt_%db_reweighted",NBTRUTH)]->Fill(top.Perp(),TTbar_corr);

	//	cout << "top pT: " << top.Perp() << endl;
		ttP+=top; 
	}
    TTbar_corr=TMath::Sqrt(TTbar_corr)/0.99; //divide by average effect
    TTbar_corrP=TMath::Sqrt(TTbar_corrP)/1.2;
    TTbar_corrM=TMath::Sqrt(TTbar_corrM)/0.8;

	hName["top_ptweight_avg"]->Fill(TTbar_corr);

	//cout << "cor: "<< TTbar_corr << endl; 

	
	zP.SetPxPyPzE(0,0,0,0);
	for(int i=0; i<z_px->size(); i++){
		TLorentzVector Z; 
		Z.SetPxPyPzE(z_px->at(i), z_py->at(i), z_pz->at(i), z_e->at(i)); 
		zP+=Z; 
	}

	double wpT=wP.Perp(); 
	double ttpT=ttP.Perp();
	double zpT=zP.Perp();
	double corr=wpT-jetPsys.Perp(); 
	
	
	hName["h_wpt"]->Fill(wP.Perp()); 
	hName["h_wrapidity"]->Fill(wP.Rapidity()); 
	
	prName["pr_deltapt_wpt"]->Fill(wpT,wpT-jetPsys.Perp());
	prName["pr_deltapt_ttpt"]->Fill(ttpT,ttpT-jetPttsys.Perp());
	prName["pr_deltapt_zpt"]->Fill(zpT,zpT-jetPsys.Perp());
	
	prName["pr_deltapt_deltaphi"]->Fill(wP.DeltaPhi(jetPsys),wpT-jetPsys.Perp()); 
	//cout << "deltaphi: " << wP.DeltaPhi(pJ) << endl; 
	prName["pr_deltaphi_wpt"]->Fill(wpT,wP.DeltaPhi(jetPsys)); 
	hName["h_deltaphi"]->Fill(wP.DeltaPhi(jetPsys)); 
	hName["h_deltapt"]->Fill(wpT-jetPsys.Perp());
	hName2D["h_wpt_jpt"]->Fill(wpT, jetPsys.Perp());
    

}

bool count_electrons(){
	TLorentzVector pLoose; 
	TLorentzVector pTight; 
	
	nLooseElectrons=0; 
	nElectrons=electron_px->size();
	
	int nMatched=0; 
	std::vector<int> matched; 
	
	for(int iloose=0; iloose<loose_electron_px->size();iloose++){
		pLoose.SetPxPyPzE(loose_electron_px->at(iloose),loose_electron_py->at(iloose),loose_electron_pz->at(iloose),loose_electron_e->at(iloose)); 
		for(int itight=0; itight<electron_px->size();itight++){
			pTight.SetPxPyPzE(electron_px->at(itight),electron_py->at(itight),electron_pz->at(itight),electron_e->at(itight)); 
			bool skip=false; 
			for(int ii=0; ii<matched.size();ii++){
				if(itight==matched[ii]) skip=true; 
			}
			if(skip) continue; 
			
			if(pLoose.DeltaR(pTight)<0.01){
				matched.push_back(itight);
				nMatched++; 
			}
		}
	}
	nLooseElectrons=loose_electron_px->size()-nMatched; 
}

bool count_muons(){
	TLorentzVector pLoose; 
	TLorentzVector pTight; 
	
	nLooseMuons=0; 
	nMuons=muon_px->size(); 
	
	int nMatched=0; 
	std::vector<int> matched; 
	
	for(int iloose=0; iloose<loose_muon_px->size();iloose++){
		pLoose.SetPxPyPzE(loose_muon_px->at(iloose),loose_muon_py->at(iloose),loose_muon_pz->at(iloose),loose_muon_e->at(iloose)); 
		for(int itight=0; itight<muon_px->size();itight++){
			pTight.SetPxPyPzE(muon_px->at(itight),muon_py->at(itight),muon_pz->at(itight),muon_e->at(itight)); 
			bool skip=false; 
			for(int ii=0; ii<matched.size();ii++){
				if(itight==matched[ii]) skip=true; 
			}
			if(skip) continue; 
			
			if(pLoose.DeltaR(pTight)<0.01){
				matched.push_back(itight);
				nMatched++;
			}
		}
	}
	nLooseMuons=loose_muon_px->size()-nMatched;
}

void control_regions(){
    
	bool print=false;
	if(nLooseMuons>=1) return; 
	if(nLooseElectrons>=1) return; 
	//&& mt>75 && mt<250
	if(nMuons==1 && nElectrons==0){
		W_enriched=true; 
		if(nBJets==0)W_enriched_bveto=true;
	}
	//change 75-105 to 87-95
    if(nMuons==2 && nElectrons==0 && muon_charge->at(0)*muon_charge->at(1)==1){
        if(dimuon_mass>50)dimuon_enriched_ss=true;
        if(dimuon_mass>87 && dimuon_mass<95)Z_enriched_ss=true;
    }
    
	if(nMuons==2 && nElectrons==0 && muon_charge->at(0)*muon_charge->at(1)==-1){
		if(dimuon_mass>50) dimuon_enriched=true;
		if(dimuon_mass>81 && dimuon_mass<101) Z_enriched_loose=true;
		if( dimuon_mass>87 && dimuon_mass<95 ) {
			Z_enriched=true;
			if(nBJets==0)Z_enriched_bveto=true;
		}
		if( st>300  && ((dimuon_mass>50 && dimuon_mass<81) || dimuon_mass>101) ) sig_2mu=true;
	}
	if(nMuons==1 && nElectrons==1){
		if(muon_charge->at(0)*electron_charge->at(0)==-1){
			tt_enriched=true; 
			if(nBJets>=1)tt_enriched_btag=true;
		}
		else tt_enriched_ss=true;
	}
	if(st>300 && nMuons==1 && nElectrons==0 ) sig_1mu=true;
	
	TLorentzVector pM2;
	TLorentzVector pM1;

	if(Z_enriched){
		pM1.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0)); 
		pM2.SetPxPyPzE(muon_px->at(1),muon_py->at(1),muon_pz->at(1),muon_e->at(1)); 
		pM1+=pM2;
		pTmuon2=pM2.Perp();
		dilepton_pt=pM1.Perp();
	}
	if(tt_enriched || tt_enriched_btag){
		pM1.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0)); 
		pM2.SetPxPyPzE(electron_px->at(0),electron_py->at(0),electron_pz->at(0),electron_e->at(0));
        pTElectron1=pM2.Perp();
		pM1+=pM2;
		pTmuon2=pM2.Perp();
		dilepton_pt=pM1.Perp();
		dilepton_mass=pM1.Perp();
	}
	
	TLorentzVector pMtmp;
	if(print){
		cout << "Loose Muons: " << endl; 
		for(int im=0; im<loose_muon_px->size();im++){
			pMtmp.SetPxPyPzE(loose_muon_px->at(im),loose_muon_py->at(im),loose_muon_pz->at(im),loose_muon_e->at(im)); 
			cout << "pt: " << pMtmp.Perp() << " eta: " << pMtmp.Eta() << " phi: " << pMtmp.Phi() << endl; 
			
		}
		
		cout << "Tight muons: " << endl; 
		
		for(int im=0; im<muon_px->size();im++){
			pMtmp.SetPxPyPzE(muon_px->at(im),muon_py->at(im),muon_pz->at(im),muon_e->at(im)); 
			cout << "pt: " << pMtmp.Perp() << " eta: " << pMtmp.Eta() << " phi: " << pMtmp.Phi() << endl; 
			
		}
		cout << endl;
	}
	
	if(print){
		cout << "Loose Electrons: " << endl; 
		for(int im=0; im<loose_electron_px->size();im++){
			pMtmp.SetPxPyPzE(loose_electron_px->at(im),loose_electron_py->at(im),loose_electron_pz->at(im),loose_electron_e->at(im)); 
			cout << "pt: " << pMtmp.Perp() << " eta: " << pMtmp.Eta() << " phi: " << pMtmp.Phi() << endl; 
			
		}
		
		cout << "Tight Electrons: " << endl; 
		
		for(int im=0; im<electron_px->size();im++){
			pMtmp.SetPxPyPzE(electron_px->at(im),electron_py->at(im),electron_pz->at(im),electron_e->at(im)); 
			cout << "pt: " << pMtmp.Perp() << " eta: " << pMtmp.Eta() << " phi: " << pMtmp.Phi() << endl; 
			
		}
		
	}
	if(print){
		cout << "nLoose muons: " << nLooseMuons << " nMuons: " << nMuons << endl; 
		cout << "nLoose Electrons: " << nLooseElectrons << " nElectrons: " << nElectrons << endl; 
		cout << endl; 
	}
	

	
}
double SF_light(TString Atagger, TString mode, TLorentzVector jetP){
    //stopwatch["SF_light"].Start(kFALSE);
	double eta=TMath::Abs(jetP.Eta()); 
	double x=jetP.Perp();
	if(x>850)x=850; 
	double SF=0; 
	if(Atagger=="CSVM" && eta>=0 && eta<=0.8 ){
		double mean=((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
		if(mode=="mean")SF=((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x))); 
		if(mode=="min")SF=mean-((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
		if(mode=="max")SF=((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)))-mean;
	}
	
	if(Atagger=="CSVM" && eta>0.8 && eta<=1.6 ){
		double mean=((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
		if(mode=="mean")SF=((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
		if(mode=="min")SF=mean-((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
		if(mode=="max")SF=((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)))-mean;
	}
	
	
	if(Atagger=="CSVM" && eta>1.6 && eta<=2.4 ){
		double mean=((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
		if(mode=="mean")SF=((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
		if(mode=="min")SF=mean-((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
		if(mode=="max")SF=((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)))-mean; 
	}
	if(Atagger=="CSVT"){
		double mean=((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
		if(mode=="mean")SF=((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
		if(mode=="min")SF=mean-((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
		if(mode=="max")SF=((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)))-mean;
	}
	
	if((mode=="min" || mode=="max") && x>850) SF=SF*2; //double systematic for high pt
    //stopwatch["SF_light"].Stop();

	return SF;
}

void load_btag_sys(){
	float ptbin_b[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600,800};
	int Nb=sizeof(ptbin_b)/sizeof(double)-1;
	
	CreateHistogram("btag_systematic_medium","","jet p_{T} [GeV]","",Nb,ptbin_b); 
	CreateHistogram("btag_systematic_tight","","jet p_{T} [GeV]","",Nb,ptbin_b); 

	
	double SFb_error_M[] = {
		0.0415707,
		0.0204209,
		0.0223227,
		0.0206655,
		0.0199325,
		0.0174121,
		0.0202332,
		0.0182446,
		0.0159777,
		0.0218531,
		0.0204688,
		0.0265191,
		0.0313175,
		0.0415417,
		0.0740446,
		0.0596716 };
	
	double SFb_error_T[] = {
		0.0515703,
		0.0264008,
		0.0272757,
		0.0275565,
		0.0248745,
		0.0218456,
		0.0253845,
		0.0239588,
		0.0271791,
		0.0273912,
		0.0379822,
		0.0411624,
		0.0786307,
		0.0866832,
		0.0942053,
		0.102403 };
	
	for(int i=1; i<=hName["btag_systematic_medium"]->GetNbinsX();i++){
		hName["btag_systematic_medium"]->SetBinContent(i,SFb_error_M[i-1]); 
		hName["btag_systematic_tight"]->SetBinContent(i,SFb_error_T[i-1]); 

	}
	
	
}

double transfer_factor(TString Atagger){
   // cout << "transfer factor: " << endl;
    double p0=1;
    double T10=0;

    double p0PL=1;
    double p0ML=1;

    double p0PH=1;
    double p0MH=1;
    
    double T10PH=0;
    double T10MH=0;
    
    double T10PL=0;
    double T10ML=0;
    
    
    for(int ijet=0; ijet<jet_px->size(); ijet++){
        bool bc=false;
        int flv=TMath::Abs(jet_algFlavor->at(ijet));
        if(flv==4 || flv==5) bc=true;
        TLorentzVector pJi;
        pJi.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
        double s=getSF(ijet,Atagger,"");
        double sP=getSF(ijet,Atagger,"max");
        double sM=getSF(ijet,Atagger,"min");

        double eff=read_btag_efficiency(pJi,flv ,"T");
        double x=s*eff;
        double xP=sP*eff;
        double xM=sM*eff;

       // cout << "x: " << x << endl;
        p0*=1-x;
        if(bc) {
            p0PH*=1-xP;
            T10PH+=xP/(1-xP);

            p0PL*=1-x;
            T10PL+=x/(1-x);
        }
        else {
            p0PH*=1-x;
            T10PH+=x/(1-x);
            
            p0PL*=1-xP;
            T10PH+=xP/(1-xP);
        }
        
        T10+=x/(1-x);

    }
    //T10=1/sum(...)
    T10=1/T10;
    T10PH=1/T10PH;
    T10PL=1/T10PL;
    
    double p0P=TMath::Sqrt(TMath::Power(p0PH-p0,2)+TMath::Power(p0PL-p0,2));
    double T10P=TMath::Sqrt(TMath::Power(T10PH-T10,2)+TMath::Power(T10PL-T10,2));

    //cout << endl;
    double p1=p0/T10;
    double p2=1-p0-p1;
    double T20=1/((1/p0)-(1+(1/T10)));
    double T20PH=1/((1/p0PH)-(1+(1/T10PH)));
    double T20PL=1/((1/p0PL)-(1+(1/T10PL)));
    double T20P=TMath::Sqrt(TMath::Power(T20PH-T20,2)+TMath::Power(T20PL-T20,2));


    //now evaulate systematic
    /*
    double T20_sys=0;
    double T20_syscov=0;
    
    for(int ijet=0; ijet<jet_px->size(); ijet++){
        TLorentzVector pJi;
        pJi.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
        double s=getSF(ijet,Atagger,"");
        double eff=read_btag_efficiency(pJi,jet_algFlavor->at(ijet) ,"T");
        double x=s*eff;
        
        double sigma_s=getSF(ijet,Atagger,"")-getSF(ijet,Atagger,"min");
        double dT20_S=T20*T20*(-TMath::Power(p0,-2)*p0/(1-x)+(1-2*x)/TMath::Power(1-x,2))*eff;
        T20_sys+=TMath::Power(dT20_S*sigma_s,2);
      
        for(int kjet=0; kjet<jet_px->size();kjet++){
            if(ijet==kjet) continue;
            bool bc=false;
            if(jet_algFlavor->at(ijet)==4 || jet_algFlavor->at(ijet)==5)bc=true;
            bool bc_k=false;
            if(jet_algFlavor->at(kjet)==4 || jet_algFlavor->at(kjet)==5)bc_k=true;
            if(bc && !bc_k) continue; // jet i is b,c and jet k is light
            if(!bc && bc_k) continue; //jet i is light and jet k is b,c
            //only evaluate term when both jets are b,c or both jets are light.
            // that means no heavy/light correlation.
            pJi.SetPxPyPzE(jet_px->at(kjet),jet_py->at(kjet),jet_pz->at(kjet),jet_e->at(kjet));
             s=getSF(kjet,Atagger,"");
             eff=read_btag_efficiency(pJi,jet_algFlavor->at(kjet) ,"T");
             x=s*eff;
            double sigma_sk=getSF(ijet,Atagger,"")-getSF(ijet,Atagger,"min");
            double dT20_Sk=T20*T20*(-TMath::Power(p0,-2)*p0/(1-x)+(1-2*x)/TMath::Power(1-x,2))*eff;
            T20_sys+=dT20_S*dT20_Sk*sigma_s*sigma_sk;
            T20_syscov+=dT20_S*dT20_Sk*sigma_s*sigma_sk;
        }//k loop
        
    }//i loop
    T20_sys=TMath::Sqrt(T20_sys); //take sqrt for actual systematic;
    
    cout << "nJets: " << nJets << " tight tags: " << nTightBJets << " nb: " << nb_truth << " nc " << nc_truth << " nl " << nl_truth << endl;
    cout << "p0: " << p0 << " p1: " << p1 << " p2: "<< p2 << endl;
    cout << "T20: " << T20 << " sys " << T20P << endl;
    */
    int NBTRUTH=nb_truth;
    if(nb_truth>2)NBTRUTH=2;
  //  cout << "nJets: " << nJets << " nb " << NBTRUTH << endl;
    hName[Form("h_T20_nJets%d_nbtruth%d",nJets,NBTRUTH)]->Fill(T20);
    hName["T20_nJets_num"]->Fill(nJets,p0);
    hName["T20_nJets_den"]->Fill(nJets,p2);

}

double prob_1tag(int ijet, TString Atagger, TString mode){
    //copmute jet by jet probability
    //it will be one term, like e1*(1-e2) where e may or may not have a SF
    //mode=="" will compute the probability for MC
    //mode=="SF" will compute the probability with nominal SF
    //mode=="SFP/M will compute the +/- systematic event bye event

    TLorentzVector pJi;
    pJi.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));

    string TAG="";
    if(Atagger=="CSVM")TAG="M";
    if(Atagger=="CSVT")TAG="T";
    
	double eff=read_btag_efficiency(pJi,jet_algFlavor->at(ijet) ,TAG);
    double SF=1;
    if(mode=="SF") SF=getSF(ijet,Atagger,"");
    if(mode=="SFP") SF=getSF(ijet,Atagger,"max");
    if(mode=="SFM") SF=getSF(ijet,Atagger,"min");

    double probTMP=1;
    for(int k=0; k<jet_px->size(); k++){

        if(k!=ijet)probTMP*=(1-SF*eff);
        else probTMP*=SF*eff;
    }
    return probTMP;
}

double getSF(int ijet, TString Atagger, TString sys){
    //stopwatch["getSF_fillLor"].Start(kFALSE);
    TLorentzVector pJi;
    pJi.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
    //stopwatch["getSF_fillLor"].Stop();

    double x=pJi.Perp();
    int pt_bin=hName["btag_systematic_medium"]->FindBin(x);
    if(pt_bin>hName["btag_systematic_medium"]->GetNbinsX())pt_bin=hName["btag_systematic_medium"]->GetNbinsX();
 /*
    double SFb = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));// Scale factor for medium CSV b-tags
    double SFb_sys=hName["btag_systematic_medium"]->GetBinContent(pt_bin);
    double SFb_sysT=hName["btag_systematic_tight"]->GetBinContent(pt_bin);
    
    double SFb_tight = (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x)); // scale factor for tight CSV b-tags
    */
    double SF=0;
    double SFsys=0;
    int flv=TMath::Abs(jet_algFlavor->at(ijet));
    if(flv==4 || flv==5){
        if(Atagger=="CSVM") SF=(0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
        if(Atagger=="CSVT") SF=(0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));
        if(Atagger=="CSVM" && sys=="min")SFsys=SF-hName["btag_systematic_medium"]->GetBinContent(pt_bin);
        if(Atagger=="CSVM" && sys=="max")SFsys=SF+hName["btag_systematic_medium"]->GetBinContent(pt_bin);
        
        if(Atagger=="CSVT" && sys=="min")SFsys=SF-hName["btag_systematic_tight"]->GetBinContent(pt_bin);
        if(Atagger=="CSVT" && sys=="max")SFsys=SF+hName["btag_systematic_tight"]->GetBinContent(pt_bin);
        
    }
    else{
        if(Atagger=="CSVM")SF=SF_light(Atagger,"mean",pJi);
        if(Atagger=="CSVT")SF=SF_light(Atagger,"mean",pJi);
        
        if(Atagger=="CSVM" && sys=="min")SFsys=SF-SF_light(Atagger,"min",pJi);
        if(Atagger=="CSVM" && sys=="max")SFsys=SF+SF_light(Atagger,"max",pJi);
        
        if(Atagger=="CSVT" && sys=="min")SFsys=SF-SF_light(Atagger,"min",pJi);
        if(Atagger=="CSVT" && sys=="max")SFsys=SF+SF_light(Atagger,"max",pJi);
    }
    
    if(sys=="") return SF;
    if(sys=="min" || sys=="max")return SFsys;
    
}

void kinematics(){
  //  stopwatch["kinematics"].Start(kFALSE);
    
	//SFb comes from
	//https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_WITHttbar_payload_EPS13.txt
	//initialize
	W_enriched=false; 
	dimuon_enriched=false;
    dimuon_enriched_ss=false;
	Z_enriched=false;
    Z_enriched_ss=false;
	Z_enriched_loose=false;
	tt_enriched=false;
	tt_enriched_ss=false;

	tt_enriched_btag=false;
	W_enriched_bveto=false;
	Z_enriched_bveto=false;
	
	sig_1mu=false;
	sig_2mu=false;
	
	jetPsys.SetPxPyPzE(0,0,0,0);//set to 0
	jetPttsys.SetPxPyPzE(0,0,0,0);
	
    st_top=0;
    
	ht=0; 
	deltaPt=0;
	B_weight=1; 
	B_weightP=1; 
	B_weightM=1;
	
	B_weight_tight=1; 
	B_weight_tightP=1; 
	B_weight_tightM=1;
	double H[]={1,1};
	double L[]={1,1};
	double HT[]={1,1};
	double LT[]={1,1};
    probMC0=1;
    probMC1=0;
    probMC2=0;
    probMC0_SF=1;
    probMC1_SF=0;
    probMC2_SF=0;
    
    probMC0_SFP=1;
    probMC1_SFP=0;
    probMC2_SFP=0;
    
    probMC0_SFM=1;
    probMC1_SFM=0;
    probMC2_SFM=0;
	
    if(jet_px->size()>=2) {
        TLorentzVector p1;
        TLorentzVector p2;
        p1.SetPxPyPzE(jet_px->at(0),jet_py->at(0),jet_pz->at(0),jet_e->at(0));
        p2.SetPxPyPzE(jet_px->at(1),jet_py->at(1),jet_pz->at(1),jet_e->at(1));
        deltaPt=p1.Perp()-p2.Perp();
    }
    
		for(int ijet=0; ijet<jet_px->size(); ijet++){
            string label=jet_label(jet_algFlavor->at(ijet));
            if(weight>0 && weight<100)hName["jet_flavor"]->Fill(label.c_str(),weight);
         //   stopwatch["kinematics_fillLor"].Start(kFALSE);
			TLorentzVector pJi;
			pJi.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
        //    stopwatch["kinematics_fillLor"].Stop();

            // cout << "flavor: "<< jet_algFlavor->at(ijet) << " pt: "<< pJi.Perp() << " weight: " << weight << endl;
           if(weight>0 && weight<100) hName[Form("jet_pt_flavor%d",TMath::Abs(jet_algFlavor->at(ijet)) )]->Fill(pJi.Perp(),weight);
            
			jetPsys+=pJi;
			if(jet_bTagT->at(ijet)==0) jetPttsys+=pJi; //for ttbar control, only include jets that are not b-tags in the sum
            if(jet_bTagT->at(ijet)==1) st_top+=pJi.Perp();
			ht+=pJi.Perp(); 
			if(MakeEfficiency) fill_btag_efficiency( pJi, jet_algFlavor->at(ijet),jet_bTagM->at(ijet),jet_bTagT->at(ijet)); 
			double x=pJi.Perp();
			int pt_bin=hName["btag_systematic_medium"]->FindBin(x);
			if(pt_bin>hName["btag_systematic_medium"]->GetNbinsX())pt_bin=hName["btag_systematic_medium"]->GetNbinsX();
			//stopwatch["kinematics_SF"].Start(kFALSE);
            /*
            double SFb = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));// Scale factor for medium CSV b-tags
			double SFb_sys=hName["btag_systematic_medium"]->GetBinContent(pt_bin);
			double SFb_sysT=hName["btag_systematic_tight"]->GetBinContent(pt_bin);

			double SFb_tight = (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x)); // scale factor for tight CSV b-tags
			double SF_lM=SF_light("CSVM","mean",pJi); 
			double SF_lT=SF_light("CSVT","mean",pJi); 
			
			double SF_lM_sysmin=SF_light("CSVM","min",pJi); 
			double SF_lM_sysmax=SF_light("CSVM","max",pJi); 

			double SF_lT_sysmin=SF_light("CSVT","min",pJi);
			double SF_lT_sysmax=SF_light("CSVT","max",pJi);
*/
        //    stopwatch["kinematics_SF"].Stop();

			
			//cout << "btag pt: " <<pJi.Perp() << " flavor: "<< jet_algFlavor->at(ijet) <<  " tag: "<< jet_bTagM->at(ijet) <<  " effM: " << read_btag_efficiency(pJi,jet_algFlavor->at(ijet) ,"M") << endl;
			//cout << "SF: " << SFb << " SFb+: " << SFb_sys << " SFl: " << SF_lM << " " << SF_lM_sysmin << " " << SF_lM_sysmax << endl;
			if(_MC){
             //   stopwatch["kinematics_SF"].Start(kFALSE);

				double SF=getSF(ijet,"CSVM","");
				double SF_tight=getSF(ijet,"CSVT","");
				double SFsysP=getSF(ijet,"CSVM","max");
				double SFsysM=getSF(ijet,"CSVM","min");

				double SFsysP_tight=getSF(ijet,"CSVT","max");
				double SFsysM_tight=getSF(ijet,"CSVT","min");
          //      stopwatch["kinematics_SF"].Stop();

                
				bool heavy=false;
				if(TMath::Abs(jet_algFlavor->at(ijet))==4 || TMath::Abs(jet_algFlavor->at(ijet))==5) {
					heavy=true;
				}
             //   stopwatch["kinematics_eff"].Start(kFALSE);

				double eff=read_btag_efficiency(pJi,jet_algFlavor->at(ijet) ,"M"); 
				double eff_tight=read_btag_efficiency(pJi,jet_algFlavor->at(ijet) ,"T"); 
             //   stopwatch["kinematics_eff"].Stop();

                probMC0*=1-eff_tight;
                probMC0_SF*=1-SF_tight*eff_tight;
                probMC0_SFP*=1-SFsysP_tight*eff_tight;

               // probMC1+=prob_1tag(ijet,"CSVT",""); //compute nominal probability for each term
               // probMC1_SF+=prob_1tag(ijet,"CSVT","SF"); //compute data scaled porbability for each term.
                //probMC1_SFP+=prob_1tag(ijet,"CSVT","SFP");
                
				//cout << " tag: "<< jet_bTagM->at(ijet) << " pt: " << pJi.Perp() << " eff: "<< eff << " SF: " << SF; 
				//cout << "SF+: " << SFsysP << " " << SFsysM << endl; 
				if(jet_bTagM->at(ijet)==1) {
					B_weight*=SF*eff/eff;
					if(heavy){
						H[0]*=SFsysP;
						H[1]*=SFsysM;
						L[0]*=SF;
						L[1]*=SF;
					}
					else {
						H[0]*=SF;
						H[1]*=SF;
						L[0]*=SFsysP;
						L[1]*=SFsysM;
					}

				}
				else {
					B_weight*=(1-SF*eff)/(1-eff); 
					if(heavy){
						H[0]*=(1-SFsysP*eff)/(1-eff);
						H[1]*=(1-SFsysM*eff)/(1-eff);
						L[0]*=(1-SF*eff)/(1-eff);
						L[1]*=(1-SF*eff)/(1-eff);
					}
					else {
						H[0]*=(1-SF*eff)/(1-eff);
						H[1]*=(1-SF*eff)/(1-eff);
						L[0]*=(1-SFsysP*eff)/(1-eff);
						L[1]*=(1-SFsysM*eff)/(1-eff);
					}

				}
				
				if(jet_bTagT->at(ijet)==1){
                    
					B_weight_tight*=SF_tight*eff_tight/eff_tight; 
					if(heavy){
						HT[0]*=SFsysP_tight;
						HT[1]*=SFsysM_tight;
						LT[0]*=SF_tight;
						LT[1]*=SF_tight;
					}
					else {
						HT[0]*=SF_tight;
						HT[1]*=SF_tight;
						LT[0]*=SFsysP_tight;
						LT[1]*=SFsysM_tight;
					}
					
				}
				else {
					B_weight_tight*=(1-SF_tight*eff_tight)/(1-eff_tight);
					if(heavy){
						H[0]*=(1-SFsysP_tight*eff_tight)/(1-eff_tight);
						H[1]*=(1-SFsysM*eff_tight)/(1-eff_tight);
						L[0]*=(1-SF_tight*eff_tight)/(1-eff_tight);
						L[1]*=(1-SF_tight*eff_tight)/(1-eff_tight);
					}//heavy flavor
					else {
						H[0]*=(1-SF_tight*eff_tight)/(1-eff_tight);
						H[1]*=(1-SF_tight*eff_tight)/(1-eff_tight);
						L[0]*=(1-SFsysP_tight*eff_tight)/(1-eff_tight);
						L[1]*=(1-SFsysM_tight*eff_tight)/(1-eff_tight);
					}//light flavor 

				}
			}
			
		}
	probMC2=TMath::Abs(1-probMC0-probMC1); // because of round off, you can get numbers like -1 * 10^-17
    probMC2_SF=TMath::Abs(1-probMC0_SF-probMC1_SF);
    probMC2_SFP=TMath::Abs(1-probMC0_SFP-probMC1_SFP);
    
    probMC0_2=probMC0/probMC2;
    probMC0_2_SF=probMC0_SF/probMC2_SF;
    probMC0_2_SFP=probMC0_SFP/probMC2_SFP;

    //cout << "sum: " << probMC0+probMC1+probMC2 << endl;
    //cout << "sum SF: " << probMC0_SF+probMC1_SF+probMC2_SF << endl;
    //cout << "sum SF+: " << probMC0_SFP+probMC1_SFP+probMC2_SFP << endl;

    //cout << "nJets: " << nJets << " probMC(0) " << probMC0 << " probMC(1) " << probMC1 << " probMC(2) " << probMC2 << endl;
    //cout << "SF nJets: " << nJets << " probMC(0) " << probMC0_SF << " probMC(1) " << probMC1_SF << " probMC(2) " << probMC2_SF << endl;
    //cout << "SF+ nJets: " << nJets << " probMC(0) " << probMC0_SFP << " " << probMC1_SFP << " " << probMC2_SFP << endl;
    
    //if(nJets>=2)cout << "P(0/2): " << probMC0_2 << " " << probMC0_2_SF << " " <<probMC0_2_SFP << endl;
   // if(nJets>=2)    transfer_factor("CSVT");

	for(int i=0; i<=2; i++) {
		if(i==0)B_weightP=B_weight+TMath::Sqrt(TMath::Power(H[i]-B_weight,2)+TMath::Power(L[i]-B_weight,2)); 
		if(i==1)B_weightM=B_weight-TMath::Sqrt(TMath::Power(H[i]-B_weight,2)+TMath::Power(L[i]-B_weight,2));
		if(i==0)B_weight_tightP=B_weight_tight+TMath::Sqrt(TMath::Power(HT[i]-B_weight_tight,2)+TMath::Power(LT[i]-B_weight_tight,2)); 
		if(i==1)B_weight_tightM=B_weight_tight-TMath::Sqrt(TMath::Power(HT[i]-B_weight_tight,2)+TMath::Power(LT[i]-B_weight_tight,2));
	}
	
    B_weight_tightPbc=HT[0];
    B_weight_tightMbc=HT[1];

    B_weight_tightPl=LT[0];
    B_weight_tightMl=LT[1];
    
    B_weightPbc=H[0];
    B_weightMbc=H[1];
    
    B_weightPl=L[0];
    B_weightMl=L[1];
    
	//cout << "B_weight: " << B_weight << " + " << B_weight_sys[0] << " - " << B_weight_sys[1] << endl;
	//cout << "B_weight_tight: " << B_weight_tight << " + " << B_weight_sys_tight[0] << " - " << B_weight_sys_tight[1] << endl;

		//cout << "N b-tag: " << weighti.size() << endl; 
		calcMt();
	dimuon_mass=-1; 
	pTmuon1=0; 
	pTmuon2=0;
	dilepton_pt=0;
	dilepton_mass=0;
	
	TLorentzVector pM;
	if(muon_px->size()>=1){
		pM.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0)); 
		pTmuon1=pM.Perp(); 
	}
	 
	if(muon_px->size()==2){
		TLorentzVector pM1;
		TLorentzVector pM2; 
		pM1.SetPxPyPzE(muon_px->at(0),muon_py->at(0),muon_pz->at(0),muon_e->at(0)); 
		pM2.SetPxPyPzE(muon_px->at(1),muon_py->at(1),muon_pz->at(1),muon_e->at(1)); 
		pM1+=pM2;
		dimuon_mass=pM1.M();
	}
	st=calcSt();
	control_regions();
    st_top+=pTmuon1;
    st_top+=pTmuon2;
    st_top+=met;

//stopwatch["kinematics"].Stop();

}

void fill_jets(){
	bool print=false; 
	if(print) cout << "fill jets" << endl; 
	hName2D[Form("st_ptsys_nJets%d",nJets)]->Fill(st,jetPsys.Perp(),weight);
	hName2D[Form("st_met_nJets%d",nJets)]->Fill(st,met,weight);

	
	hName2D["jetpt_nJets"]->Fill(nJets,jetPsys.Perp(),weight); 
	if(W_enriched){
		hName[Form("jetpTsys%d_W",nJets)]->Fill(jetPsys.Perp(),weight); 
		for(int ijet=0; ijet<jet_px->size(); ijet++) {
			if(ijet>=nJetmax) continue;
			TLorentzVector pJ;
			pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet), jet_pz->at(ijet), jet_e->at(ijet)); 
			//cout << "pt jet: " << pJ.Perp() << endl;
            hName[Form("W_jetpt_jet%d_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight);
			if(st>300 && st<=600) hName[Form("W_jetpt_jet%d_st1_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			if(st>600 && st<=1000) hName[Form("W_jetpt_jet%d_st2_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			if(st>1000) hName[Form("W_jetpt_jet%d_st3_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 

		}
	}
	
	if(Z_enriched){
		hName[Form("jetpTsys%d_Z",nJets)]->Fill(jetPsys.Perp(),weight); 
		
		for(int ijet=0; ijet<jet_px->size(); ijet++) {
			if(ijet>=nJetmax) continue;
			TLorentzVector pJ;
			pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet), jet_pz->at(ijet), jet_e->at(ijet));
            hName[Form("Z_jetpt_jet%d_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight);
			if(st>300 && st<=600) hName[Form("Z_jetpt_jet%d_st1_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			if(st>600 && st<=1000) hName[Form("Z_jetpt_jet%d_st2_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			if(st>1000) hName[Form("Z_jetpt_jet%d_st3_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			
		}
	}
	
	
	if(tt_enriched){
        int NBT=nTightBJets;
        if(NBT>2)NBT=2;
        int Nl=0;
        int Ng=0;
        int Nb=0;
        int Nc=0;
		if(NBT>2)hName[Form("jetpTsys%d_1Mu_1El",nJets)]->Fill(jetPttsys.Perp(),weight);
       // cout << "nJets: "<< nJets << " nBT: " << NBT << endl;
       // cout << "jet_px size: "<< jet_px->size() << endl;
		for(int ijet=0; ijet<jet_px->size(); ijet++) {
			if(ijet>=nJetmax) continue;
            if(jet_algFlavor->at(ijet)<=3) Nl++;
            if(jet_algFlavor->at(ijet)==21)Ng++;
            if(jet_algFlavor->at(ijet)==4)Nc++;
            if(jet_algFlavor->at(ijet)==21)Nb++;
            
            TLorentzVector pJ;
         //   cout << "ij: "<< ijet << endl;
			pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet), jet_pz->at(ijet), jet_e->at(ijet));
           // cout <<"ijet: " << ijet << " nJets " << nJets << " nB " << NBT << endl;
            hName[Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight",ijet,nJets,NBT)]->Fill(pJ.Perp(),weight);
            hName[Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight_SF",ijet,nJets,NBT)]->Fill(pJ.Perp(),weight*B_weight_tight);
            hName[Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight_SFP",ijet,nJets,NBT)]->Fill(pJ.Perp(),weight*B_weight_tightP);
            hName[Form("tt_jetpt_jet%d_nJets%d_%dbtags_tight_SFM",ijet,nJets,NBT)]->Fill(pJ.Perp(),weight*B_weight_tightM);

            hName[Form("tt_jetpt_jet%d_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight);
			if(st>300 && st<=600) hName[Form("tt_jetpt_jet%d_st1_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			if(st>600 && st<=1000) hName[Form("tt_jetpt_jet%d_st2_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
			if(st>1000) hName[Form("tt_jetpt_jet%d_st3_nJets%d",ijet,nJets)]->Fill(pJ.Perp(),weight); 
		}
        hName[Form("nb_nJets%d_%dbtag",nJets,NBT)]->Fill(Nb);
        hName[Form("nc_nJets%d_%dbtag",nJets,NBT)]->Fill(Nc);
        hName[Form("nl_nJets%d_%dbtag",nJets,NBT)]->Fill(Nl);
        hName[Form("ng_nJets%d_%dbtag",nJets,NBT)]->Fill(Ng);

        prName[Form("jet_flavor_nJets%d",nJets)]->Fill("n_{b}",Nb);
        prName[Form("jet_flavor_nJets%d",nJets)]->Fill("n_{c}",Nc);
        prName[Form("jet_flavor_nJets%d",nJets)]->Fill("n_{g}",Ng);
        prName[Form("jet_flavor_nJets%d",nJets)]->Fill("n_{l}",Nl);

        prName[Form("jet_flavor_nJets%d_%dbtag",nJets,NBT)]->Fill("n_{b}",Nb);
        prName[Form("jet_flavor_nJets%d_%dbtag",nJets,NBT)]->Fill("n_{c}",Nc);
        prName[Form("jet_flavor_nJets%d_%dbtag",nJets,NBT)]->Fill("n_{g}",Ng);
        prName[Form("jet_flavor_nJets%d_%dbtag",nJets,NBT)]->Fill("n_{l}",Nl);

	}
	TLorentzVector pJ;
	//if(nJets<1 || nJets>3) return; 
	for(int ijet=0; ijet<jet_px->size(); ijet++){
		int jet=ijet+1;
		if(jet>nJetmax)jet=nJetmax; 
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet)); 
		hName[Form("jetpT%d_%dMu",jet,muon_px->size())]->Fill(pJ.Perp(),weight); 
		
		if(tt_enriched){
			hName[Form("jetpT%d_1Mu_1El",jet)]->Fill(pJ.Perp(),weight);
            if(jet_bTagT->at(ijet)) hName["jet_pt_1Mu_1El_btag_tight_SF"]->Fill(pJ.Perp(),weight*B_weight );
			if(jet_bTagM->at(ijet))hName[Form("bjetpT%d_1Mu_1El",jet)]->Fill(pJ.Perp(),weight*B_weight); 
			if(jet_bTagT->at(ijet))hName[Form("bjetpT%d_1Mu_1El_tight",jet)]->Fill(pJ.Perp(),weight*B_weight_tight);
			if(jet_bTagT->at(ijet)==0)hName[Form("0bjetpT%d_1Mu_1El_tight",jet)]->Fill(pJ.Perp(),weight*B_weight_tight);

		}
		
		if(W_enriched){
			hName[Form("jetpT%d_W",jet)]->Fill(pJ.Perp(),weight); 
		}
		if(Z_enriched){
			hName[Form("jetpT%d_Z",jet)]->Fill(pJ.Perp(),weight); 
			if(jet_bTagM->at(ijet))hName[Form("bjetpT%d_Z",jet)]->Fill(pJ.Perp(),weight*B_weight_tight); 
			if(jet_bTagT->at(ijet))hName[Form("bjetpT%d_Z_tight",jet)]->Fill(pJ.Perp(),weight*B_weight_tight); 

		}
		
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


double read_btag_efficiency(TLorentzVector pJ, int flavor,string Tagger){
	int ieta=0; 
	flavor=TMath::Abs(flavor);
	if(flavor<4) flavor=3; //udgs light flavor
	//if(flavor==21) flavor=3; //udgs light flavor
	double etaBin[]={0,0.8,1.6,2.4};
	
	double eta=TMath::Abs(pJ.Eta()); 
	
	for(int i=0; i<3; i++){
		if(eta>etaBin[i] && eta<=etaBin[i+1]) ieta=i; 
	}
	
	TH1F *h=(TH1F*)effName[Form("btag%s_eff_eta%d_flavor%d",Tagger.c_str(),ieta,flavor)]->GetPassedHistogram();
	int ptbin=h->FindBin(pJ.Perp()); 
	if(ptbin > h->GetNbinsX()) ptbin=h->GetNbinsX(); 
	
	double eff=effName[Form("btag%s_eff_eta%d_flavor%d",Tagger.c_str(),ieta,flavor)]->GetEfficiency(ptbin); 
	return eff;
	
}

void fill_btag_efficiency(TLorentzVector pJ, int flavor, bool btagM, bool btagT){
	int ieta=0; 
	flavor=TMath::Abs(flavor);
	if(flavor<4) flavor=3; //udgs light flavor
	//if(flavor==21) flavor=3; //udgs light flavor
	double etaBin[]={0,0.8,1.6,2.4};
	
	double eta=TMath::Abs(pJ.Eta()); 
	
	for(int i=0; i<3; i++){
		if(eta>etaBin[i] && eta<=etaBin[i+1]) ieta=i; 
	}
	
	effName[Form("btagM_eff_eta%d_flavor%d",ieta,flavor)]->Fill(btagM,pJ.Perp());
	effName[Form("btagT_eff_eta%d_flavor%d",ieta,flavor)]->Fill(btagT,pJ.Perp());
	
	if(st>300){
		effName[Form("btagM_eff_eta%d_flavor%d_stmin",ieta,flavor)]->Fill(btagM,pJ.Perp());
		effName[Form("btagT_eff_eta%d_flavor%d_stmin",ieta,flavor)]->Fill(btagT,pJ.Perp());
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
/*
void check_jet_cleaning(){
    std::vector<int> jetOverlap_mu;
    std::vector<int> jetOverlap_e;
    
    for(int ij=0; ij<jet_px->size();ij++){
        TLorentzVector pJ;
        pJ.SetPxPyPzE(jet_px->at(ij),jet_py->at(ij),jet_pz->at(ij),jet_e->at(ij));
        for(int imu=0; imu<muon_px->size();imu++){
            TLorentzVector pMu;
            pMu.SetPxPyPzE(muon_px->at(imu),muon_py->at(imu),muon_pz->at(imu),muon_e->at(imu));
            double dR=pJ.DeltaR(pMu);
            if(dR<0.5) jetOverlap_mu.push_back(ij);
            
        }
        for(int imu=0; imu<muon_px->size();imu++){
            //yes... the loop variabe is mu. But notice that the 4 vector is filled by the electron variable
            TLorentzVector pMu;
            pMu.SetPxPyPzE(electron_px->at(imu),electron_py->at(imu),electron_pz->at(imu),electron_e->at(imu));
            double dR=pJ.DeltaR(pMu);
            if(dR<0.5) jetOverlap_e.push_back(ij);
        }
        
    }
    
    cout << "bad jets: " << jetOverlap_mu.size() << " " << jetOverlap_e.size() << endl;
}
*/
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
	/*
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
	
	std::vector<float>::size_type ijet = 0;
	
	while ( ijet < jet_px->size() ) {
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
		if(TMath::Abs(pJ.DeltaPhi(pMet))<DeltaPhiCut_jetMet) remove_jet(ijet);		
		else {
			++ijet;
		}
	}
	 
	jetCut_pt1=75; // remove first jet if less than 75 GeV
	if(nJets>=1) {
		pJ.SetPxPyPzE(jet_px->at(0),jet_py->at(0),jet_pz->at(0),jet_e->at(0)); 
		
		if(pJ.Perp()<jetCut_pt1){
			remove_jet(0);
		}
	}
	
	if(print)cout << "after delta Phi cut: " << endl; 
	if(print)print_jets();
	
	nJets=jet_px->size();
	
	int Nmu=muon_px->size(); 

	if(Nmu<1 || Nmu>2) return false; 
	if(nJets<1) return false;
	st=calcSt();
	if(st<300) return false; 
	*/
    
    std::vector<float>::size_type ijet = 0;
	
	while ( ijet < jet_px->size() ) {
		pJ.SetPxPyPzE(jet_px->at(ijet),jet_py->at(ijet),jet_pz->at(ijet),jet_e->at(ijet));
		if(pJ.Perp()<100) remove_jet(ijet);
		else {
			++ijet;
		}
	}
    
	
}

bool basic_selections(){
	if(st<300) return false; 
	//if(mt<85) return false; 
	if(muon_px->size()<1 || muon_px->size()>2) return false;
	if(muon_px->size()==2) {
		//if(muon_charge->at(0)*muon_charge->at(1)==1) return false;
		if(dimuon_mass<50) return false; 
	}
	
	if(W_enriched==0 && Z_enriched==0 && dimuon_enriched==0 && dimuon_enriched_ss==0 && tt_enriched==0 && tt_enriched_ss==0 && sig_1mu==0 && sig_2mu==0) return false;

	return true; 
}

void pileup_reweight(){
	if(nJets>=2 && nJets<=4 && muon_px->size()>=1) hName["h_before_npv"]->Fill(nPv,weight);
	//cout << "weight: "<< weight << " puWt: " << puWt_nom << endl; 
	if(_MC) weight=weight*puWt_nom;	
	//cout << "weight: " << weight<< endl; 
	
}

void event_loop(TChain *tree){
	int maxEvt=100000;
	//entries=maxEvt;
	
	TStopwatch t; 
	t.Start(kFALSE); 
	
    TStopwatch total;
	total.Start(kFALSE);
    
    TStopwatch fill_st_SW;
    TStopwatch kinematics_SW;
    TStopwatch btag_SW;
    TStopwatch gen_pt_SW;
    
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

		count_muons();
		count_electrons();
		weight=W; 
		//cout << "weight before eff: " << weight << endl; 
	
		efficiency(); 
		pileup_reweight();
        btag_SW.Start(kFALSE);
        nBtag();
        btag_SW.Stop();

		//fill_deltaR(); 
		//fill_deltaPhi();
		//fill_rapidity();
        //check_jet_cleaning();
		//if(!cuts()) continue;
		if(nJets>nJetmax)nJets=nJetmax;
        gen_pt_SW.Start(kFALSE);
		gen_pt();
        gen_pt_SW.Stop();
        
        kinematics_SW.Start(kFALSE);
		kinematics();
        kinematics_SW.Stop();
        
        //fill_bjet_mu_mass();
		if(basic_selections()==false) continue; 
		
		//fill_jets();
		fill_muon(); 
		fill_pileup();

        fill_st_SW.Start(kFALSE);
		fill_st();
        fill_st_SW.Stop();
        
        fill_event_category();
		
		/*
		kinematics();
		basic_selections(); 
		nBtag();
		acceptance(); 
		 */
        //fill_categories();
		
	}
    cout << "acceptance: " << Nacc << endl; 
    
    double Tbc=hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP_bc",0)]->GetBinContent(2)/hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP_bc",2)]->GetBinContent(2);
    double Tl=hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP_l",0)]->GetBinContent(2)/hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP_l",2)]->GetBinContent(2);
    double TE_=hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP",0)]->GetBinContent(2)/hName[Form("nJets_1Mu_1El_%dbtag_tight_SFP",2)]->GetBinContent(2);
    double T=hName[Form("nJets_1Mu_1El_%dbtag_tight_SF",0)]->GetBinContent(2)/hName[Form("nJets_1Mu_1El_%dbtag_tight_SF",2)]->GetBinContent(2);

    double TE=TMath::Sqrt(TMath::Power(Tbc-T,2)+TMath::Power(Tl-T,2));
    
   // cout << "SumW: " << sumW_EMu << "  Nevents: " << sum_EMu << endl;
    
    //cout << "T : " << T << " TE: " << TE << " My method: "<< (TE_-T) << endl;
    
    total.Stop();
//	cout << "fill st: " << fill_st_SW.CpuTime() << endl;
//	cout << "fill kin: " << kinematics_SW.CpuTime() << endl;
 /*   cout << "kinematics: " << stopwatch["kinematics"].CpuTime() << endl;
    cout << "kinematics lorentz vectors: " << stopwatch["kinematics_fillLor"].CpuTime() << endl;
    cout << "kinematics SF: " << stopwatch["kinematics_SF"].CpuTime() << endl;
    cout << "kinematics eff: " << stopwatch["kinematics_eff"].CpuTime() << endl;
    cout << "getSF fill Lor: " << stopwatch["getSF_fillLor"].CpuTime()<<endl;
    cout << "SF_light: " << stopwatch["SF_light"].CpuTime() << endl;*/

//   cout << "gen_pt: " << gen_pt_SW.CpuTime() << endl;
 //   cout << "btag: " <<btag_SW.CpuTime() << endl;
 //   cout << "total time: " << total.CpuTime() << endl;
}

void acceptance(){
	if(tt_enriched && nBJets==0 && nJets>=4){
		Nacc+=1./NG; 
	}
}

void analyze_file(TString folder){
	load_btag_sys();
	if(!MakeEfficiency) load_btagEff();
	get_muTrigEff();
	get_corrections();
	//TTree *tree = (TTree*)infile->Get("tree");
	TChain *tree = new TChain("tree"); 
	if(folder.Contains("w"))Correct_=true; 
	if(folder.Contains("dy"))Correct_=true; 
	TString FN=folder+"/hist_analysis*.root";
	cout << "FN: " << FN << endl; 
	int Nfiles=tree->Add(FN);
	int NGEN=0;
	for(int ii=0; ii<Nfiles;ii++){
		TString FN_i=folder;
		FN_i=FN_i+Form("/hist_analysis_%d.root",ii); 
		cout << "FN_i: "<< FN_i << endl; 
		TFile *f=new TFile(FN_i,"READ"); 
		cout << "Get number of events: " ; 
		TH1F *h=(TH1F*)f->Get("h_Nevents");
		NGEN+=h->GetBinContent(1); 
		//cout << "N: " << h->GetBinContent(1) << endl; 
		delete h;
		delete f; 
	}
	NG=NGEN; 
	cout << "Actual number processed: " << NGEN << endl; 
	
	if(_MC) weight=weight/NGEN;
	
	cout << "actual weight used: "<< weight << endl; 
	
	cout << "Intitialize Tree: " << endl; 
	initialize_tree(tree);
	expo_corr=new TF1("expo_cor", "expo",0,1000);
	expo_corr->SetParameters(0.156,-0.00137);
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
	
	xs["wJetsMatchingDown"]=36632; 
	Ngen["wJetsMatchingDown"]=10000; 
	fileName["wJetsMatchingDown"]="wJets_matchingdown";

	xs["wJetsht150"]=290.7;
	xs["wJetsht200"]=111.37;
	xs["wJetsht250"]=59.24;
	xs["wJetsht300"]=47.25;
	xs["wJetsht400"]=31.12;
	
	
	if(LO){
		xs["wJetsht150"]=235.6;
		xs["wJetsht200"]=90.27;
		xs["wJetsht250"]=48.01;
		xs["wJetsht300"]=38.3;
		xs["wJetsht400"]=25.22;
	}
	
	Ngen["wJetsht150"]=21686209;
	Ngen["wJetsht200"]=10039771;
	Ngen["wJetsht250"]=6575572;
	Ngen["wJetsht300"]=6840509;
	Ngen["wJetsht400"]=6619654;
	
	fileName["wJetsht150"]="wJets_ht_150To200";
	fileName["wJetsht200"]="wJets_ht_200To250";
	fileName["wJetsht250"]="wJets_ht_250To300";
	fileName["wJetsht300"]="wJets_ht_300To400";
	fileName["wJetsht400"]="wJets_ht_400Toinf";
	
	xs["dy1Jets"]=666;
	xs["dy2Jets"]=215;
	xs["dy3Jets"]=61;
	xs["dy4Jets"]=27;
	
	Ngen["dy1Jets"]=24045248-1200000;
	Ngen["dy2Jets"]=21852156;
	Ngen["dy3Jets"]=11015445;
	Ngen["dy4Jets"]=6402827;
	
	xs["dyJets_matchingdown"]=3403.08;
	xs["dyJets_matchingup"]=3462.15;
	xs["dyJets_scaledown"]=3933.66;
	xs["dyJets_scaleup"]=3508.47;

	Ngen["dyJets_matchingdown"]=1;
	Ngen["dyJets_matchingup"]=1;
	Ngen["dyJets_scaledown"]=1; 
	Ngen["dyJets_scaleup"]=1;
	
	fileName["dyJets_matchingdown"]="dyJets_matchingdown";
	fileName["dyJets_matchingup"]="dyJets_matchingup";
	fileName["dyJets_scaledown"]="dyJets_scaledown"; 
	fileName["dyJets_scaleup"]="dyJets_scaleup";
	
	xs["ttSemiLept"]=107.6;
	xs["ttFullLept"]=25.6;
	
	fileName["ttSemiLept"]="ttJetsSemiLept";
	fileName["ttFullLept"]="ttJetsFullLept";

	Ngen["ttSemiLept"]=86711159-840000;
	Ngen["ttFullLept"]=12119013-630000;
	
	
	xs["ttJets_matchingdown"]=243.44;
	xs["ttJets_matchingup"]=248.86;
	xs["ttJets_scaledown"]=411.16;
	xs["ttJets_scaleup"]=176.19;

	Ngen["ttJets_matchingdown"]=1;
	Ngen["ttJets_matchingup"]=1;
	Ngen["ttJets_scaledown"]=1;
	Ngen["ttJets_scaleup"]=1;
	
	fileName["ttJets_matchingdown"]="ttJets_matchingdown";
	fileName["ttJets_matchingup"]="ttJets_matchingup";
	fileName["ttJets_scaledown"]="ttJets_scaledown";
	fileName["ttJets_scaleup"]="ttJets_scaleup";
	

	xs["ttG"]=25.24;
	Ngen["ttG"]=1719954;
	fileName["ttG"]="ttGJets";
	
	
	//single top xs
	
	fileName["TBar_t"]="TBar_t";
	fileName["TBar_s"]="TBar_s";
	fileName["TBar_tW"]="TBar_tW";
	
	fileName["T_t"]="T_t";
	fileName["T_s"]="T_s";
	fileName["T_tW"]="T_tW"; 
	
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
	
	fileName["UDD300"]="UDD_300";
	xs["UDD300"]=1.99608; 
	Ngen["UDD300"]=98000;
	
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
	
	fileName["stealth_300_200"]="stealth_300_200"; 
	xs["stealth_300_200"]=0.543*19828.3/1000;
	Ngen["stealth_300_200"]=238886;
	
	fileName["stealth_400_200"]="stealth_400_200"; 
	xs["stealth_400_200"]=0.543*3.54338;
	Ngen["stealth_400_200"]=238886;
	
	fileName["stealth_500_300"]="stealth_500_300"; 
	xs["stealth_500_300"]=0.543*847.051/1000; 
	Ngen["stealth_500_300"]=239562;
	
	fileName["stealth_600_300"]="stealth_600_300"; 
	xs["stealth_600_300"]=0.543*0.244862; 
	Ngen["stealth_600_300"]=239562;
	
	fileName["stealth_700_400"]="stealth_700_400"; 
	xs["stealth_700_400"]=0.543*0.0799667; 
	Ngen["stealth_700_300"]=239562;
	
	fileName["stealth_500_400"]="stealth_500_400"; 
	xs["stealth_500_400"]=0.543*847.051/1000; 
	Ngen["stealth_500_400"]=239562; 
	
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
	
	xs["test"]=1.0;
	Ngen["test"]=1.0; 
	fileName["test"]="test_tt"; //change this line for the file name when running test job. 
	
}

void setup_files(int jobnumber){

	if(jobnumber==-1){
	//test job, not standard 
		sample_list.push_back("stealth_300_200"); 
		test=true;
		_MC=true;
		output_file_name="stealth_300_200_LeptonCut_15_15";
		selections="_ptmu15_pte15";
		
	}
	
	if(jobnumber==0) {
		sample_list.push_back("singleMu"); 
		_MC=false; 
		output_file_name="singleMu"; 
		selections="_trigger";
	}

	if(jobnumber==1) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("dy%dJets",nJ);
			sample_list.push_back(tag); 
		}
		_MC=true; 
		output_file_name="dy"; 
		selections="_trigger";

	}
	//dy systematics
	if(jobnumber==2) {
		sample_list.push_back("dyJets_scaleup"); 
		
		_MC=true; 
		output_file_name="dyJets_scaleup"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==3) {
		sample_list.push_back("dyJets_scaledown"); 
		
		_MC=true; 
		output_file_name="dyJets_scaledown"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==4) {
		sample_list.push_back("dyJets_matchingup"); 
		
		_MC=true; 
		output_file_name="dyJets_matchingup"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==5) {
		sample_list.push_back("dyJets_matchingdown"); 
		
		_MC=true; 
		output_file_name="dyJets_matchingdown"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==6) {
		sample_list.push_back("dyJets_scaleup"); 
		
		_MC=true; 
		output_file_name="dyJets_scaleup"; 
		selections="_trigger";
		
	}
	//ttbar semi+fully leptonic
	if(jobnumber==7) {
		sample_list.push_back("ttSemiLept"); 
		sample_list.push_back("ttFullLept"); 
		//sample_list.push_back("ttG"); 
		
		_MC=true; 
		output_file_name="ttbar"; 
		selections="_trigger";

	}
	
	if(jobnumber==8) {
		sample_list.push_back("ttSemiLept"); 
			
		_MC=true; 
		output_file_name="ttSemiLept"; 
		selections="_trigger";

	}
	if(jobnumber==9) {
		sample_list.push_back("ttFullLept"); 
		
		_MC=true; 
		output_file_name="ttFullLept"; 
		selections="_trigger";

	}
	if(jobnumber==10) {
		sample_list.push_back("ttJets_scaleup"); 
		
		_MC=true; 
		test=true;
		output_file_name="ttJets_scaleup"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==11) {
		sample_list.push_back("ttJets_scaledown"); 
		
		_MC=true; 
		output_file_name="ttJets_scaledown"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==12) {
		sample_list.push_back("ttJets_matchingup"); 
		
		_MC=true; 
		output_file_name="ttJets_matchingup"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==13) {
		sample_list.push_back("ttJets_matchingdown"); 
		
		_MC=true; 
		output_file_name="ttJets_matchingdown"; 
		selections="_trigger";
		
	}
	
	
	if(jobnumber==14) {
		sample_list.push_back("QCD_30-50"); 
		sample_list.push_back("QCD_50-80");
		sample_list.push_back("QCD_80-120");
		sample_list.push_back("QCD_120-170");
		sample_list.push_back("QCD_170-300");
		sample_list.push_back("QCD_300-470");
		sample_list.push_back("QCD_470-600");
		sample_list.push_back("QCD_600-800");
		sample_list.push_back("QCD_800-1000");

		_MC=true; 
		output_file_name="QCD"; 
		selections="_trigger";
		
	}
		

	if(jobnumber==15){
		sample_list.push_back("TBar_t");
		sample_list.push_back("TBar_s");
		sample_list.push_back("TBar_tW");
		
		sample_list.push_back("T_t");
		sample_list.push_back("T_s");
		sample_list.push_back("T_tW");
		_MC=true; 
		output_file_name = "singleTop"; 
		selections="_trigger"; 

	}
	
	if(jobnumber==16){	
		sample_list.push_back("WWJetsTo2L2Nu");
		sample_list.push_back("WZJetsTo2L2Q");
		sample_list.push_back("WZJetsTo3LNu");
		sample_list.push_back("WZJetsTo2QLNu");
		
		sample_list.push_back("ZZJetsTo2L2Q");
		sample_list.push_back("ZZJetsTo2L2Nu");
		sample_list.push_back("ZZJetsTo4L");
		
		_MC=true; 
		QCD=false;
		output_file_name = "diboson"; 
		selections="_trigger"; 
		
	}
	//stealth
	if(jobnumber==17) {
		sample_list.push_back("stealth_300_200"); 
		
		_MC=true; 
		output_file_name="stealth_300_200"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==18) {
		sample_list.push_back("stealth_400_200"); 
		
		_MC=true; 
		output_file_name="stealth_400_200"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==19) {
		sample_list.push_back("stealth_500_300"); 
		
		_MC=true; 
		output_file_name="stealth_500_300"; 
		selections="_trigger";
		
	}
	
	
	if(jobnumber==20) {
		sample_list.push_back("stealth_600_300"); 
		
		_MC=true; 
		output_file_name="stealth_600_300"; 
		selections="_trigger";
		
	}
	if(jobnumber==21) {
		sample_list.push_back("stealth_700_400"); 
		
		_MC=true; 
		output_file_name="stealth_700_400"; 
		selections="_trigger";
		
	}
	//RPV
	if(jobnumber==22) {
		sample_list.push_back("RPV300"); 
		
		_MC=true; 
		output_file_name="RPV300"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==23) {
		sample_list.push_back("RPV400"); 
		
		_MC=true; 
		output_file_name="RPV400"; 
		selections="_trigger";
		
	}
	
	
	if(jobnumber==24) {
		sample_list.push_back("RPV500"); 
		
		_MC=true; 
		output_file_name="RPV500"; 
		selections="_trigger";
		
	}
	if(jobnumber==25) {
		sample_list.push_back("RPV600"); 
		
		_MC=true; 
		output_file_name="RPV600"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==26) {
		sample_list.push_back("RPV700"); 
		
		_MC=true; 
		output_file_name="RPV700"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==27) {
		sample_list.push_back("RPV800"); 
		
		_MC=true; 
		output_file_name="RPV800"; 
		selections="_trigger";
		
	}
	
	
	if(jobnumber==28) {
		sample_list.push_back("RPV900"); 
		
		_MC=true; 
		output_file_name="RPV900"; 
		selections="_trigger";
		
	}
	
	if(jobnumber==29) {
		sample_list.push_back("RPV1000"); 
		
		_MC=true; 
		output_file_name="RPV1000"; 
		selections="_trigger";
		
	}	
	

	if(jobnumber==30){
		sample_list.push_back("UDD300");
		_MC=true; 
		QCD=false;
		output_file_name="UDD300";
		selections="_trigger"; 
	}
	
	if(jobnumber==31) {
		for(int nJ=1; nJ<=4; nJ++) {
			TString tag=Form("w%dJets",nJ);
			sample_list.push_back(tag); 
		}
		_MC=true; 
		output_file_name="wJets"; 
		selections="_trigger";
	}
	
	if(jobnumber==32) {
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
		sample_list.push_back("QCD_600-800");
		sample_list.push_back("QCD_800-1000");
		
		
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
		output_file_name="allMC"; 
		selections="_trigger";
		
	}
	
}

void open_files(){
	TString dir="/eos/uscms/store/user/btcarlso/trees/March2/";

	initialize();
		
	//TString sample_list[]="wJetsht150","wJetsht200","wJetsht250","wJetsht300","wJetsht400",
	//TString sample_list[]={"RPV500"};
	
	int N_files=sample_list.size(); 
	cout << "N Files: " << N_files << endl; 

	weight=1;
	cout << "output file: " << output_file_name << endl; 
	TString output_FILE_NAME="output_file_"+output_file_name+".root";
	output_file=new TFile(output_FILE_NAME,"RECREATE"); 
	cout << "Analyzing samples: " << endl; 
	for (int iF=0; iF<N_files; iF++){
		cout << sample_list.at(iF) << " " << fileName[sample_list.at(iF)] << endl; 
	}
	cout << endl; 
	
	
	for (int iF=0; iF<N_files; iF++) {
		cout << sample_list[iF] << endl; 
		cout << "Ngen all procssed: " << Ngen[sample_list[iF]] << endl; 
		if(_MC) weight=xs[sample_list[iF]]*19600;//(Ngen[sample_list[iF]]);
		cout << "weight: " << weight << endl; 
		TString FN=dir+fileName[sample_list.at(iF)]+selections;

		analyze_file(FN);
	
	}
	
}

void CreateEfficiency(TString name, TString title){
	double jetPt[]={30,50,75,100,200,500,1000};
	const int N=sizeof(jetPt)/sizeof(double)-1;
	TEfficiency *peff = new TEfficiency(name, title,N,jetPt);
	
	effName[name]=peff;
	
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

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t nBinsX, const Double_t* xBins,
					 Int_t nBinsY,Double_t yLow, Double_t yUp)
{
	TH2F* h = new TH2F(name, title, nBinsX, xBins, nBinsY, yLow,yUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName2D[name] = h;
}

void CreateProfile(const char* name,   const char* title,
				   const char* xTitle, const char* yTitle,
				   Int_t       nBinsX, Float_t *xBins)
{
	TProfile* pr = new TProfile(name, title, nBinsX, xBins); 
	
	pr->GetXaxis()->SetTitle(xTitle);
	pr->GetYaxis()->SetTitle(yTitle);
	
	prName[name] = pr;
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

void CreateProfile(const char* name,   const char* title,
                   const char* xTitle, const char* yTitle,
                   Int_t       nBinsX, Double_t    xLow, Double_t xUp,
                   Int_t   nBinsY, Double_t    yLow, Double_t yUp, Double_t zLow, Double_t zUp  )
{
	TProfile2D* pr = new TProfile2D(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp, zLow, zUp);
	
	pr->GetXaxis()->SetTitle(xTitle);
	pr->GetYaxis()->SetTitle(yTitle);
	
	prName2D[name] = pr;
}


void writeHisto(){
   
    
	output_file->cd();
	int N1D=0, N2D=0, NP1D=0, NE1=0; 
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++)N1D++;
	//	for (std::map<TString,TProfile*>::iterator it=prName.begin(); it!=prName.end(); it++) NP1D++;
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++)N2D++; 
	for (std::map<TString,TEfficiency*>::iterator it=effName.begin(); it!=effName.end(); it++)NE1++;
	
	//	CreateHistogram("TProfile_names","","","",NP1D,0.,10.); 

	CreateHistogram("TH1F_names","","","",N1D,0.,10.); 
	CreateHistogram("TH2F_names","","","",N2D,0.,10.); 
	CreateHistogram("TEfficiency_names","","","",NE1,0,10); 
	
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		if(it->first=="TH1F_names" || it->first=="TH2F_names" || it->first=="TEfficiency_names" )continue;
		hName["TH1F_names"]->Fill(it->first,1);
	}
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++){
		if(it->first=="TH1F_names" || it->first=="TH2F_names" || it->first=="TEfficiency_names" )continue;
		hName["TH2F_names"]->Fill(it->first,1); 
	}
	
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		if(it->first.Contains("nJets_")){
			for(int i=1; i<=hName[it->first]->GetNbinsX(); i++){
				if(i<6)hName[it->first]->GetXaxis()->SetBinLabel(i,Form("%d",i)); 
				else hName[it->first]->GetXaxis()->SetBinLabel(6,"#geq 6"); 
			}
		}
		hName[it->first]->Write();
	}
	for (std::map<TString,TH2F*>::iterator it=hName2D.begin(); it!=hName2D.end(); it++) {
		hName2D[it->first]->Write();
	}

	for (std::map<TString,TProfile*>::iterator it=prName.begin(); it!=prName.end(); it++){
		prName[it->first]->Write();
	}
    for (std::map<TString,TProfile2D*>::iterator it=prName2D.begin(); it!=prName2D.end(); it++){
		prName2D[it->first]->Write();
	}
    
	if(MakeEfficiency){
		TString EFF_NAME="efficiency_file_"+output_file_name+".root"; 
		
		TFile *efficiency_file=new TFile(EFF_NAME,"RECREATE"); 
		efficiency_file->cd(); 
		for (std::map<TString,TEfficiency*>::iterator it=effName.begin(); it!=effName.end(); it++){
			effName[it->first]->Write();
			hName["TEfficiency_names"]->Fill(it->first,1); 
		}
		hName["TEfficiency_names"]->Write(); 
		efficiency_file->Close();
		delete efficiency_file; 
	}
	output_file->Close();
	delete output_file; 
}