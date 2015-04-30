/*
 *  set_limit.c
 *  
 *
 *  Created by Benjamin Carlson on 1/15/14.
 *  Copyright 2014 Carnegie Mellon University. All rights reserved.
 *
 */

#include "set_limit.h"

void set_limit(){
	gROOT->SetBatch();
	GetHistograms();
	bookGraphs();
	create_theory();
	read_data();
	
	for(int nJ=5; nJ<=8; nJ++){
		for(int mu=1; mu<=2; mu++){
			if(mu==1 && nJ==5) continue; 
			//cout << "draw nJ: " << nJ << " mu:  " << mu << endl; 
			draw_plots( nJ,  mu);
	}
	}
	
	draw_plots();
	
	TFile *limit_plots = new TFile("limit_plots.root","RECREATE"); 
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		TString name="/uscms_data/d3/btcarlso/Figures_SUSYAN/"; 
		name=name+it->first+".pdf"; 
		CName[it->first]->Print(name);
	}
	for (std::map<TString,TGraphAsymmErrors*>::iterator it=grName.begin(); it!=grName.end(); it++) {
		grName[it->first]->Write();
		TString name="/uscms_data/d3/btcarlso/Figures_SUSYAN/"; 
		name=name+it->first+".pdf"; 
		//CName[it->first]->Print(name);
	}
	
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
}

void draw_plots(){
	CreateCanvas("Limit_Plot_combiend","",600,600); 
	CName["Limit_Plot_combiend"]->cd();
	gPad->SetLogy(); 
	grName["stop_xs"]->SetLineColor(kRed);
	
	grName["xs_limit_combined"]->SetLineColor(kBlack);
	grName["xs_obs_combined"]->SetLineColor(kBlack);
	grName["xs_obs_combined"]->SetLineStyle(kDashed);
	
	grName["xs_limit_1sigma_combined"]->SetFillStyle(1001);
	grName["xs_limit_2sigma_combined"]->SetFillStyle(1001);
	
	grName["xs_limit_1sigma_combined"]->SetFillColor(419);
	grName["xs_limit_2sigma_combined"]->SetFillColor(5);
	grName["xs_limit_1sigma_combined"]->Draw("a3");
	grName["xs_limit_1sigma_combined"]->GetXaxis()->SetRangeUser(300,1000);
	grName["xs_limit_1sigma_combined"]->GetYaxis()->SetRangeUser(5,1.5*grName["stop_xs"]->Eval(300));
	grName["xs_limit_1sigma_combined"]->GetXaxis()->SetTitle("M_{#tilde{t}} [GeV]"); 
	grName["xs_limit_1sigma_combined"]->GetYaxis()->SetTitle("#sigma [fb]"); 
	
	grName["xs_limit_1sigma_combined"]->Draw("a3");
	grName["xs_limit_combined"]->Draw("l same");
	grName["xs_obs_combined"]->Draw("l same");
	
	//grName[Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)]->Draw("3 same");
	grName["stop_xs"]->Draw("c same"); 
	
	
}

void draw_plots(int nJ, int mu){
	CreateCanvas(Form("Limit_Plot_nJets%d_%dMu",nJ,mu),"",600,600); 
	CName[Form("Limit_Plot_nJets%d_%dMu",nJ,mu)]->cd();
	gPad->SetLogy(); 
	grName["stop_xs"]->SetLineColor(kRed);
	if(grName[Form("xs_limit_nJets%d_%dMu",nJ,mu)]->GetN()==0) cout << "nJ: " << nJ << " mu: " << mu << endl; 

	grName[Form("xs_limit_nJets%d_%dMu",nJ,mu)]->SetLineColor(kBlack);
	grName[Form("xs_obs_limit_nJets%d_%dMu",nJ,mu)]->SetLineColor(kBlack);
	grName[Form("xs_obs_limit_nJets%d_%dMu",nJ,mu)]->SetLineStyle(kDashed);

	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->SetFillStyle(1001);
	grName[Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)]->SetFillStyle(1001);

	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->SetFillColor(419);
	grName[Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)]->SetFillColor(5);
	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->Draw("a3");
	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->GetXaxis()->SetRangeUser(300,1000);
	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->GetYaxis()->SetRangeUser(5,1.5*grName["stop_xs"]->Eval(300));
	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->GetXaxis()->SetTitle("M_{#tilde{t}} [GeV]"); 
	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->GetYaxis()->SetTitle("#sigma [fb]"); 

	grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->Draw("a3");
	grName[Form("xs_limit_nJets%d_%dMu",nJ,mu)]->Draw("l same");
	grName[Form("xs_obs_limit_nJets%d_%dMu",nJ,mu)]->Draw("l same");

	//grName[Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)]->Draw("3 same");
	grName["stop_xs"]->Draw("c same"); 


}

void create_theory(){
	double mass[]={300,400,500,600,700,800,900,1000};
	int N=sizeof(mass)/sizeof(double);
	double stopxs[]={1.99608,0.35683,0.0855847,0.0248009,0.0081141,0.00289588,0.00109501,0.000435488};
	for(int iP=0; iP<N; iP++){
		grName["stop_xs"]->SetPoint(iP,mass[iP],stopxs[iP]*1000); 
	}
}

void CreateGraph(TString name){
	TGraphAsymmErrors *gr = new TGraphAsymmErrors(); 
	gr->SetName(name); 
	grName[name]=gr;
}

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	h->SetStats(kFALSE);
	hName[name] = h;
}

void bookGraphs(){
	CreateHistogram("Mass_Bins","","","",7,250,1050);
	for(int nJ=5; nJ<=8; nJ++){
		for(int mu=1; mu<=2; mu++){
			if(mu==1 && nJ==5) continue; 
			CreateGraph(Form("xs_limit_nJets%d_%dMu",nJ,mu)); 
			CreateGraph(Form("xs_obs_limit_nJets%d_%dMu",nJ,mu));
			CreateGraph(Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)); 
			CreateGraph(Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)); 
		}
	}
	CreateGraph("xs_limit_1sigma_combined");
	CreateGraph("xs_limit_2sigma_combined");

	CreateGraph("xs_limit_combined");
	CreateGraph("xs_obs_combined");

	CreateGraph("stop_xs"); 

	
	
}

void read_data(){
	ifstream input; 
	input.open("/uscms_data/d3/btcarlso/LimitSetting/CMSSW_6_1_1/src/HiggsAnalysis/CombinedLimit/limit_summary.txt"); 
	int mass; 
	int nJ;
	int mu;
	
	double obs;
	double exp;
	double expp1;
	double expp2;
	double expm1;
	double expm2;
	
	
	bool print=false;
	while(!input.eof()){
		input >> mass >> nJ >> mu >> obs >> exp >> expm1 >> expm2 >> expp1 >> expp2; 
		string line; 
		getline(input,line);
		
		if(print){
			cout << "mass: " << mass << " nJ : " << nJ << " mu: " << mu; 
			cout << " obs: " << obs << " exp: " << exp << " expm1: " << expm1 << " expm2: " << expm2; 
			cout << " expp1: " << expp1 << " expp2: " << expp2 << endl; 
		}
		
		if(mass>1000 || mass<300)continue;
		if(nJ>8 || nJ<5) continue;
		if(mu>2 || mu<1) continue; 
		
		int iP=hName["Mass_Bins"]->FindBin(mass)-1; 
		if(print)cout << "iP: " << iP << endl; 
		double ex=50;
	
		if(exp>0){
			double stopXS=grName["stop_xs"]->GetY()[iP];
			double xs=grName["stop_xs"]->GetY()[iP]*exp;
			double obsXS=stopXS*obs;
			if(print) cout << "ip: " << iP << " " << xs << endl;

			double eyl_1S=stopXS*(exp-expm1);
			double eyh_1S=stopXS*(expp1-exp);
			
			double eyl_2S=stopXS*(exp-expm2);
			double eyh_2S=stopXS*(expp2-exp);
			if(print)cout << "mass: " << mass << " nJ " << nJ << " mu " << mu << endl;  
			if(print)cout << "xs: " << xs << " +/- " << eyl_1S << " " << eyh_1S << endl; 

			grName[Form("xs_obs_limit_nJets%d_%dMu",nJ,mu)]->SetPoint(iP,mass,obsXS); 
			grName[Form("xs_limit_nJets%d_%dMu",nJ,mu)]->SetPoint(iP,mass,xs); 
			grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->SetPoint(iP,mass,xs); 
			grName[Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)]->SetPoint(iP,mass,xs); 
			grName[Form("xs_limit_1sigma_nJets%d_%dMu",nJ,mu)]->SetPointError(iP,ex,ex,eyl_1S,eyh_1S); 
			grName[Form("xs_limit_2sigma_nJets%d_%dMu",nJ,mu)]->SetPointError(iP,ex,ex,eyl_2S,eyh_2S); 
		}
		
	}
	
	input.close();
	input.open("/uscms_data/d3/btcarlso/LimitSetting/CMSSW_6_1_1/src/HiggsAnalysis/CombinedLimit/limit_summary_combined.txt"); 

	
	while(!input.eof()){
		input >> mass >> obs >> exp >> expm1 >> expm2 >> expp1 >> expp2; 
		string line; 
		getline(input,line);
		
		if(print){
			cout << "mass: " << mass;
			cout << " obs: " << obs << " exp: " << exp << " expm1: " << expm1 << " expm2: " << expm2; 
			cout << " expp1: " << expp1 << " expp2: " << expp2 << endl; 
		}
		
		if(mass>1000 || mass<300)continue;
		
		int iP=hName["Mass_Bins"]->FindBin(mass)-1; 
		if(print)cout << "iP: " << iP << endl; 
		double ex=50;
		
		if(exp>0){
			double stopXS=grName["stop_xs"]->GetY()[iP];
			double xs=grName["stop_xs"]->GetY()[iP]*exp;
			double obsXS=stopXS*obs;
			if(print) cout << "ip: " << iP << " " << xs << endl;
			
			double eyl_1S=stopXS*(exp-expm1);
			double eyh_1S=stopXS*(expp1-exp);
			
			double eyl_2S=stopXS*(exp-expm2);
			double eyh_2S=stopXS*(expp2-exp);
			if(print)cout << "mass: " << mass << " nJ " << nJ << " mu " << mu << endl;  
			if(print)cout << "xs: " << xs << " +/- " << eyl_1S << " " << eyh_1S << endl; 
		
			grName["xs_obs_combined"]->SetPoint(iP,mass,obsXS); 
			grName["xs_limit_combined"]->SetPoint(iP,mass,xs); 
			grName["xs_limit_1sigma_combined"]->SetPoint(iP,mass,xs); 
			grName["xs_limit_2sigma_combined"]->SetPoint(iP,mass,xs); 
			grName["xs_limit_1sigma_combined"]->SetPointError(iP,ex,ex,eyl_1S,eyh_1S); 
			grName["xs_limit_2sigma_combined"]->SetPointError(iP,ex,ex,eyl_2S,eyh_2S); 
		}
		
	}
	
}

void GetHistograms(){
	TFile *input_file = new TFile("susy_histograms.root","READ"); 
	TH1F *TH1F_names=(TH1F*)input_file->FindObjectAny("TH1F_names"); 
	
	
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names") continue; 
		//cout <<" filename: "<< file <<  " hist_name: "<< name << endl; 
		
		TH1F *h=(TH1F*)input_file->FindObjectAny(name); 
		h->SetStats(kFALSE); 
		h->SetName(h->GetName()); 
		hName[h->GetName()]=h; 
	}
	delete TH1F_names; 
	//cout << "end of GetHistograms() " << endl; 
}
