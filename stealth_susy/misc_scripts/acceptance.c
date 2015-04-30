/*
 *  acceptance.c
 *  
 *
 *  Created by Benjamin Carlson on 10/29/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "acceptance.h"

void acceptance(){
	gROOT->SetBatch(); 
	open_files();
	load_histograms();
	compute_acceptance(1); 
	compute_acceptance(2); 

	string signal[]={"RPV300","RPV400","RPV500","RPV600","RPV700","RPV800","RPV900","RPV1000"};
	for(int nJ=3; nJ<=9; nJ++){
		for(int i=0; i<8; i++){
			compute_sensitivity(1,nJ,signal[i]); 
			if(nJ>7) continue; 
			compute_sensitivity(2,nJ,signal[i]); 
		}
	}
	sensitivity_table(1);
	sensitivity_table(2);
	plot_sensitivity();
	plot_contamination();
	compute_contamination_Z();
	TFile *output_file = new TFile("plots_acceptance.root","RECREATE"); 
	output_file->cd();
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		CName[it->first]->Print("/uscms_data/d3/btcarlso/Figures_SUSYAN/"+it->first+".png");
	}
}

void open_file(TString name){
	TString file_name="output_file_"+name+".root"; 
	
	TFile *f = new TFile(file_name,"READ");
	if(f->IsOpen()!=1) {
		cout << "File: " << file_name << " Failed to open!" << endl; 
		return; 
	}
	fName[name]=f; 
	
}

void sensitivity_table(int mu){
	cout << "sensitivity table: " << mu << endl; 
	ofstream output; 
	output.open(Form("sensitivity_%dMu.tex",mu)); 
	
	output <<"\\begin{tabular}{"; 
	int nJetmax=0; 
	if(mu==1) nJetmax=9; 
	if(mu==2) nJetmax=7; 
	for(int i=4; i<=nJetmax; i++) output << "c";
	output << "}\\\\ \\hline" << endl; 
	output << "& ";
	for(int i=4; i<=nJetmax; i++) output << Form(" $\\N_{jets}=%d$ &",i);
	output << "\\\\ \\hline" << endl; 
	for(int irpv=300; irpv<=1000; irpv+=100){
		string signal=Form("RPV%d",irpv); 
		output << irpv << " & "; 
		for(int nJ=4; nJ<=nJetmax; nJ++){
			double st_thresh=hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->GetBinCenter(hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->GetMaximumBin()); 
			double s_b=hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->GetMaximum(); 
			output << st_thresh << " & " << s_b << " & "; 
		}
		output << "\\\\" << endl; 
	}
	output << "\\hline" << endl; 
	output.close(); 
}

void compute_acceptance(int mu){
	cout << "mu: " << mu << endl; 
	ofstream output; 
	output.open(Form("acceptance_%dMu.tex",mu)); 
	int N=6; 
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "Mass [GeV] & 4 jets & 5 jets & 6 jets & >=7 & >=8 & >=9 \\\\ \\hline " << endl; 
	output << "300 & "; 
	for(int nJ=4; nJ<=9; nJ++){
		string mass_point="RPV300"; 
		TString hist_name=Form("st_nJets%d_%dMu_uW_50GeV%s",nJ,mu,mass_point.c_str()); 
		double acc=100*hName[hist_name]->Integral(1,hName[hist_name]->FindBin(3000))/39198;
		if(nJ<9)output << acc << " & "; 
		else output << acc; 
	}
	output << "\\\\ " << endl; 
	
	output << "500 & "; 
	for(int nJ=4; nJ<=9; nJ++){
		string mass_point="RPV500"; 
		TString hist_name=Form("st_nJets%d_%dMu_uW_50GeV%s",nJ,mu,mass_point.c_str()); 
		double acc=100*hName[hist_name]->Integral(1,hName[hist_name]->FindBin(3000))/34397;
		if(nJ<9)output << acc << " & "; 
		else output << acc; 
	}
	output << "\\\\ " << endl; 
	
	output << "1000 & "; 
	for(int nJ=4; nJ<=9; nJ++){
		string mass_point="RPV1000"; 
		TString hist_name=Form("st_nJets%d_%dMu_uW_50GeV%s",nJ,mu,mass_point.c_str()); 
		double acc=100*hName[hist_name]->Integral(1,hName[hist_name]->FindBin(3000))/39996;
		if(nJ<9)output << acc << " & "; 
		else output << acc; 
	}
	output << "\\\\ \\hline" << endl; 
	output << "\\end{tabular}" << endl; 
	output.close();
}

void compute_contamination_Z(){

	
	cout << "RPV300: " << endl; 
	string signal="RPV300"; 
	for(int j=1; j<10; j++){
		double jet=0.5+static_cast<double>(j); 
		int nJ=hName[Form("nJets_3Mass%s",signal.c_str())]->FindBin(jet);
		double sig=hName[Form("nJets_3Mass%s",signal.c_str())]->GetBinContent(nJ);
		double bkg=hName[Form("nJets_3Mass%s","allMC")]->GetBinContent(nJ);
		cout << "jet: " << j << " s/s+b " << sig/(sig+bkg) << endl; 

	}
	
	signal="RPV500"; 
	cout << "RPV500: " << endl; 
	for(int j=1; j<10; j++){
		double jet=0.5+static_cast<double>(j); 
		int nJ=hName[Form("nJets_3Mass%s",signal.c_str())]->FindBin(jet);
		double sig=hName[Form("nJets_3Mass%s",signal.c_str())]->GetBinContent(nJ);
		double bkg=hName[Form("nJets_3Mass%s","allMC")]->GetBinContent(nJ);
		cout << "jet: " << j << " s/s+b " << sig/(sig+bkg) << endl; 

	}
	
}

void plot_contamination(){

	CreateCanvas("contamination","",Cx,Cy); 
	CName["contamination"]->cd();
	
	TLegend Leg(0.2,0.6,0.6,0.9);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetLineWidth(0);
	
	int nJ=5; 
	int mu=2; 
	string signal="RPV300"; 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineWidth(2); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kBlack); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetAxisRange(0,2500);
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetMinimum(0);
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw(); 	
	
	Leg.AddEntry(hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,300),"L"); 
	
	nJ=4; 
	mu=2; 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kBlue); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 	
		
	Leg.AddEntry(hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,300),"L"); 

	nJ=3; 
	mu=2; 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kRed); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 
	Leg.AddEntry(hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,300),"L"); 

	nJ=4; 
	mu=2; 
	signal="RPV500"; 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kGreen); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 
	Leg.AddEntry(hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,500),"L"); 

	nJ=4; 
	mu=2; 
	signal="RPV1000"; 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kBlue); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineStyle(kDashed); 
	hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 	
	Leg.AddEntry(hName[Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,1000),"L"); 

	Leg.DrawClone("same"); 
}

void plot_sensitivity(){
	CreateCanvas("s_b_plot","",Cx,Cy);
	CName["s_b_plot"]->cd();
	
	TLegend Leg(0.4,0.5,0.9,0.85);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetLineWidth(0);
	
	int nJ=7; 
	int mu=2; 
	string signal="RPV300"; 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineWidth(2); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kBlack); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetAxisRange(0,2500);
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw(); 	
	
	Leg.AddEntry(hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,300),"L"); 

	
	nJ=7;
	mu=1; 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kRed); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 	
	Leg.AddEntry(hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,300),"L"); 

	
	nJ=4; 
	mu=2; 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kBlue); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 	
	
	Leg.AddEntry(hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,300),"L"); 

	
	nJ=7; 
	mu=2; 
	signal="RPV500"; 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kGreen); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineStyle(kDashed); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 
	Leg.AddEntry(hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("%d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,500),"L"); 

	nJ=7; 
	mu=2; 
	signal="RPV1000"; 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineColor(kBlue); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->SetLineStyle(kDashed); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Scale(10); 
	hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())]->Draw("same"); 	
	Leg.AddEntry(hName[Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())],Form("#times 10, %d-jets,#mu=%d,m_{#tilde{t}}=%d GeV",nJ,mu,1000),"L"); 

	Leg.DrawClone("same"); 
	
	

}

void compute_sensitivity(int mu, int nJ, string signal){
	cout << "sensitivity (mu): " << mu << " nJ: " << nJ << " file: " <<signal << endl; 
	string file="allMC"; 
	TH1F *s_b=(TH1F*)hName[Form("st_nJets%d_%dMu_50GeV%s",nJ,mu,file.c_str())]->Clone(Form("sb_nJets%d_%dMu%s",nJ,mu,signal.c_str())); 
	TH1F *c=(TH1F*)hName[Form("st_nJets%d_%dMu_50GeV%s",nJ,mu,file.c_str())]->Clone(Form("contamination_nJets%d_%dMu%s",nJ,mu,signal.c_str())); 
	c->GetYaxis()->SetTitle("Signal Contamination"); 
	c->GetXaxis()->SetTitle("S_{T}^{thresh} [GeV]"); 
	for(int i=1; i<=s_b->GetNbinsX();i++){
		s_b->SetBinContent(i,0); 
		s_b->SetBinError(i,0);
		
		c->SetBinContent(i,0);
		c->SetBinError(i,0);
		
	}
	s_b->GetYaxis()->SetTitle("S/#sqrt{B}");
	s_b->GetXaxis()->SetTitle("S_{T}^{thresh} [GeV]"); 
	int Nbins=hName[Form("st_nJets%d_%dMu_50GeV%s",nJ,mu,file.c_str())]->GetNbinsX(); 
	for(int ist=1; ist<=Nbins; ist++){
		if(s_b->GetBinCenter(ist)>4000) continue; 
		double bkg=hName[Form("st_nJets%d_%dMu_50GeV%s",nJ,mu,file.c_str())]->Integral(ist,Nbins); 
		double sig=hName[Form("st_nJets%d_%dMu_50GeV%s",nJ,mu,signal.c_str())]->Integral(ist,Nbins);
		
		if((sig+bkg)<0.0000000000000001) continue; 
		
		c->SetBinContent(ist,sig/(sig+bkg)); 
		c->SetBinError(ist,0);
		
		if(bkg<0.000000000000001) continue; 
		
		//cout << "st: " << s_b->GetBinCenter(ist) << " " << sig/TMath::Sqrt(bkg) << endl; 
		
		s_b->SetBinContent(ist,sig/TMath::Sqrt(bkg)); 
		s_b->SetBinError(ist,0); 
	}
	double st_opt = s_b->GetBinCenter(s_b->GetMaximumBin());
	
	//cout << "Optimum: " << st_opt << " S/sqrt(B): " << s_b->GetBinContent(s_b->GetMaximumBin())<< endl; 
	
	hName[s_b->GetName()]=s_b; 
	hName[c->GetName()]=c; 
	
}


void open_files(){
	
	for(int irpv=300; irpv<=1000; irpv+=100)open_file(Form("RPV%d",irpv)); 

	open_file("allMC"); 
}

void load_histograms(){
	TString name="nJets_3Mass"; 
	for(int irpv=300; irpv<=1000; irpv+=100)GetHistogram(name,Form("RPV%d",irpv)); 
	
	GetHistogram(name,"allMC");
	
	for(int nJ=3; nJ<=9; nJ++){
		name=Form("st_nJets%d_2Mu_uW_50GeV",nJ);
		for(int irpv=300; irpv<=1000; irpv+=100)GetHistogram(name,Form("RPV%d",irpv)); 

		
		name=Form("st_nJets%d_1Mu_uW_50GeV",nJ);
		for(int irpv=300; irpv<=1000; irpv+=100)GetHistogram(name,Form("RPV%d",irpv)); 

		
		name=Form("st_nJets%d_2Mu_50GeV",nJ);
		for(int irpv=300; irpv<=1000; irpv+=100)GetHistogram(name,Form("RPV%d",irpv)); 

		GetHistogram(name,"allMC");

		
		name=Form("st_nJets%d_1Mu_50GeV",nJ);
		for(int irpv=300; irpv<=1000; irpv+=100)GetHistogram(name,Form("RPV%d",irpv)); 

		GetHistogram(name,"allMC");

	}
}

void GetHistogram(TString hist_name, TString file){
	TH1F *h=(TH1F*)fName[file]->FindObjectAny(hist_name);
	h->SetStats(kFALSE);
	hName[hist_name+file]=h; 
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	createC->SetFillColor(10); 
	CName[Name]=createC;
}