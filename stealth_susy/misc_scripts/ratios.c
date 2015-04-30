/*
 *  ratios.c
 *  
 *
 *  Created by Benjamin Carlson on 9/20/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "ratios.h"

void ratios(){
	gROOT->SetBatch(); 
	open_files(); 
	get_ratios();
	
	for(int mu=1; mu<=2; mu++){
		plot_ratios("wJets",mu); 
		plot_ratios("ttbar",mu); 
		plot_ratios("dy",mu); 
		plot_ratios("singleMu",mu); 
		plot_ratios("allMC",mu);
	}
	int mu=-1;
	plot_ratios("wJets",mu); 
	plot_ratios("ttbar",mu); 
	plot_ratios("dy",mu); 
	plot_ratios("singleMu",mu); 
	plot_ratios("allMC",mu);
	
	cout << "Compare Ratios: " << endl; 
	compare_ratios(3,2,-1); 
	compare_ratios(4,3,-1);
	compare_ratios(5,4,-1);
	
	compare_ratios(4,3,1);
	compare_ratios(4,3,2);

	compare_ratios(5,4,1);
	compare_ratios(5,4,2);
	
	compare_ratios(5,3,1);
	compare_ratios(5,3,2);
	
	compare_ratios(6,4,1);
	compare_ratios(6,4,2);

	compare_ratios_muon(5,4);
	compare_ratios_muon(6,5);
	
	compare_ratios(6,5,-1); 
	compare_ratios(6,5,1); 
	compare_ratios(6,5,2); 

	//compare_ratios(7,6,1);
	
	cout << "Write canvases: "<< endl; 
	
	TFile *output_file = new TFile("plots_ratios.root","RECREATE"); 
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

void get_ratios(){
	for(int num=3; num<=6; num++){
		TString file="wJets";
		for(int den = 2; den<=6; den++){
			TString name=Form("Ratio%d-%d",num,den); 
			GetGraph(name,file); 
		}
		 file="dy";
		for(int den = 2; den<=6; den++){
			TString name=Form("Ratio%d-%d",num,den); 
			GetGraph(name,file); 
		}
		
		
		file="ttbar"; 
		for(int den = 2; den<=6; den++){
			TString name=Form("Ratio%d-%d",num,den); 
			GetGraph(name,file); 
		}
		file="singleMu"; 
		for(int den = 2; den<=6; den++){
			TString name=Form("Ratio%d-%d",num,den); 
			GetGraph(name,file); 
		}
		file="allMC"; 
		for(int den = 2; den<=6; den++){
			TString name=Form("Ratio%d-%d",num,den); 
			GetGraph(name,file); 
		}
		
	}
	for(int mu=1; mu<=2; mu++){
		for(int num=3; num<=6; num++){
			TString file="wJets";
			for(int den = 2; den<=6; den++){
				TString name=Form("Ratio%d-%d_%dMu",num,den,mu); 
				GetGraph(name,file); 
			}
			file="dy"; 
			for(int den = 2; den<=6; den++){
				TString name=Form("Ratio%d-%d_%dMu",num,den,mu); 
				GetGraph(name,file); 
			}
			
			file="ttbar"; 
			for(int den = 2; den<=6; den++){
				TString name=Form("Ratio%d-%d_%dMu",num,den,mu); 
				GetGraph(name,file); 
			}
			file="singleMu"; 
			for(int den = 2; den<=6; den++){
				TString name=Form("Ratio%d-%d_%dMu",num,den,mu); 
				GetGraph(name,file); 
			}
			file="allMC"; 
			for(int den = 2; den<=6; den++){
				TString name=Form("Ratio%d-%d_%dMu",num,den,mu); 
				GetGraph(name,file); 
			}
			
		}
	}
}

void setcolor(TGraphAsymmErrors *gr, int num){

	gr->SetMarkerStyle(8); 
	gr->SetMarkerSize(0.5);
	gr->SetLineColor(kBlack);
	if(num==6) gr->SetMarkerColor(kRed);
	if(num==4) gr->SetMarkerColor(kBlue);
	if(num==5) gr->SetMarkerColor(kGreen); 
	if(num==3) gr->SetMarkerColor(kBlack); 
	
}

double avg (TGraphAsymmErrors *gr){

	double average=0; 
	double w=0; 
	
	for(int i=0; i<gr->GetN(); i++){
		if(gr->GetX()[i]>1000) {
			double sigma=(gr->GetEYhigh()[i]+gr->GetEYlow()[i])/2; 
			average+=gr->GetY()[i]/sigma; 
			w+=1./sigma;  
		}
	}
	average=average/w; 
	return average; 
}

void scale_graph(TGraphAsymmErrors *gr, double SF){
	//Scale TGraphAsymmErrors 
	
	for (int i=0;i<gr->GetN();i++){
		gr->GetY()[i] *= SF;
		gr->GetEYhigh()[i]*=SF; 
		gr->GetEYlow()[i]*=SF;
	}
}

void compare_ratios_muon(int num, int den){
	cout << "Ratio comparison of 1 and 2 muons: " << num << "-" << den << endl; 
	TString canvas_name=Form("Compare_R%d-%d_muon",num,den);
	
	CreateCanvas(canvas_name,canvas_name,600,600);
	CName[canvas_name]->cd();

	int mu=1; 
	TString name_mu1=Form("Ratio%d-%d_%dMuallMC",num,den,mu);
	mu=2; 
	TString name_mu2=Form("Ratio%d-%d_%dMuallMC",num,den,mu);
	cout << "created canvas: "<< endl; 
	
	TGraphAsymmErrors *grmu1 = (TGraphAsymmErrors*)grName[name_mu1]->Clone(name_mu1+"_norm"); 
	TGraphAsymmErrors *grmu2 = (TGraphAsymmErrors*)grName[name_mu2]->Clone(name_mu2+"_norm"); 
	cout << "clone graphs: "<< endl; 
	double SMC=1./avg(grmu1); 
	double Sdata=1./avg(grmu2); 
	
	
	setcolor(grmu1,6); 
	setcolor(grmu2,4);
	
	TLegend *L = new TLegend(0.15,0.15,0.5,0.35);
	L->SetFillColor(10); 
	L->AddEntry(grmu1,"#mu=1", "LEP"); 
	L->AddEntry(grmu2,"#mu=2", "LEP"); 
	
	cout << "Draw graphs: "<< endl; 
	grmu1->Draw("ap "); 
	grmu1->GetXaxis()->SetRangeUser(300,3000); 
	grmu1->GetYaxis()->SetRangeUser(0,1.5*grmu1->Eval(1200)); 
	grmu1->GetYaxis()->SetTitle(Form("%djets/%djets",num,den)); 
	grmu1->DrawClone("ap"); 
	grmu2->DrawClone("p same");
	L->DrawClone("same");
	
}

void compare_ratios(int num, int den,int mu){
	cout << "Ratio comparison: " << num << "-" << den << endl; 
	TString canvas_name=Form("Compare_R%d-%d_%dMu",num,den,mu);
	if(mu==-1)canvas_name=Form("Compare_R%d-%d",num,den);
	
	CreateCanvas(canvas_name,canvas_name,600,600);
	CName[canvas_name]->cd();
	
	TString name=Form("Ratio%d-%d_%dMu",num,den,mu);
	if(mu==-1)name=Form("Ratio%d-%d",num,den);
	TString namedata=name+"singleMu"; 
	TString nameMC=name+"allMC"; 
	

	TGraphAsymmErrors *grMC = (TGraphAsymmErrors*)grName[nameMC]->Clone(nameMC+"_norm"); 
	TGraphAsymmErrors *grdata = (TGraphAsymmErrors*)grName[namedata]->Clone(namedata+"_norm"); 
	double SMC=1./avg(grMC); 
	double Sdata=1./avg(grdata); 
	
	//scale_graph(grMC,SMC); 
	//scale_graph(grdata,Sdata); 
	

	setcolor(grMC,6); 
	setcolor(grdata,4);

	TLegend *L = new TLegend(0.15,0.15,0.5,0.35);
	L->SetFillColor(10); 
	L->AddEntry(grdata,"Single #mu: data", "LEP"); 
	L->AddEntry(grMC,"W+Jets + ttbar + dy", "LEP"); 
	
	grdata->Draw("ap "); 
	grdata->GetXaxis()->SetRangeUser(300,3000); 
	grdata->GetYaxis()->SetRangeUser(0,1.5*grdata->Eval(1200)); 
	grdata->GetYaxis()->SetTitle(Form("%djets/%djets",num,den)); 
	grdata->DrawClone("ap"); 
	grMC->DrawClone("p same");
	L->DrawClone("same");

}

void plot_ratios(TString file, int mu){

	TString canvas_name=file+Form("_Ratio_%dMu",mu);
	if(mu==-1)canvas_name=file+"_Ratio"; 
	CreateCanvas(canvas_name,canvas_name,600,600); 
	CName[canvas_name]->Divide(2,2); 
	
	for(int den=2; den<=5; den++){
		CName[canvas_name]->cd(den-1); 
		bool first=true; 
		for(int num=3; num<=6; num++){
			TString name=Form("Ratio%d-%d_%dMu",num,den,mu);
			if(mu==-1) name=Form("Ratio%d-%d",num,den);
			name=name+file; 
			setcolor(grName[name],num); 

			if(num>den){
				if(first){
					grName[name]->Draw("ap"); 
					grName[name]->GetXaxis()->SetRangeUser(300,3000); 
					grName[name]->GetYaxis()->SetRangeUser(0,2*grName[name]->Eval(2000)); 
					grName[name]->GetYaxis()->SetTitle(Form("njets/%d",den)); 
					grName[name]->Draw("ap");
					first=false; 
				}
				else grName[name]->Draw("p same");
			}
		}
	}
	
}

void open_files(){
	//Open MC files
	
	//	open_file("QCD");
	open_file("wJets");
	open_file("dy");
	open_file("ttbar");
	open_file("singleMu");
	open_file("allMC"); 
	
}

void GetGraph(TString graph_name, TString file){
	TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fName[file]->FindObjectAny(graph_name);
	grName[graph_name+file]=gr; 
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