/*
 *  st_stack.c
 *  
 *
 *  Created by Benjamin Carlson on 7/18/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "st_stack.h"

void st_stack(){
	gROOT->SetBatch();
	open_files();
	plot_npv();
	//fill_stack("st", 1);
	
	for(int nJ=1;nJ<=10; nJ++){
		bool blind=false; 
		
		fill_stack(Form("st_nJets%d_3Mass_2Mu",nJ),blind); 
		
		if(nJ>=5)blind=true; 
		fill_stack(Form("st_nJets%d_4Mass_2Mu",nJ),blind); 
		blind=false; 
		
		if(nJ<=6) blind=true; 
		fill_stack(Form("st_nJets%d_1Mu_GeV",nJ),blind); 
		
		blind=false; 
		if(nJ<=5) blind=true; 
		fill_stack(Form("st_nJets%d_2Mu_GeV",nJ),blind); 
	

	}

	//fill_stack("ht");
	
	fill_stack_inclusive("nJets_1Mu",true);
	fill_stack_inclusive("nJets_2Mu",false);

	fill_stack_inclusive("nJets_3Mass",false);
	fill_stack_inclusive("nJets_4Mass",true);

	
	fill_stack_inclusive("nJets_st500_1000_1Mu",true);
	fill_stack_inclusive("nJets_st500_1000_2Mu",true);
	fill_stack_inclusive("h_met",false);
	fill_stack_inclusive("h_mumu",false);
	for(int nJ=1; nJ<=8; nJ++){
		bool blind=false; 
		if(nJ>=5)blind=true; 
		fill_stack_inclusive(Form("h_mumu_nJet%d",nJ),blind);
	}
	
	
	fill_stack_inclusive("ht_total",false);
	fill_stack_inclusive("st_total",false);
	plot_signal();
	
	TFile *output_file = new TFile("plots_stacks.root","RECREATE"); 
	output_file->cd();
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		CName[it->first]->Print("/uscms_data/d3/btcarlso/Figures_SUSYAN/"+it->first+".png");
	}
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
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

void open_files(){
	//Open MC files
	
	open_file("QCD");
	for(int irpv=300; irpv<=1000; irpv+=100){
		open_file(Form("RPV%d",irpv));
	}
	open_file("stealth_ww500_300");
	open_file("stealth_ww1400_300");
	open_file("wJets");
	open_file("dy");
	open_file("ttbar");
	open_file("singleTop");
	open_file("diboson"); 
	open_file("singleMu");
	open_file("allMC"); 

}

void CreateStack(TString Name,TString Title){
	THStack *hstack = new THStack(Name,Title); 
	stackName[Name]=hstack;
}

void marker_style(TH1F *h){
	h->SetMarkerStyle(7); 
	h->SetMarkerSize(1); 
	
}

void plot_signal(){
	
	TString hist_name="h_mumu"; 
	GetHistogram(hist_name,"RPV300"); 
	GetHistogram(hist_name,"RPV500"); 
	GetHistogram(hist_name,"RPV1000"); 
	
	TLegend Lhmumu(0.71,0.71,0.87,0.87);
	Lhmumu.SetFillColor(10);
	Lhmumu.SetLineColor(10);
	Lhmumu.SetLineWidth(0); 
	
	CreateCanvas("signal_mumu","",700,600);
	CName["signal_mumu"]->cd();
	
	hName[hist_name+"RPV300"]->SetLineColor(kRed);
	hName[hist_name+"RPV500"]->SetLineColor(kBlue);
	hName[hist_name+"RPV1000"]->SetLineColor(kBlack);
	
	hName[hist_name+"RPV500"]->Scale(5);
	hName[hist_name+"RPV1000"]->Scale(20);
	
	Lhmumu.AddEntry(hName[hist_name+"RPV300"],"m_{#tilde{t}=300 GeV}","LEP");
	Lhmumu.AddEntry(hName[hist_name+"RPV500"],"5 #times m_{#tilde{t}=500 GeV}","LEP");
	//Lhmumu.AddEntry(hName[hist_name+"RPV1000"],"20 #times m_{#tilde{t}=1000 GeV}","LEP");
	
	hName[hist_name+"RPV300"]->Draw();
	hName[hist_name+"RPV500"]->Draw("same");
	//hName[hist_name+"RPV1000"]->Draw("");
	Lhmumu.DrawClone(); 
	
	 hist_name="nJets"; 
	GetHistogram(hist_name,"RPV300"); 
	GetHistogram(hist_name,"RPV500"); 
	GetHistogram(hist_name,"RPV1000"); 
	
	TLegend L1(0.71,0.71,0.87,0.87);
	L1.SetFillColor(10);
	L1.SetLineColor(10);
	L1.SetLineWidth(0); 
	
	CreateCanvas("signal_njets","",700,600);
	CName["signal_njets"]->cd();
	gPad->SetLogy(); 
	
	hName[hist_name+"RPV300"]->SetLineColor(kRed);
	hName[hist_name+"RPV500"]->SetLineColor(kBlue);
	hName[hist_name+"RPV1000"]->SetLineColor(kBlack);
	
	L1.AddEntry(hName[hist_name+"RPV300"],"m_{#tilde{t}=300 GeV}","LEP");
	L1.AddEntry(hName[hist_name+"RPV500"],"m_{#tilde{t}=500 GeV}","LEP");
	L1.AddEntry(hName[hist_name+"RPV1000"],"m_{#tilde{t}=1000 GeV}","LEP");

	hName[hist_name+"RPV300"]->SetMinimum(0.001);
	hName[hist_name+"RPV300"]->SetMaximum(20000);
	hName[hist_name+"RPV300"]->Draw();
	hName[hist_name+"RPV500"]->Draw("same");
	hName[hist_name+"RPV1000"]->Draw("same");
	L1.DrawClone(); 

	hist_name="st_nJets4"; 
	GetHistogram(hist_name,"RPV300"); 
	GetHistogram(hist_name,"RPV500"); 
	GetHistogram(hist_name,"RPV1000"); 
	
	hist_name="st_nJets7"; 
	GetHistogram(hist_name,"RPV300"); 
	GetHistogram(hist_name,"RPV500"); 
	GetHistogram(hist_name,"RPV1000"); 
	
	CreateCanvas("signal_st","",700,600);
	CName["signal_st"]->cd();
	TLegend L2(0.55,0.55,0.85,0.85);
	L2.SetFillColor(10);
	L2.SetLineColor(10);
	L2.SetLineWidth(0); 
	
	hist_name="st_nJets4"; 
	hName[hist_name+"RPV300"]->SetLineColor(kRed);
	hName[hist_name+"RPV300"]->SetFillStyle(21); 
	hName[hist_name+"RPV500"]->SetLineColor(kBlue);
	hName[hist_name+"RPV500"]->SetFillStyle(22); 
	hName[hist_name+"RPV1000"]->SetLineColor(kBlack);
	hName[hist_name+"RPV1000"]->SetLineWidth(2);
	hName[hist_name+"RPV1000"]->SetFillStyle(23); 
	hName[hist_name+"RPV1000"]->Scale(500); 
	
	
	L2.AddEntry(hName[hist_name+"RPV300"],"m_{#tilde{t}=300 GeV}, 4-jets", "LEP");
	L2.AddEntry(hName[hist_name+"RPV500"],"m_{#tilde{t}=500 GeV}, 4-jets", "LEP");
	L2.AddEntry(hName[hist_name+"RPV1000"],"500 #times m_{#tilde{t}=1000 GeV}, 4-jets", "LEP");

	hName[hist_name+"RPV300"]->Draw("hist");
	hName[hist_name+"RPV500"]->Draw("hist same");
	hName[hist_name+"RPV1000"]->Draw("hist same");
	
	hist_name="st_nJets7"; 
	hName[hist_name+"RPV300"]->SetLineColor(kRed);
	hName[hist_name+"RPV500"]->SetLineColor(kBlue);
	hName[hist_name+"RPV1000"]->SetLineColor(kBlack);
	hName[hist_name+"RPV1000"]->SetLineWidth(2);
	hName[hist_name+"RPV300"]->SetLineStyle(kDashed);
	hName[hist_name+"RPV500"]->SetLineStyle(kDashed);
	hName[hist_name+"RPV1000"]->SetLineStyle(kDashed);
	

	L2.AddEntry(hName[hist_name+"RPV300"],"#times m_{#tilde{t}=300 GeV}, 7-jets", "LEP");
	L2.AddEntry(hName[hist_name+"RPV500"],"10 #times m_{#tilde{t}=500 GeV}, 7-jets", "LEP");
	L2.AddEntry(hName[hist_name+"RPV1000"],"100 #times m_{#tilde{t}=1000 GeV}, 7-jets", "LEP");
	
//	hName[hist_name+"RPV300"]->Scale(10); 
	hName[hist_name+"RPV500"]->Scale(10); 
	hName[hist_name+"RPV1000"]->Scale(100); 

	hName[hist_name+"RPV300"]->Draw("hist same");
	hName[hist_name+"RPV500"]->Draw("hist same");
	hName[hist_name+"RPV1000"]->Draw("hist same");
	L2.DrawClone(); 

	
}

void plot_npv(){

	CreateCanvas("npv","",600,600); 
	
	CName["npv"]->cd();
	TString hist_name="h_npv"; 
	
	TLegend L(0.6,0.5,0.8,0.8); 
	
	GetHistogram(hist_name,"singleMu");
	GetHistogram(hist_name,"allMC");
	GetHistogram(hist_name,"RPV500");
	GetHistogram(hist_name,"stealth_ww500_300");
	
	marker_style(hName[hist_name+"stealth_ww500_300"]); 
	hName[hist_name+"stealth_ww500_300"]->SetMarkerColor(kGreen); 
	hName[hist_name+"stealth_ww500_300"]->DrawNormalized(); 
	
	marker_style(hName[hist_name+"singleMu"]); 
	hName[hist_name+"singleMu"]->DrawNormalized("same"); 	


	marker_style(hName[hist_name+"allMC"]);
	hName[hist_name+"allMC"]->SetMarkerColor(kRed);
	hName[hist_name+"allMC"]->DrawNormalized("same");
	
	marker_style(hName[hist_name+"RPV500"]); 
	hName[hist_name+"RPV500"]->SetMarkerColor(kBlue); 
	hName[hist_name+"RPV500"]->DrawNormalized("same"); 



	L.AddEntry(hName[hist_name+"singleMu"],"Data","LEP");
	L.AddEntry(hName[hist_name+"allMC"],"MC Bkg","LEP");
	L.AddEntry(hName[hist_name+"RPV500"],"RPV Signal","LEP");
	L.AddEntry(hName[hist_name+"stealth_ww500_300"],"Stleath Signal","LEP");

	L.DrawClone("same"); 
}


void fill_stack_inclusive(TString variable, bool blind){
		//add blinding/unblinding control 
	TLegend *L = new TLegend(0.65,0.45,0.85,0.85); 
	
	TString jet_canvas=variable+"_stack_comp";
	TString jet_stack=variable+"_hstack";
	
	CreateCanvas(jet_canvas,"stacked comparison",600,600); 
	CreateStack(jet_stack,"stack"); 

	TString file="QCD";
	TString name=variable; 
	GetHistogram(name,file);
	hName[name+file]->SetFillColor(7); 
	hName[name+file]->SetMarkerColor(7); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	L->AddEntry(hName[name+file],"QCD");
	
	stackName[jet_stack]->Add(hName[name+file]); 

	TH1F *h=(TH1F*)hName[name+file]->Clone("MC_stack"); 
	hName[h->GetName()]=h;
	
	
	file="diboson";
	GetHistogram(name,file);
	hName[name+file]->SetFillColor(9); 
	hName[name+file]->SetMarkerColor(2); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	
	L->AddEntry(hName[name+file],"Diboson"); 	
	stackName[jet_stack]->Add(hName[name+file]);
	hName["MC_stack"]->Add(hName[name+file]);
	
	file="singleTop";
	GetHistogram(name,file);
	hName[name+file]->SetFillColor(5); 
	hName[name+file]->SetMarkerColor(2); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	
	L->AddEntry(hName[name+file],"Single Top"); 	
	stackName[jet_stack]->Add(hName[name+file]);
	hName["MC_stack"]->Add(hName[name+file]);
	
	
	file="dy";
	GetHistogram(name,file);
	hName[name+file]->SetFillColor(4); 
	hName[name+file]->SetMarkerColor(2); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	
	L->AddEntry(hName[name+file],"Drell-Yan"); 	
	stackName[jet_stack]->Add(hName[name+file]);
	hName["MC_stack"]->Add(hName[name+file]);
	
	file="ttbar";
	GetHistogram(name,file);
	hName[name+file]->SetFillColor(3); 
	hName[name+file]->SetMarkerColor(3); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	L->AddEntry(hName[name+file],"TTbar"); 
	
	stackName[jet_stack]->Add(hName[name+file]); 
	hName["MC_stack"]->Add(hName[name+file]);
	
	file="wJets";
	GetHistogram(name,file);
	hName[name+file]->SetFillColor(2); 
	hName[name+file]->SetMarkerColor(4); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	//hName[name+file]->SetAxisRange(0,st_max);
	L->AddEntry(hName[name+file],"W+Jets"); 
	
	stackName[jet_stack]->Add(hName[name+file]); 
	hName["MC_stack"]->Add(hName[name+file]);
	 	
	file="RPV500";
	GetHistogram(name,file);
	hName[name+file]->SetMarkerStyle(20);
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetLineColor(kRed);
	L->AddEntry(hName[name+file],"m_{#tilde{t}}=500 GeV"); 
	
	file="RPV300";
	GetHistogram(name,file);
	hName[name+file]->SetMarkerStyle(20);
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetLineColor(kRed);
	hName[name+file]->SetLineStyle(kDashed);
	L->AddEntry(hName[name+file],"m_{#tilde{t}}=300 GeV"); 
	
	file="stealth_ww500_300";
	GetHistogram(name,file);
	hName[name+file]->SetLineStyle(2);
	hName[name+file]->SetLineWidth(2);
	hName[name+file]->SetLineColor(kOrange);
	L->AddEntry(hName[name+file],"stealth, M_{squark}=500 GeV"); 
	
	file="singleMu";
	GetHistogram(name,file);
	hName[name+file]->SetMarkerStyle(20);
	hName[name+file]->SetMarkerSize(0.5);
	if(blind) hName[name+file]->SetAxisRange(1,5);
	L->AddEntry(hName[name+file],"data"); 
	

	TPad *shape_pad = new TPad("stpad","st-shape",0.0,0.3,1,1);
	TPad *ratio_pad = new TPad("ratiopad","Ratio",0,0.,1,0.3);
	
	
	CName[jet_canvas]->cd(); 
	shape_pad->Draw();
	ratio_pad->Draw();
	
	shape_pad->cd();
	gPad->SetLogy();
	
	float min=hName[name+file]->GetBinContent(hName[name+file]->FindBin(2750));
	if(min<0.001) min=0.001; 
//	stackName[jet_stack]->SetMinimum(min);
	stackName[jet_stack]->Draw("hist");
	hName[name+file]->DrawCopy("same");
	hName[name+file]->DrawCopy("same");
	
	hName[name+"RPV500"]->DrawCopy("same");
	hName[name+"RPV300"]->DrawCopy("same");
	hName[name+"stealth_ww500_300"]->DrawCopy("same");

	L->Draw("same"); 

	ratio_pad->cd(); 
	cout << "divide: " << endl; 
	hName[name+file]->Divide(hName["MC_stack"]); 
	hName[name+file]->GetYaxis()->SetTitle("Data/MC"); 
	hName[name+file]->SetMinimum(0.5);
	hName[name+file]->SetMaximum(1.5);
	
	hName[name+file]->Draw(); 

	CName[jet_canvas]->Update(); 
	
	
}


void fill_stack(TString variable, bool blind){
	//modify this code to be given a string with the name
	//create just one stack per call
	//fix y axis min thing
	//plot h_mumu in bins
	bool print=true; 
	TLegend *L = new TLegend(0.6,0.5,0.8,0.7); 

	if(print) cout << "variable: " << variable << endl; 
	
		TString jet_canvas=variable+"_stack_comp"; 
		TString jet_stack=variable+"_hstack"; 
		CreateCanvas(jet_canvas,"stacked comparison",600,600); 
		CreateStack(jet_stack,variable);  

		if(print)cout << "variable: " << variable << endl; 
		float st_max=3000;
				
		if(print) cout << "QCD" << endl; 
		 TString file="QCD";
		 TString name;
		 name=variable;
		
		 GetHistogram(name,file);
		 hName[name+file]->SetAxisRange(0,st_max);
		 hName[name+file]->SetFillColor(7); 
		 hName[name+file]->SetMarkerColor(7); 
		 hName[name+file]->SetMarkerStyle(21); 
		 hName[name+file]->SetMarkerSize(0.5);
		 L->AddEntry(hName[name+file],"QCD");
		 
		 stackName[jet_stack]->Add(hName[name+file]); 
		TString tot_MC="MC_stack_"+variable; 

		TH1F *h=(TH1F*)hName[name+file]->Clone(tot_MC);
		hName[h->GetName()]=h;
		 
		if(print)cout << "diboson" << endl; 
		file="diboson";
		//TString name=variable+Form("_nJets%d_GeV",nJets);
		
		GetHistogram(name,file);
		hName[name+file]->SetFillColor(9); 
		hName[name+file]->SetMarkerColor(2); 
		hName[name+file]->SetMarkerStyle(21); 
		hName[name+file]->SetMarkerSize(0.5);
		hName[name+file]->SetAxisRange(0,st_max);
		L->AddEntry(hName[name+file],"Diboson"); 
		
		stackName[jet_stack]->Add(hName[name+file]);
		hName[tot_MC]->Add(hName[name+file]);
		
		if(print) cout <<"singleTop" << endl;
		file="singleTop";
		//TString name=variable+Form("_nJets%d_GeV",nJets);
		
		GetHistogram(name,file);
		hName[name+file]->SetFillColor(5); 
		hName[name+file]->SetMarkerColor(2); 
		hName[name+file]->SetMarkerStyle(21); 
		hName[name+file]->SetMarkerSize(0.5);
		hName[name+file]->SetAxisRange(0,st_max);
		L->AddEntry(hName[name+file],"Single Top"); 
		
		stackName[jet_stack]->Add(hName[name+file]);
		hName[tot_MC]->Add(hName[name+file]);

		
		if(print) cout << "dy"  << endl; 
		file="dy";
		//TString name=variable+Form("_nJets%d_GeV",nJets);

		GetHistogram(name,file);
		hName[name+file]->SetFillColor(4); 
		hName[name+file]->SetMarkerColor(2); 
		hName[name+file]->SetMarkerStyle(21); 
		hName[name+file]->SetMarkerSize(0.5);
		hName[name+file]->SetAxisRange(0,st_max);
		L->AddEntry(hName[name+file],"Drell-Yan"); 

		stackName[jet_stack]->Add(hName[name+file]);
		hName[tot_MC]->Add(hName[name+file]);
		
		if(print) cout << "ttbar" << endl; 
		file="ttbar";
		GetHistogram(name,file);
		hName[name+file]->SetFillColor(3); 
		hName[name+file]->SetMarkerColor(3); 
		hName[name+file]->SetMarkerStyle(21); 
		hName[name+file]->SetMarkerSize(0.5);
		hName[name+file]->SetAxisRange(0,st_max);
		L->AddEntry(hName[name+file],"TTbar"); 

		stackName[jet_stack]->Add(hName[name+file]); 
		hName[tot_MC]->Add(hName[name+file]);
		if(print) cout << "wJets" << endl; 
		file="wJets";
		GetHistogram(name,file);
		hName[name+file]->SetFillColor(2); 
		hName[name+file]->SetMarkerColor(4); 
		hName[name+file]->SetMarkerStyle(21); 
		hName[name+file]->SetMarkerSize(0.5);
		hName[name+file]->SetAxisRange(0,st_max);
		L->AddEntry(hName[name+file],"W+Jets"); 

		stackName[jet_stack]->Add(hName[name+file]); 
		hName[tot_MC]->Add(hName[name+file]);
		
		//stackName[jet_stack]->Add(hName[name+file]); 
		
		if(print) cout << "RPV500" << endl; 
			file="RPV500";
			GetHistogram(name,file);
			hName[name+file]->SetLineStyle(1);
			hName[name+file]->SetLineWidth(2);
			hName[name+file]->SetLineColor(kRed);
			hName[name+file]->SetAxisRange(0,st_max);
			L->AddEntry(hName[name+file],"RPV: M_{#tilde{t}}=500 GeV"); 
		if(print) cout << "stealth500/300" << endl; 
		file="stealth_ww500_300";
			GetHistogram(name,file);
			hName[name+file]->SetLineStyle(2);
			hName[name+file]->SetLineWidth(2);
			hName[name+file]->SetLineColor(kRed);
			hName[name+file]->SetAxisRange(0,st_max);
			L->AddEntry(hName[name+file],"M_{squark}=500 GeV"); 
		if(print) cout << "stelath 1400/300" << endl; 
			file="stealth_ww1400_300";
			GetHistogram(name,file);
			hName[name+file]->SetLineStyle(1);
			hName[name+file]->SetLineWidth(2);
			hName[name+file]->SetLineColor(kOrange);
			hName[name+file]->SetAxisRange(0,st_max);
			L->AddEntry(hName[name+file],"M_{squark}=1400 GeV #times 10,000"); 
		
		if(print) cout << "single mu: " << endl; 
		file="singleMu";
		GetHistogram(name,file);
		hName[name+file]->SetMarkerStyle(20);
		hName[name+file]->SetMarkerSize(0.5);
		L->AddEntry(hName[name+file],"data"); 
	
	TString shapepad_name="stpad_"+variable;
	TString ratiopad_name="ratiopad_"+variable; 
		
		TPad *shape_pad = new TPad(shapepad_name,"st-shape",0.0,0.3,1,1);
		TPad *ratio_pad = new TPad(ratiopad_name,"Ratio",0,0.,1,0.3);
	
		CName[jet_canvas]->cd(); 
		shape_pad->Draw();
		ratio_pad->Draw();

		shape_pad->cd();
		gPad->SetLogy();
		//if(nJets<=5){
			hName[name+file]->GetXaxis()->SetTitle("S_{T} [GeV]"); 
			hName[name+file]->GetYaxis()->SetTitle("Events/GeV"); 
			hName[name+file]->DrawCopy();
			stackName[jet_stack]->Draw("hist same");
			if(!blind)hName[name+file]->DrawCopy("same");

			float min=hName[name+file]->GetBinContent(hName[name+file]->FindBin(2750));
			
			float max1=hName[name+file]->GetMaximum(); 
			float max2=stackName[jet_stack]->GetMaximum(); 
			float max=TMath::Max(max1,max2); 
			/*
			if(min<0.0001) min=0.0001;
			stackName[jet_stack]->SetMinimum(min);
			stackName[jet_stack]->SetMaximum(max); 
			stackName[jet_stack]->Draw("hist");
			stackName[jet_stack]->GetXaxis()->SetTitle("S_{T} [GeV]");
			stackName[jet_stack]->GetYaxis()->SetTitle("Events/GeV"); 
			stackName[jet_stack]->Draw("hist");
			hName[name+"RPV500"]->DrawCopy("same");
			*/
			hName[name+"stealth_ww1400_300"]->Scale(10000);
			hName[name+"RPV500"]->DrawCopy("same"); 
			hName[name+"stealth_ww1400_300"]->DrawCopy("same"); 
			hName[name+"stealth_ww500_300"]->DrawCopy("same"); 

		
		L->DrawClone("same"); 
		
	
		ratio_pad->cd(); 
	
		TLatex BL(0.2,0.5,"Blinded!!!"); 
		BL.SetNDC(kTRUE); 
		BL.SetTextSize(0.1); 
		hName[name+file]->Divide(hName[tot_MC]); 
		hName[name+file]->GetYaxis()->SetTitle("Data/MC"); 
		hName[name+file]->SetMinimum(0.5);
		hName[name+file]->SetMaximum(1.5);
		hName[name+file]->DrawCopy();
		if(!blind)hName[name+file]->DrawCopy(); 
		//else BL.DrawClone(); 
	
}

void GetHistogram(TString hist_name, TString file){
	//cout <<" filename: "<< file <<  " hist_name: "<< hist_name << endl; 
	TH1F *h=(TH1F*)fName[file]->FindObjectAny(hist_name);
	h->SetStats(kFALSE);
	hName[hist_name+file]=h; 
}
