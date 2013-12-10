/*
 *  susy_plots.c
 *  
 *
 *  Created by Benjamin Carlson on 11/15/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "susy_plots.h"

void susy_plots(){
	gROOT->SetBatch();
	open_files(); // opens all files and loads all histograms. 
	combine_histograms(); // combines several histograms 
	compute_ratios(); 
	
	compare_ratios(3,2,1);
	compare_ratios(3,2,2);

	compare_ratios(4,2,1);
	compare_ratios(4,2,2);
	
	compare_ratios(4,3,1);
	compare_ratios(4,3,2);
	
	compare_ratios(5,3,1);
	compare_ratios(5,3,2);
	
	plot_ratios(2,1);
	plot_ratios(3,1);
	plot_ratios(2,2);
	plot_ratios(3,2);
	
	compute_correction(2,300,600);
	compute_correction(3,300,600);
	compute_correction(4,300,600);
	
	compute_correction(2,600,1000);
	compute_correction(3,600,1000);
	compute_correction(4,600,1000);
	
	compute_correction(2,1000,-1);
	compute_correction(3,1000,-1);
	compute_correction(4,1000,-1);
	print_numbers();
	
	corrections_table();
	 
	for(int nJ=1;nJ<=7; nJ++){
		bool blind=false; 
		if(nJ>5) blind=true; 
		for(int mu=1; mu<=2; mu++){	
			for (std::map<TString,TFile*>::iterator it=fName.begin(); it!=fName.end(); it++) {
				if(it->first!="_singleMu")frac_bkg(Form("st_nJets%d_%dMu_GeV",nJ,mu),it->first); 
			}
		}
		
		draw_frac(Form("st_nJets%d_1Mu_GeV",nJ),Form("bkg_frac_nJets%d_st1mu",nJ)); 
		draw_frac(Form("st_nJets%d_2Mu_GeV",nJ),Form("bkg_frac_nJets%d_st2mu",nJ)); 

		fill_stack(Form("st_nJets%d_1Mu_GeV",nJ),blind); 
		fill_stack(Form("st_nJets%d_2Mu_GeV",nJ),blind); 
		fill_stack(Form("st_nJets%d_1Mu_1El",nJ),blind);
		
		fill_stack(Form("st_nJets%d_W",nJ),blind); 
		fill_stack(Form("st_nJets%d_Z",nJ),blind); 

		
		fill_stack(Form("h_mumu_nJet%d",nJ),0); 
		fill_stack(Form("mt_nJets%d_1Mu",nJ),0);		
		fill_stack(Form("mt_nJets%d_2Mu",nJ),0);
		fill_stack(Form("met_nJets%d_1Mu",nJ),0);
		fill_stack(Form("met_nJets%d_2Mu",nJ),0);
		fill_stack(Form("DeltaR_Jet%d_%dMu",nJ,1),0); 
		fill_stack(Form("DeltaPhi_Jet%d_%dMu",nJ,1),0); 

		fill_stack(Form("DeltaR_Jet%d_%dMu",nJ,2),0); 
		
		fill_stack(Form("DeltaPhiJetMet_Jet%d_%dMu",nJ,1),0);
		fill_stack(Form("DeltaPhiJetMet_Jet%d_%dMu",nJ,2),0);
		
		fill_stack(Form("nJets%d_DeltaPhiMuonMet_1Mu",nJ),0);
		fill_stack(Form("nJets%d_DeltaPhiMuonMet_2Mu",nJ),0); 
		
		fill_stack(Form("nJets%d_muonpT1",nJ),0);
		fill_stack(Form("nJets%d_muonpT2",nJ),0); 
	}
	

	fill_stack("jetpT1_1Mu",0);
	fill_stack("jetpT2_1Mu",0);
	fill_stack("jetpT3_1Mu",0);
	fill_stack("jetpT4_1Mu",0);
	
	fill_stack("jetpT1_2Mu",0);
	fill_stack("jetpT2_2Mu",0);
	fill_stack("jetpT3_2Mu",0);
	fill_stack("jetpT4_2Mu",0);
	
	fill_stack("jet1_deltaeta_RF",0);
	fill_stack("jet1_deltaphi_RF",0);
	
	fill_stack("jet1_muMass",0); 
	fill_stack("jet1_muMass_pt1",0); 
	fill_stack("jet1_muMass_pt2",0); 
	fill_stack("jet1_muMass_pt3",0); 

	fill_stack("jet1_munuMass",0); 
	fill_stack("jet1_nuMass",0); 
	fill_stack("jet1_mujetptM1",0); 
	fill_stack("jet1_mujetptM2",0); 
	fill_stack("jet1_mujetptM3",0); 



	fill_stack("jet1_ptmu35_muMass",0); 
	fill_stack("jet1_ptmu40_muMass",0); 
	fill_stack("jet1_ptmu45_muMass",0); 

	fill_stack("jet1_ptjet35_muMass",0); 
	fill_stack("jet1_ptjet40_muMass",0); 
	fill_stack("jet1_ptjet45_muMass",0); 
	
	fill_stack("bjet1_muMass",0); 

	fill_stack("jet2_muMass",0); 
	fill_stack("bjet2_muMass",0); 

	fill_stack("nJets_W",0);
	fill_stack("nJets_Z",0);
	fill_stack("nJets_1Mu_1El",0);
	
	fill_stack("nJets_1Mu",0);		
	fill_stack("nJets_2Mu",0);	
	
	for(int bkg_nJ=2; bkg_nJ<=4; bkg_nJ++){
		for(int mu=1; mu<=2; mu++){
			fractional_uncertainty(bkg_nJ,mu,"_singleMu");
		}
	}
	for(int nJ=4; nJ<=7; nJ++){
		for(int mu=1; mu<=2; mu++){
			for(int bkg_nJ=2; bkg_nJ<=4; bkg_nJ++){
				if(nJ<=bkg_nJ)continue; 
				bkg_est(nJ,bkg_nJ,mu,"_singleMu");
				plot_expected(nJ,bkg_nJ,mu,"_singleMu"); 
			}
			overlay_bkg_expectations(nJ,mu,"_singleMu");
		}//mu
	}//nJ

	
		for(int nJ=4; nJ<=7; nJ++){
			for(int mu=1; mu<=2; mu++){
				CreateHistogram(acc_name(nJ,mu),"","M_{#tilde{t}} [GeV]","acc #times #epsilon",8,250,1050); 
				for(int irpv=300; irpv<=1000;irpv+=100){
				TString signal=Form("_RPV%d",irpv);
				compute_sensitivity(nJ,mu,signal); 
				compute_acceptance(nJ,mu,irpv); 
			}//irpv loop
		}//nJ loop
	}//irpv loop 
	
	
	
	write();
}

void write(){
	TFile *output_file = new TFile("susy_plots.root","RECREATE"); 
	output_file->cd();
	output_file->mkdir("histograms"); 
	output_file->cd("histograms"); 
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		//if(it->first.Contains("23"))cout << it->first << endl; 
		//cout << hName[it->first]->GetName() << endl; 
		if(it->first.Contains("st") && it->first.Contains("Mu"))hName[it->first]->Write();
		if(it->first.Contains("acc") || it->first.Contains("sb"))hName[it->first]->Write();

	}
	output_file->mkdir("graphs"); 
	output_file->cd("graphs"); 
	for (std::map<TString,TGraphAsymmErrors*>::iterator it=grName.begin(); it!=grName.end(); it++) {
		grName[it->first]->Write();
	}
	output_file->cd(); 
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		TString name="/uscms_data/d3/btcarlso/Figures_SUSYAN/"; 
		name=name+it->first+".png"; 
		CName[it->first]->Print(name);
	}
}

void plot_ratios(int den, int mu){
	CreateCanvas(Form("ratios_%dden_%dMu",den,mu),"",Cx,Cy); 
	CName[Form("ratios_%dden_%dMu",den,mu)]->cd(); 

	TLegend L(0.19,0.72,0.35,0.88);
	L.SetFillColor(10);
	L.SetLineColor(10);
	L.SetBorderSize(0);
	
	L.AddEntry(grName[ratio_name(4,den,mu,"_allMC")],grName[ratio_name(4,den,mu,"_allMC")]->GetYaxis()->GetTitle(),"LEP"); 

	grName[ratio_name(4,den,mu,"_allMC")]->Draw("ap"); 
	grName[ratio_name(4,den,mu,"_allMC")]->GetYaxis()->SetTitle(Form("#alpha(n-jets,%d-jets)",den)); 
	grName[ratio_name(4,den,mu,"_allMC")]->GetXaxis()->SetTitle("S_{T} [GeV]");

	grName[ratio_name(4,den,mu,"_allMC")]->Draw("ap"); 
	
	for(int num=5; num<=7; num++){
		grName[ratio_name(num,den,mu,"_allMC")]->Draw("p same"); 
		L.AddEntry(grName[ratio_name(num,den,mu,"_allMC")],grName[ratio_name(num,den,mu,"_allMC")]->GetYaxis()->GetTitle(),"LEP"); 
	}
	draw_headersim();
	L.DrawClone(); 
	
}

void open_files(){
	//Open MC files
	
	open_file("QCD");
	for(int irpv=300; irpv<=1000; irpv+=100){
		open_file(Form("RPV%d",irpv));
	}
//	open_file("stealth_ww500_300");
//	open_file("stealth_ww1400_300");
	open_file("wJets");
	open_file("dy");
	open_file("ttbar");
	open_file("singleTop");
	open_file("diboson"); 
	open_file("singleMu");
	open_file("allMC"); 
	
	Ngen["_RPV300"]=39198;
	Ngen["_RPV400"]=39999;
	Ngen["_RPV500"]=34397;
	Ngen["_RPV600"]=39193;
	Ngen["_RPV700"]=39998;
	Ngen["_RPV800"]=39998;
	Ngen["_RPV900"]=39997;
	Ngen["_RPV1000"]=39996;

	
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

void compare_ratios(int num, int den,int mu){
	cout << "Ratio comparison: " << num << "-" << den << endl; 
	TString canvas_name=Form("Compare_R%d-%d_%dMu",num,den,mu);
	if(mu==-1)canvas_name=Form("Compare_R%d-%d",num,den);
	
	CreateCanvas(canvas_name,canvas_name,600,600);
	CName[canvas_name]->cd();
	
	TString name=Form("Ratio%d-%d_%dMu",num,den,mu);
	if(mu==-1)name=Form("Ratio%d-%d",num,den);
	TString namedata=ratio_name(num,den,mu,"_singleMu"); 
	TString nameMC=ratio_name(num,den,mu,"_allMC"); 
	
	TLatex txt(0.2,0.8,Form("#mu=%d",mu)); 
	txt.SetNDC(kTRUE);
	
	
	TGraphAsymmErrors *grMC = (TGraphAsymmErrors*)grName[nameMC]->Clone(nameMC+"_norm"); 
	TGraphAsymmErrors *grdata = (TGraphAsymmErrors*)grName[namedata]->Clone(namedata+"_norm"); 
	//double SMC=1./avg(grMC); 
	//double Sdata=1./avg(grdata); 
	
	//scale_graph(grMC,SMC); 
	//scale_graph(grdata,Sdata); 
	
	setcolor(grMC,6); 
	setcolor(grdata,4);
	
	TLegend *L = new TLegend(0.45,0.15,0.85,0.35);
	L->SetFillColor(10); 
	L->SetLineColor(10);
	L->SetLineWidth(0);
	L->AddEntry(grdata,"data", "LEP"); 
	L->AddEntry(grMC,"Combined MC", "LEP"); 
	
	grdata->Draw("ap "); 
	grdata->GetXaxis()->SetRangeUser(300,3000); 
	grdata->GetYaxis()->SetRangeUser(0,2*grdata->Eval(1200)); 
	grdata->GetYaxis()->SetTitle(Form("#alpha(%djets/%djets)",num,den)); 
	grdata->DrawClone("ap"); 
	grMC->DrawClone("p same");
	txt.DrawClone();
	L->DrawClone("same");
	draw_header();
	
}

TString ratio_name(int num, int den, int mu,TString file){
	TString name=Form("Ratio%d-%d_%dMu",num,den,mu); 

	name=name+file; 
	return name; 
}

TString expected_bkg(int nJ, int bkg_nJ, int mu, TString mode, TString file){
	TString name=Form("stexp_nJets%d_bkg%d_%dMu"); 
	name=name+mode+file; 
	return name; 
}

TString bkg_systematic(int nJ, int bkg_nJ, int mu, TString mode, TString file){
	TString name=Form("stsys_nJets%d_bkg%d_%dMu",nJ,bkg_nJ,mu); 
	name=name+mode+file; 
	
	if(mode=="stat"){
		name=Form("ststat_nbkg%d_%dMu",bkg_nJ,mu);
		name=name+file;
	}
	
	return name; 
}

void overlay_bkg_expectations(int nJ, int mu, TString file){
	TString canvas_name=Form("Overlay_Plot_bkgexp_nJets%d_%dMu",nJ,mu);
	canvas_name=canvas_name+file; 
	CreateCanvas(canvas_name,"",Cx,Cy); 
	gPad->SetLogy(); 
	int bkg_nJ=2; 
	TLegend L(0.65,0.65,0.85,0.85);
	L.SetFillColor(10);
	L.SetLineColor(10);
	L.SetBorderSize(0);

	
	TLatex txt(0.2,0.8,Form("#mu=%d",mu));
	txt.SetNDC(kTRUE);
	
	L.AddEntry(hName[expected_bkg(nJ,bkg_nJ,mu,"",file)],Form("%d-jets, %d-jets bkg",nJ,bkg_nJ),"L"); 
	hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->DrawCopy(); 
	
	bkg_nJ=3; 
	hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->SetLineColor(kBlue); 
	L.AddEntry(hName[expected_bkg(nJ,bkg_nJ,mu,"",file)],Form("%d-jets, %d-jets bkg",nJ,bkg_nJ),"L"); 
	hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->DrawCopy("same"); 
	if(nJ>4){
		bkg_nJ=4;
		hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->SetLineColor(kBlack); 
		L.AddEntry(hName[expected_bkg(nJ,bkg_nJ,mu,"",file)],Form("%d-jets, %d-jets bkg",nJ,bkg_nJ),"L"); 
		hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->DrawCopy("same"); 
	}
	txt.DrawClone();
	L.DrawClone(); 
	draw_header();		   
}

void plot_expected(int nJ, int bkg_nJ, int mu, TString file){
	cout << "plot expected: " << endl; 
	TString canvas_name=Form("Plot_bkgexp_nJets%d_%dbkg_%dMu",nJ,bkg_nJ,mu);
	canvas_name=canvas_name+file; 
	CreateCanvas(canvas_name,"",Cx,Cy); 
	
	TPad *shape_pad = new TPad(Form("stpad_%d",nJ),"st-shape",0.0,0.3,1,1);
	TPad *ratio_pad = new TPad(Form("ratiopad_%d",nJ),"Ratio",0,0.,1,0.3);
	CName[canvas_name]->cd();
	
	shape_pad->Draw();
	ratio_pad->Draw();
	
	shape_pad->cd(); 
	gPad->SetLogy();
	cout << "compute ratios: "<< endl; 
	
	TString obs_exp=Form("Obs_exp_nJets%d_%dMu",nJ,mu); 
	obs_exp=obs_exp+file; 
	TH1F *ratio=(TH1F*)hName[st_name(nJ,mu,"",file)]->Clone(obs_exp); 
	
	TString data_MC=Form("data_MC_nJets%d_%dMu",nJ,mu); 
	data_MC=data_MC+file; 
	TH1F *ratioD_MC=(TH1F*)hName[st_name(nJ,mu,"",file)]->Clone(data_MC); 
	
	ratio->Divide(hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]); 
	ratio->SetAxisRange(500,3000);
	
	ratioD_MC->Divide(hName[st_name(nJ,mu,"","_allMC")]); 
	ratioD_MC->SetAxisRange(500,3000);
	cout << "data/MC ratios done. " << endl; 
	double chi2_MC=0;
	double chi2_datadriven=0; 
	
	for(int i=1; i<=ratio->GetNbinsX();i++){
		chi2_MC+=TMath::Power(ratioD_MC->GetBinContent(i),2); 
		chi2_datadriven+=TMath::Power(ratio->GetBinContent(i),2); 
	}
	
	TGraphAsymmErrors *grRatio=new TGraphAsymmErrors(ratio); 
	TGraphAsymmErrors *grRatio_stat = new TGraphAsymmErrors(ratio); 
	
	for(int i=1; i<=ratio->GetNbinsX(); i++){
		double st=hName[st_name(nJ,mu,"",file)]->GetBinCenter(i);
		double sigmaEp=hName[bkg_systematic(nJ,bkg_nJ,mu,"Ep",file)]->GetBinContent(i); 
		double sigmaEm=hName[bkg_systematic(nJ,bkg_nJ,mu,"Em",file)]->GetBinContent(i); 
		double stat=hName[bkg_systematic(nJ,bkg_nJ,mu,"stat",file)]->GetBinContent(i); 

	//	cout << "sigmaEp: " << sigmaEp << " stat: " << stat	 << endl; 
		
		double R=ratio->GetBinContent(i);
		grRatio->SetPoint(i-1,st,R); 
		grRatio_stat->SetPoint(i-1,st,R); 
		grRatio->SetPointError(i-1,0,0,R*sigmaEm,R*sigmaEp); 
		grRatio_stat->SetPointError(i-1,0,0,R*stat,R*stat); 
		
	}
	TLegend L(0.6,0.65,0.8,0.85); 
	L.SetFillColor(10);
	L.SetLineColor(10);
	L.SetBorderSize(0);
	
	L.AddEntry(hName[st_name(nJ,mu,"","_allMC")],"MC prediction");
	L.AddEntry(hName[st_name(nJ,mu,"","_singleMu")],"Data","LEP");
	L.AddEntry(hName[expected_bkg(nJ,bkg_nJ,mu,"",file)],"Data-driven prediction","L");
	L.AddEntry(hName[expected_bkg(nJ,bkg_nJ,mu,"Ep",file)],"Data-driven prediction systematic","L");

	
	hName[st_name(nJ,mu,"","_allMC")]->SetFillColor(kBlue); 
	hName[st_name(nJ,mu,"","_allMC")]->SetLineColor(kBlue); 
	hName[st_name(nJ,mu,"",file)]->SetMarkerStyle(22); 
	hName[st_name(nJ,mu,"",file)]->SetLineColor(kBlack); 
//	hName[st_name(nJ,mu,"","_allMC")]->SetTitle(""); 
	
	hName[st_name(nJ,mu,"","_allMC")]->SetAxisRange(500,3000);
	hName[st_name(nJ,mu,"","_allMC")]->DrawCopy("histo"); 
	hName[st_name(nJ,mu,"",file)]->DrawCopy("E1 same"); 
	hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->DrawCopy("same");
	hName[expected_bkg(nJ,bkg_nJ,mu,"Ep",file)]->DrawCopy("same");
	hName[expected_bkg(nJ,bkg_nJ,mu,"Em",file)]->DrawCopy("same");
	L.DrawClone();
	draw_header();
	
	TLatex Lat(0.3,0.8,Form("%d-jets, %d-jets bkg, #mu=%d",nJ,bkg_nJ,mu)); 
	Lat.SetNDC(kTRUE);
	Lat.DrawClone(); 
	
	
	//h->DrawCopy("same");
	
	TF1 *line= new TF1("line","pol0",500,3000); 
	line->SetParameter(0,1); 
	line->SetLineColor(kRed); 
	ratio_pad->cd();


	grRatio->SetMarkerStyle(7);
	grRatio->SetMarkerSize(1); 
	grRatio_stat->SetLineColor(kRed); 
	grRatio->Draw("ap"); 
	grRatio->GetXaxis()->SetRangeUser(500,3000); 
	grRatio->GetYaxis()->SetTitle("Data/Pred"); 
	grRatio->GetXaxis()->SetTitle("S_{T} [GeV]"); 
//	grRatio->SetTitle(""); 
	grRatio->DrawClone("ap"); 
	grRatio_stat->DrawClone("|| same"); 
	
	ratioD_MC->SetMarkerColor(kBlue); 
	ratioD_MC->SetLineColor(kBlue); 

	ratioD_MC->DrawCopy("histo same"); 
	line->DrawClone("same"); 



	
}

void fractional_uncertainty(int bkg_nJ, int mu, TString file){
	int nJ=0; 
	float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	
	CreateHistogram(bkg_systematic(nJ,bkg_nJ,mu,"stat",file),"","S_{T}","Fractional Uncertainty, statistical",N,st_bins);

	for(int i=1; i<=hName[st_name(bkg_nJ,mu,"",file)]->GetNbinsX(); i++){
		
		double NE=hName[st_name(bkg_nJ,mu,"",file)]->GetBinError(i); 
		double N=hName[st_name(bkg_nJ,mu,"",file)]->GetBinContent(i); 
		hName[bkg_systematic(nJ,bkg_nJ,mu,"stat",file)]->SetBinContent(i,NE/N); 
	}
}


TString sig_bkg(int nJ,int mu,TString file){
	TString name= Form("sb_nJets%d_%dMu",nJ,mu);
	name=name+file; 
	return name; 
}

void compute_sensitivity(int nJ, int mu, TString signal){
	//cout << "sensitivity (mu): " << mu << " nJ: " << nJ << " file: " <<signal << endl; 
	TString file = "_singleMu"; 
	int bkg_nJ=3; 
	float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	
	CreateHistogram(sig_bkg(nJ,mu,signal),"","S_{T}^{thresh} [GeV]","S/#sqrt(B)",N,st_bins); 
	
	
	int Nbins=hName[st_name(nJ,mu,"",signal)]->GetNbinsX(); 
	for(int ist=1; ist<=Nbins; ist++){
		
		double bkg=hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->Integral(ist,Nbins); 
		double sig=hName[st_name(nJ,mu,"",signal)]->Integral(ist,Nbins);
		
		if(isnan(sig/TMath::Sqrt(bkg))==1 || isinf(sig/TMath::Sqrt(bkg))==1) continue; 
		//cout << "ist: " << ist << sig/TMath::Sqrt(bkg) << endl; 
		hName[sig_bkg(nJ,mu,signal)]->SetBinContent(ist,sig/TMath::Sqrt(bkg)); 
		hName[sig_bkg(nJ,mu,signal)]->SetBinError(ist,0); 
	}
	
}

TString acc_name(int nJ, int mu){
	TString name=Form("acc_nJets%d_%dMu",nJ,mu);
	return name; 
}

void compute_acceptance(int nJ,int mu, int irpv){
	TString signal=Form("_RPV%d",irpv); 
	int ist_max=hName[sig_bkg(nJ,mu,signal)]->GetMaximumBin(); 
	int min_ist=hName[sig_bkg(nJ,mu,signal)]->FindBin(500);
	if(ist_max>min_ist)ist_max=min_ist; 
	int N=hName[sig_bkg(nJ,mu,signal)]->GetNbinsX(); 
	double acc=hName[st_name(nJ,mu,"_uW",signal)]->Integral(ist_max,N)/Ngen[signal]; 
	hName[acc_name(nJ,mu)]->Fill(irpv,acc); 
	hName[acc_name(nJ,mu)]->SetBinError(hName[acc_name(nJ,mu)]->FindBin(irpv),0); 
}


void bkg_est(int nJ, int bkg_nJ, int mu, TString file){
	//cout << "nJ: " << nJ << " bkg_nJ: " << bkg_nJ << " mu: "<< mu  << " " << file << endl; 
	//cout << "Get ratios and create exp. histograms: " << endl; 
	float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	TString RN=ratio_name(nJ,bkg_nJ,mu,"_allMC");//Ratio name, always get the transfer factor (ratio) from MC
	
	CreateHistogram(expected_bkg(nJ,bkg_nJ,mu,"",file),"","S_{T}","Expected Events",N,st_bins); 
	CreateHistogram(expected_bkg(nJ,bkg_nJ,mu,"Ep",file),"","S_{T}","Fractional Uncertainty, +#sigma",N,st_bins);
	CreateHistogram(expected_bkg(nJ,bkg_nJ,mu,"Em",file),"","S_{T}","Fractioanl Uncertainty, -#sigma",N,st_bins);

	CreateHistogram(bkg_systematic(nJ,bkg_nJ,mu,"Ep",file),"","S_{T}","Fractional Uncertainty, +#sigma",N,st_bins);
	CreateHistogram(bkg_systematic(nJ,bkg_nJ,mu,"Em",file),"","S_{T}","Fractioanl Uncertainty, -#sigma",N,st_bins);


	//	double st_threshold[]={500,500,875,1100,1300,1300,2000};
	
	//cout << "Fill expected histogram: " << endl; 
	
	for(int ist=1; ist<=hName[st_name(nJ,mu,"",file)]->GetNbinsX(); ist++){
		double st=hName[st_name(nJ,mu,"",file)]->GetBinCenter(ist);
		double N=hName[st_name(bkg_nJ,mu,"",file)]->GetBinContent(ist); 
		int iB=hName[st_name(nJ,mu,"",file)]->FindBin(st)-1; 
		double Ncor=grName[RN]->GetY()[iB]*N; 
		//cout << "iB: " << iB << endl; 
		
		double fracEp=0; 
		double fracEm=0; 
		
		fracEp+=TMath::Power((grName[RN]->GetErrorYhigh(iB))/grName[RN]->GetY()[iB],2); 
		fracEm+=TMath::Power((grName[RN]->GetErrorYlow(iB))/grName[RN]->GetY()[iB],2); 
		
		if(N>0){
			fracEp+=TMath::Power(hName[bkg_systematic(nJ,bkg_nJ,mu,"stat",file)]->GetBinContent(ist),2); 
			fracEm+=TMath::Power(hName[bkg_systematic(nJ,bkg_nJ,mu,"stat",file)]->GetBinContent(ist),2); 
		}
		
		fracEp=TMath::Sqrt(fracEp); 
		fracEm=TMath::Sqrt(fracEm);
	
		
		double NcorEp=(1+fracEp)*Ncor; 
		double NcorEm=(1-fracEm)*Ncor; 
		
		
		if(st>=500){
			//cout << "st: " << st << " NE/N: " << NE/N << endl; 
			//cout << "st: " << st << " N: " << N << " Ncor : " << Ncor << " Ncor+: " << NcorEp << " Ncor-: " << NcorEm << endl; 
			
			hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->SetBinContent(ist,Ncor); 
			hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->SetBinError(ist,0);
			
			hName[expected_bkg(nJ,bkg_nJ,mu,"Ep",file)]->SetBinContent(ist, NcorEp); 
			hName[expected_bkg(nJ,bkg_nJ,mu,"Em",file)]->SetBinContent(ist, NcorEm);
			
			hName[bkg_systematic(nJ,bkg_nJ,mu,"Ep",file)]->SetBinContent(ist, fracEp); 
			hName[bkg_systematic(nJ,bkg_nJ,mu,"Em",file)]->SetBinContent(ist, fracEm);
			
			hName[bkg_systematic(nJ,bkg_nJ,mu,"Ep",file)]->SetBinError(ist,0);
			hName[bkg_systematic(nJ,bkg_nJ,mu,"Em",file)]->SetBinError(ist,0);
			
			hName[expected_bkg(nJ,bkg_nJ,mu,"Ep",file)]->SetBinError(ist,0); 
			hName[expected_bkg(nJ,bkg_nJ,mu,"Em",file)]->SetBinError(ist,0);
		}//fill for st>500 (where there should always be events)
		
	}//end of loop over ist
	
	hName[expected_bkg(nJ,bkg_nJ,mu,"",file)]->SetLineColor(kRed); 
	hName[expected_bkg(nJ,bkg_nJ,mu,"Ep",file)]->SetLineColor(kOrange); 
	hName[expected_bkg(nJ,bkg_nJ,mu,"Em",file)]->SetLineColor(kOrange);
	
}//end of bkg_est function

void compute_ratios(){
	for (std::map<TString,TFile*>::iterator it=fName.begin(); it!=fName.end(); it++){
		TString file=it->first; 
		if(file.Contains("RPV") || file.Contains("singleTop") || file.Contains("stealth") ) continue;
		if(file.Contains("diboson") ) continue;

		for(int mu=1; mu<=2; mu++){
			for(int num=3;num<=7; num++){

				//compute ratio to exclusive channels of 2-6 jets 
				
				for(int den=2; den<=6; den++){
					if(den>=num) continue; 
					TString RN=ratio_name(num,den,mu,file);
					quantile_ratio(hName[st_name(num,mu,"",file)],hName[st_name(den,mu,"",file)],
								   hName[st_name(num,mu,"_uW",file)],hName[st_name(den,mu,"_uW",file)],RN); 
					setcolor(grName[RN],num); 
					grName[RN]->GetYaxis()->SetTitle(Form("#alpha(%d-jets/%d-jets)",num,den));

				}//den loop
				
				
				//compute ratio to 2&3 jets combined 
				if(num<4) continue; 
				int den_comb=23; 
				TString RN=ratio_name(num,den_comb,mu,file); 
				//cout << "ratio name: "<< RN << endl; 
				quantile_ratio(hName[st_name(num,mu,"",file)],hName[st_name(den_comb,mu,"",file)],hName[st_name(num,mu,"_uW",file)],hName[st_name(den_comb,mu,"_uW",file)],RN); 
				setcolor(grName[RN],num); 
				grName[RN]->GetYaxis()->SetTitle(Form("#alpha(%d-jets/%d-jets)",num,den_comb));
				
				if(num<5) continue; 
				den_comb=234; 
				RN=ratio_name(num,den_comb,mu,file); 
				quantile_ratio(hName[st_name(num,mu,"",file)],hName[st_name(den_comb,mu,"",file)],hName[st_name(num,mu,"_uW",file)],hName[st_name(den_comb,mu,"_uW",file)],RN); 
				setcolor(grName[RN],num); 
				grName[RN]->GetYaxis()->SetTitle(Form("#alpha(%d-jets/%d-jets)",num,den_comb));
				
				
			}//num loop
			
		}//mu loop 
	}//file loop 
}

void quantile_ratio(TH1F *h_num, TH1F *h_den, TH1F *h_num_uW, TH1F *h_den_uW, TString ratio_name){
	//pass 3 histogrmas, numerator, denominator, and ratio
	//cout << "ratio name: " << ratio_name << endl; 
	TGraphAsymmErrors *g_rat = new TGraphAsymmErrors(h_den); 
	g_rat->SetName(ratio_name); 
	g_rat->GetXaxis()->SetTitle("S_{T} [GeV]"); 
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

	for(int i=0; i<g_rat->GetN(); i++){
		double x=g_rat->GetX()[i];
		double Num=h_num->GetBinContent(i+1);
		double Den=h_den->GetBinContent(i+1); 

		double y=Num/Den; 
		if(isnan(y)==1) y=0; 
		
		double exh=g_rat->GetErrorXhigh(i); 
		double exl=g_rat->GetErrorXlow(i); 
		
		double eyh=g_rat->GetErrorYhigh(i); 
		double eyl=g_rat->GetErrorYlow(i); 
		
		double yR=g_rat->GetY()[i]; 
		
		g_rat->SetPoint(i,x,y);
		g_rat->SetPointError(i,exl,exh,(eyl/yR)*y,(eyh/yR)*y); 
		
	}
	g_rat->SetTitle(""); 
	grName[ratio_name]=g_rat; 
	
}

TString st_name(int nJ, int Mu, TString weight, TString file){
	TString name=Form("st_nJets%d_%dMu",nJ,Mu); 
	name=name+weight+file; 
	return name; 
}

void combine_histograms(){
	//combine 2,3, and 4 jets weighted and un weighted for all signal and bkg samples 
	for(int iw=0; iw<=1; iw++){
		for(int Mu=1; Mu<=2; Mu++){
			for (std::map<TString,TFile*>::iterator it=fName.begin(); it!=fName.end(); it++){
				TString file=it->first; 
				TString weight="";
				if(iw==1)weight="_uW"; 
				
				TString hist_name=st_name(23,Mu,weight,file);
				
				//cout << "hist name: "<< hist_name << endl; 
				//cout << "2-jet: " << st_name(2,Mu,weight,file) << endl; 
				TH1F *h23=(TH1F*)hName[st_name(2,Mu,weight,file)]->Clone(hist_name); 
				h23->SetName(hist_name);
				h23->Add(hName[st_name(3,Mu,weight,file)]); 
				hName[h23->GetName()]=h23; 		
				
				hist_name=st_name(234,Mu,weight,file);
				TH1F *h234=(TH1F*)hName[st_name(2,Mu,weight,file)]->Clone(hist_name); 
				h234->SetName(hist_name);
				h234->Add(hName[st_name(3,Mu,weight,file)]); 
				h234->Add(hName[st_name(4,Mu,weight,file)]); 
				
				hName[h234->GetName()]=h234; 	
			}//loop over files
		}//loop over mu 
	}//weight
	
}

void open_file(TString name){
	cout << "open file:" << name << endl; 
	TString file_name="output_file_"+name+".root"; 
	//if(name=="wJets") file_name="output_file_wJets_inclusive.root";
	TFile *f = new TFile(file_name,"READ");
	if(f->IsOpen()!=1) {
		cout << "File: " << file_name << " Failed to open!" << endl; 
		return; 
	}
	name="_"+name;
	fName[name]=f; 
	GetHistograms(name);
	
}

void frac_bkg(TString variable, TString file){
	TString fraction_bkg=variable+file;
	fraction_bkg=fraction_bkg+"_frac"; 
	TH1F *h=(TH1F*)hName[variable+file]->Clone(fraction_bkg);
	hName[fraction_bkg]=h; 
	TString allMCname=variable+"_allMC";
	
	hName[fraction_bkg]->Divide(hName[allMCname]); 
	
}

void draw_frac(TString variable, TString canvasname){

	CreateCanvas(canvasname,"",600,600);
	CName[canvasname]->cd();
	
	hName[variable+"_dy_frac"]->SetLineColor(kBlue); 
	hName[variable+"_ttbar_frac"]->SetLineColor(kGreen);
	hName[variable+"_wJets_frac"]->SetLineColor(kRed); 
	
	hName[variable+"_ttbar_frac"]->SetMinimum(0);
	hName[variable+"_ttbar_frac"]->SetMaximum(1);
	hName[variable+"_ttbar_frac"]->DrawClone("histo");
	hName[variable+"_wJets_frac"]->DrawClone("histo same"); 
	hName[variable+"_dy_frac"]->DrawClone("histo same"); 
	
}

double get_N(double st1, double st2, TH1F *h){
	//cout << "name: " << h->GetName() << endl; 
	if(st2<0) st2=h->GetBinCenter(h->GetNbinsX()); 
	
	double N=0; 
	for(int ipt=1; ipt<=h->GetNbinsX(); ipt++){
		double st=h->GetBinCenter(ipt); 
		if(st>st1 && st<=st2) N+=h->GetBinContent(ipt); 
	}
	return N; 
}

void compute_correction(int nJ, double st1, double st2){
	cout << "nJ: " << nJ << " st: " << st1 << "-" << st2 << endl; 
	TString name=Form("st_nJets%d",nJ); 
	TString file="_singleMu"; 
	
	double Ndata[3];
	double Nw[3]; 
	double Ndy[3]; 
	double Nttbar[3];
	double NsingleTop[3]; 
	double Ndiboson[3];
	double NQCD[3]; 
	double Nother[3];
	
	get_column(name,file,st1,st2,Ndata); 
	file="_wJets"; 
	get_column(name,file,st1,st2,Nw); 
	file="_dy";
	get_column(name,file,st1,st2,Ndy); 
	file="_ttbar";
	get_column(name,file,st1,st2,Nttbar); 
	file="_singleTop";
	get_column(name,file,st1,st2,NsingleTop); 
	file="_diboson";
	get_column(name,file,st1,st2,Ndiboson); 
	file="_QCD";
	get_column(name,file,st1,st2,NQCD); 
	Nother[0]=NsingleTop[0]+Ndiboson[0]+NQCD[0];//contribution from "other" in each region 
	Nother[1]=NsingleTop[1]+Ndiboson[1]+NQCD[1];
	Nother[2]=NsingleTop[2]+Ndiboson[2]+NQCD[2];
	
	const int N=3; 
	
	TMatrixD NMC(N,N); 
	TVectorD aOther; aOther.Use(N,Nother);

	TVectorD aW; aW.Use(N,Nw);
	TVectorD aDY; aDY.Use(N,Ndy);
	TVectorD aTTBAR; aTTBAR.Use(N,Nttbar);
	TVectorD aDATA; aDATA.Use(N,Ndata);
	aDATA-=aOther;
	
	TVectorD Corr=aDATA;
	
	TMatrixDColumn(NMC,0)=aW;
	TMatrixDColumn(NMC,1)=aDY;
	TMatrixDColumn(NMC,2)=aTTBAR;

	TMatrixD invNMC=NMC.Invert();
	Corr*=invNMC;
	//Corr.Print();
	
	cout << "R(W): " << Corr(0) << " R(Z) " << Corr(1) << " R(ttbar) " << Corr(2) << endl; 
	
	MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)]=Corr(0); 
	MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)]=Corr(1); 
	MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)]=Corr(2); 
	for(int i=0; i<3; i++){
		Nevent[Form("N_%d_data_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Ndata[i];
		Nevent[Form("N_%d_W_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Nw[i];
		Nevent[Form("N_%d_DY_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Ndy[i];
		Nevent[Form("N_%d_ttbar_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Nttbar[i];
		Nevent[Form("N_%d_other_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Nother[i];

	}
	/*
	NMC.Print();
	
	cout << "Nw, Ndy, Nttbar" << endl; 
	cout << "{" << Nw[0] << "," << Ndy[0] << "," << Nttbar[0] << "}"<< endl; 
	cout << "{" << Nw[1] << "," << Ndy[1] << "," << Nttbar[1] << "}"<< endl; 
	cout << "{" << Nw[2] << "," << Ndy[2] << "," << Nttbar[2] << "}"<< endl; 
*/
	
}

void print_numbers(){

	int nJ=2;
	float st1=300;
	float st2=600; 
	
	ofstream output; 
	output.open("numbers_table.tex");
	output << "\\begin{tabular}{";
	int N=5; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}" << endl;  
	output << "Control & Data & W+Jets & DY & TTBar & Other \\\\ \\hline" << endl; 
	for(int i=0; i<3; i++){
		output << i+1 << " & "; 
		output << Nevent[Form("N_%d_data_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 
		output << Nevent[Form("N_%d_W_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 
		output << Nevent[Form("N_%d_DY_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 
		output << Nevent[Form("N_%d_ttbar_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 
		output << Nevent[Form("N_%d_other_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " \\\\ " << endl; 

	}
	
	output << "\\hline" << endl; 
	output.close(); 
	
}

void corrections_table(){

	ofstream output; 
	output.open("correction_table.tex");
	output << "\\begin{tabular}{";
	int N=5; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}" << endl;  
	output << "$N_{jets}$ & $S_{T}$ [GeV] & $R_{W}$ & $R_{DY}$ & $R_{tt}$ \\\\ \\hline" << endl; 
	for(int nJ=2; nJ<=4; nJ++){
		double st1=300;
		double st2=600;

		output << nJ << " & ";
		if(st2>0)output << st1 << "-" << st2 << " & "; 
		else output << "$\\ge$ " << st1 << " & "; 
		output << MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " \\\\ ";
		output << endl; 
		
		st1=600;
		st2=1000;
		
		output << nJ << " & ";
		if(st2>0)output << st1 << "-" << st2 << " & "; 
		else output << "$\\ge$ " << st1 << " & ";  
		output << MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " \\\\ ";
		output << endl; 
		
		st1=1000;
		st2=-1;
		
		output << nJ << " & ";
		if(st2>0)output << st1 << "-" << st2 << " & "; 
		else output << "$\\ge$ " << st1 << " & "; 
		
		output << MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " \\\\ \\hline";
		output << endl; 
		
	}
	output << "\\hline" << endl; 
	output.close(); 
}

void get_column(TString name, TString file,double st1,double st2, double *N){
	//name=st_nJets%d
	TString hist_name=name+"_W"+file; 
	double N1=get_N(st1,st2,hName[hist_name]); 
	hist_name=name+"_Z"+file;
	double N2=get_N(st1,st2,hName[hist_name]); 
	hist_name=name+"_1Mu_1El"+file; 
	double N3=get_N(st1,st2,hName[hist_name]); 
	N[0]=N1;//W enriched
	N[1]=N2; // Z enriched
	N[2]=N3; //ttbar enriched
}


void fill_stack(TString variable, bool blind){
	//modify this code to be given a string with the name
	//create just one stack per call
	//fix y axis min thing
	//plot h_mumu in bins
	cout << "variable: " << variable << endl; 

	bool print=false; 
	TLegend *L = new TLegend(0.6,0.5,0.8,0.7); 
	L->SetFillColor(10);
	L->SetLineColor(10);
	L->SetLineWidth(0);
	
	TString jet_canvas=variable+"_stack_comp"; 
	TString jet_stack=variable+"_hstack"; 
	CreateCanvas(jet_canvas,"stacked comparison",600,600); 
	CreateStack(jet_stack,variable);  
	
	float st_max=3000;

	if(print) cout << "QCD" << endl; 
	TString file="_QCD";
	TString name;
	name=variable;
	double max_x=hName[name+file]->GetBinCenter(hName[name+file]->GetNbinsX())+hName[name+file]->GetBinWidth(hName[name+file]->GetNbinsX())/2; 
	double min_x=hName[name+file]->GetBinLowEdge(1); 
	
	//hName[name+file]->SetAxisRange(0,st_max);
	hName[name+file]->SetFillColor(7); 
	hName[name+file]->SetMarkerColor(7); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"_QCD");
	
	stackName[jet_stack]->Add(hName[name+file]); 
	TString tot_MC="MC_stack_"+variable; 
	
	TH1F *h=(TH1F*)hName[name+file]->Clone(tot_MC);
	hName[h->GetName()]=h;
	
	if(print)cout << "diboson" << endl; 
	file="_diboson";
	//TString name=variable+Form("_nJets%d_GeV",nJets);
	hName[name+file]->SetFillColor(9); 
	hName[name+file]->SetMarkerColor(2); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"Diboson"); 
	
	stackName[jet_stack]->Add(hName[name+file]);
	hName[tot_MC]->Add(hName[name+file]);
	
	if(print) cout <<"singleTop" << endl;
	file="_singleTop";
	//TString name=variable+Form("_nJets%d_GeV",nJets);
	
	hName[name+file]->SetFillColor(5); 
	hName[name+file]->SetMarkerColor(2); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"Single Top"); 
	
	stackName[jet_stack]->Add(hName[name+file]);
	hName[tot_MC]->Add(hName[name+file]);
	
	
	if(print) cout << "dy"  << endl; 
	file="_dy";
	//TString name=variable+Form("_nJets%d_GeV",nJets);
	
	hName[name+file]->SetFillColor(4); 
	hName[name+file]->SetMarkerColor(2); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"Drell-Yan"); 
	
	stackName[jet_stack]->Add(hName[name+file]);
	hName[tot_MC]->Add(hName[name+file]);
	
	if(print) cout << "ttbar" << endl; 
	file="_ttbar";

	hName[name+file]->SetFillColor(3); 
	hName[name+file]->SetMarkerColor(3); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"TTbar"); 
	
	stackName[jet_stack]->Add(hName[name+file]); 
	hName[tot_MC]->Add(hName[name+file]);
	if(print) cout << "wJets" << endl; 
	file="_wJets";
	
	hName[name+file]->SetFillColor(2); 
	hName[name+file]->SetMarkerColor(4); 
	hName[name+file]->SetMarkerStyle(21); 
	hName[name+file]->SetMarkerSize(0.5);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"W+Jets"); 
	
	stackName[jet_stack]->Add(hName[name+file]); 
	hName[tot_MC]->Add(hName[name+file]);
	
	//stackName[jet_stack]->Add(hName[name+file]); 
	
	if(print) cout << "RPV300" << endl; 
	file="_RPV300";
	
	hName[name+file]->SetLineStyle(1);
	hName[name+file]->SetLineWidth(2);
	hName[name+file]->SetLineColor(kOrange);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"RPV: M_{#tilde{t}}=300 GeV"); 
	
	if(print) cout << "RPV500" << endl; 
	file="_RPV500";
	
	hName[name+file]->SetLineStyle(1);
	hName[name+file]->SetLineWidth(2);
	hName[name+file]->SetLineColor(kRed);
	hName[name+file]->SetAxisRange(min_x,max_x);
	L->AddEntry(hName[name+file],"RPV: M_{#tilde{t}}=500 GeV"); 
	/*
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
	*/ 
	if(print) cout << "single mu: " << endl; 
	file="_singleMu";
	
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
	if(variable.Contains("mt")==0)	gPad->SetLogy();
	//if(nJets<=5){
//	hName[name+file]->GetXaxis()->SetTitle("S_{T} [GeV]"); 
//	hName[name+file]->GetYaxis()->SetTitle("Events/GeV"); 
	
	TLatex BL(0.36,0.815,hName[name+file]->GetTitle()); 
	hName[name+file]->SetTitle(""); 
	BL.SetNDC(kTRUE); 
	BL.SetTextSize(0.05); 
	//cout << variable << endl; 
	if(!blind){		
		//cout << "not blinded: "<< endl; 
		hName[name+file]->DrawCopy(" E1 ");
		stackName[jet_stack]->Draw("hist same");
		hName[name+file]->DrawCopy(" E1 same");

	}
	else{
		//cout << "blinded: " << endl; 
		stackName[jet_stack]->Draw("hist");
	}
	
	hName[name+"_RPV500"]->DrawCopy("hist same"); 
	hName[name+"_RPV300"]->DrawCopy("hist same"); 

	
	L->DrawClone("same"); 
	BL.DrawClone(); 

	
	ratio_pad->cd(); 
	

	hName[name+file]->Divide(hName[tot_MC]); 
	hName[name+file]->GetYaxis()->SetTitle("Data/MC"); 
	hName[name+file]->SetMinimum(0.5);
	hName[name+file]->SetMaximum(1.5);
	if(!blind)hName[name+file]->DrawCopy(); 
	//else BL.DrawClone(); 
	
}
void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
}
void CreateStack(TString Name,TString Title){
	THStack *hstack = new THStack(Name,Title); 
	stackName[Name]=hstack;
}

void draw_headersim(){
	
	L1sim.SetNDC(kTRUE); 
	
	L1sim.DrawLatex(0.15,0.92, cms_sim); 
	
}

void draw_header(){
	
	L1.SetNDC(kTRUE); 
	L2.SetNDC(kTRUE); 
	L1.SetTextSize(0.03);
	L2.SetTextSize(0.03); 
	L1.DrawLatex(0.15,0.92, cms_pre); 
	L2.DrawLatex(0.43,0.92, lumi); 
	
}

void GetHistograms(TString file){
	TH1F *TH1F_names=(TH1F*)fName[file]->FindObjectAny("TH1F_names"); 
	
	
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names") continue; 
		//cout <<" filename: "<< file <<  " hist_name: "<< name << endl; 

		TH1F *h=(TH1F*)fName[file]->FindObjectAny(name); 
		h->SetStats(kFALSE); 
		h->SetName(h->GetName()+file); 
		hName[h->GetName()]=h; 
	}
	delete TH1F_names; 
	//cout << "end of GetHistograms() " << endl; 
}

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, const Float_t* xBins)
{
	TH1F* h = new TH1F(name, title, nBinsX, xBins);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	h->SetStats(kFALSE); 
	hName[name] = h;
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