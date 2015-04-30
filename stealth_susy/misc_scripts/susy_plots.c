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

	initialze_correction_bins();
	bookHisto(); 
	

	for(int nJ=2; nJ<=4; nJ++){
		for(int ist=0; ist<NbinsCor; ist++){
			float st1=correction_bins[ist];
			float st2=correction_bins[ist+1];
			compute_correction(nJ,st1,st2);
		}
		
		//print_numbers(nJ,300,500);
		//print_numbers(nJ,500,800);
		//print_numbers(nJ,800,-1);

	}
	fit_corrections();

	table_contamination(1);
	table_contamination(2);

//	corrections_table();

	draw_corrections();
	
	/*
	rescaleMC(3,1);
	rescaleMC(3,2);
	
	rescaleMC(4,1);
	rescaleMC(4,2);
	*/
	rescaleMC(5,1);
	rescaleMC(5,2);
	
	rescaleMC(6,1);
	rescaleMC(6,2);
	
	rescaleMC(7,1);
	rescaleMC(7,2);
	
	rescaleMC(8,1);
	rescaleMC(8,2);
	
	for(int nJ=2;nJ<=8; nJ++){
		for(int mu=1; mu<=2; mu++){	
			add_corrected_bkg(nJ,mu); 
			relative_systematic(nJ,mu);
		}
	}
	
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
	
	for(int nJ=1;nJ<=8; nJ++){
		bool blind=false; 
		//if(nJ>6) blind=true; 
		for(int mu=1; mu<=2; mu++){	
			for (std::map<TString,TFile*>::iterator it=fName.begin(); it!=fName.end(); it++) {
				if(it->first!="_singleMu")frac_bkg(Form("st_nJets%d_%dMu_GeV",nJ,mu),it->first); 
			}
		}
		
		draw_frac(Form("st_nJets%d_1Mu_GeV",nJ),Form("bkg_frac_nJets%d_st1mu",nJ)); 
		draw_frac(Form("st_nJets%d_2Mu_GeV",nJ),Form("bkg_frac_nJets%d_st2mu",nJ)); 

		fill_stack(Form("st_nJets%d_1Mu",nJ),blind); 
		fill_stack(Form("st_nJets%d_2Mu",nJ),blind); 
		fill_stack(Form("st_nJets%d_1Mu_1El",nJ),blind);
		
		fill_stack(Form("st_nJets%d_W",nJ),blind); 
		fill_stack(Form("st_nJets%d_Z",nJ),blind); 

		for(int ist=1;ist<=3;ist++){
			for(int ij=0; ij<nJ; ij++) {
				fill_stack(Form("W_jetpt_jet%d_st%d_nJets%d",ij,ist,nJ),0); 
				fill_stack(Form("Z_jetpt_jet%d_st%d_nJets%d",ij,ist,nJ),0); 
				fill_stack(Form("tt_jetpt_jet%d_st%d_nJets%d",ij,ist,nJ),0); 
				
			}
		}
		
		fill_stack(Form("Z_nJets%d_dileptonpT",nJ),0);
		fill_stack(Form("tt_nJets%d_dileptonpT",nJ),0);
		fill_stack(Form("W_nJets%d_muonpT",nJ),0);
		for(int mu=1; mu<=2; mu++){
			fill_stack(Form("Z_nJets%d_muonpT%d",nJ,mu),0);
			fill_stack(Form("tt_nJets%d_muonpT%d",nJ,mu),0);
		}
		
		fill_stack(Form("met_nJets%d_1Mu_1El",nJ),0);
		fill_stack(Form("met_nJets%d_W",nJ),0);
		fill_stack(Form("met_nJets%d_Z",nJ),0);

		fill_stack(Form("ht_nJets%d_1Mu_1El",nJ),0);
		fill_stack(Form("ht_nJets%d_W",nJ),0);
		fill_stack(Form("ht_nJets%d_Z",nJ),0);
		
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
	for(int nJ=1; nJ<=nJetsmax; nJ++){
		fill_stack(Form("jetpTsys%d_W",nJ),0);
		fill_stack(Form("jetpTsys%d_Z",nJ),0);
		fill_stack(Form("jetpTsys%d_1Mu_1El",nJ),0);
	}
	
	fill_stack("jetpT1_1Mu",0);
	fill_stack("jetpT2_1Mu",0);
	
	fill_stack("jetpT1_W",0);
	fill_stack("jetpT1_Z",0);
	fill_stack("jetpT2_W",0);
	fill_stack("jetpT2_Z",0);
	
	fill_stack("jetpT1_1Mu_1El",0);
	fill_stack("jetpT2_1Mu_1El",0);
	
	//fill_stack("jetpT3_1Mu",0);
	//fill_stack("jetpT4_1Mu",0);
	
	fill_stack("jetpT1_2Mu",0);
	fill_stack("jetpT2_2Mu",0);
	//fill_stack("jetpT3_2Mu",0);
	//fill_stack("jetpT4_2Mu",0);
	
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
	jet_meanpt();

		for(int nJ=4; nJ<=nJetsmax; nJ++){
			for(int mu=1; mu<=2; mu++){
				CreateHistogram(acc_name(nJ,mu),"","M_{#tilde{t}} [GeV]","acc #times #epsilon",8,250,1050); 
				for(int irpv=300; irpv<=1000;irpv+=100){
				TString signal=Form("_RPV%d",irpv);
				compute_sensitivity(nJ,mu,signal); 
				compute_acceptance(nJ,mu,irpv); 
				compute_expected_obs(nJ, mu, irpv);
			}//irpv loop
		}//nJ loop
	}//irpv loop 
	
	//now make inclusive datacard
	for(int irpv=300; irpv<=1000;irpv+=100) compute_expected_obs(irpv);
	
	for(int nJ=5;nJ<=nJetsmax;nJ++){
		for(int mu=1; mu<=2; mu++){
			if(nJ==5 && mu==1) continue; 
			table_obs(nJ,mu);
		}
	}
	
	write();
}

void initialze_correction_bins(){
	NbinsCor=5; 
	correction_bins[0]=300;
	correction_bins[1]=400;
	correction_bins[2]=500;
	correction_bins[3]=600;
	correction_bins[4]=700;
	correction_bins[5]=3000;
	
}

void close(){
	for (std::map<TString,TFile*>::iterator it=fName.begin(); it!=fName.end(); it++) {
		fName[it->first]->Close();
	}
}

void write(){
	TFile *output_file = new TFile("susy_plots.root","RECREATE"); 
	TFile *output_histogram = new TFile("susy_histograms.root","RECREATE"); 
	output_file->cd();
	//output_file->mkdir("histograms"); 
	//output_file->cd("histograms"); 
	bool write_=true; 
	output_histogram->cd(); 
	
	TH1F *hist_names = new TH1F("TH1F_names","",1,0,1); 
	
	for (std::map<TString,TH1F*>::iterator it=hName.begin(); it!=hName.end(); it++) {
		//if(it->first.Contains("23"))cout << it->first << endl; 
		//cout << hName[it->first]->GetName() << endl; 
		//if(write_ && it->first.Contains("st") && it->first.Contains("Mu"))hName[it->first]->Write();
		if(it->first.Contains("acc") || it->first.Contains("sb") || it->first.Contains("observed") || it->first.Contains("expected")){
			hist_names->Fill(it->first,1);
			hName[it->first]->Write();
		}
		
	}
	hist_names->Write();
	output_file->mkdir("graphs"); 
	output_file->cd("graphs"); 
	for (std::map<TString,TGraphAsymmErrors*>::iterator it=grName.begin(); it!=grName.end(); it++) {
		if(write_)grName[it->first]->Write();
	}
	output_file->cd(); 
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		if(write_)CName[it->first]->Write();
		TString name="/uscms_data/d3/btcarlso/Figures_SUSYAN/"; 
		name=name+it->first+".pdf"; 
		//CName[it->first]->Print(name);
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
	
	for(int num=5; num<=nJetsmax; num++){
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
	open_file("UDD300"); 
//	open_file("stealth_ww500_300");
//	open_file("stealth_ww1400_300");
	open_file("wJets");
	open_file("dy");
	open_file("ttbar");
	open_file("singleTop");
	open_file("diboson"); 
	open_file("singleMu");
	open_file("allMC"); 
	
	Ngen["_UDD300"]=98000;
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
	if(num==5) gr->SetMarkerColor(419); 
	if(num==3) gr->SetMarkerColor(kBlack); 
	
}

void compare_ratios(int num, int den,int mu){
	//cout << "Ratio comparison: " << num << "-" << den << endl; 
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
TString sig_bkg(int nJ,int mu,TString file){
	TString name= Form("sb_nJets%d_%dMu",nJ,mu);
	name=name+file; 
	return name; 
}

void compute_expected_obs(int iRPV){
	//compute expected/observed number of events in a combined datacard for each mass point
	
	TString signal=Form("_RPV%d",iRPV); 
	ofstream datacard; 
	string file_name=Form("datacard_combined_RPV%d.txt",iRPV); 
	datacard.open(file_name.c_str());
	
	datacard << "imax 7 number of channels " << endl; 
	datacard << "jmax 1 number of backgrounds " << endl; 
	datacard << "kmax 1 number of nuisance parameters (sources of systematical uncertainties) " << endl; 
	datacard << "------------" << endl; 
	
	//fill an array with the total number of expected,observed and signal events for each bin
	
	double Nexp[7];
	double Nobs[7];
	double Nsig[7]; 
	double Nexp_sys[7];
	
	int ibin=0; 
	
	for(int mu=1; mu<=2; mu++){
		for(int nJ=5; nJ<=8; nJ++){
			if(nJ==5 && mu==1) continue; 
			//initialize arrays
			Nexp[ibin]=0;
			Nobs[ibin]=0;
			Nsig[ibin]=0;
			Nexp_sys[ibin]=0;
			double Nexp_systmp=0;
			
			TString name_data=st_name(nJ,mu,"","_singleMu");
			TString name_MC=st_name(nJ,mu,"","_correctedMC");
			TString name_sys=Form("st_nJets%d_%dMu_systematic",nJ,mu);
			int Nbins=hName[name_data]->GetNbinsX(); //number of bins
			int isto=hName[sig_bkg(nJ,mu,signal)]->GetMaximumBin(); //st threshold bin
			//if(isto=Nbins)isto=Nbins-1;
			
			for(int ist=isto; ist<=Nbins; ist++){
				Nexp[ibin]+=hName[name_MC]->GetBinContent(ist); 
				Nobs[ibin]+=hName[name_data]->GetBinContent(ist); 
				Nsig[ibin]+=hName[st_name(nJ,mu,"",signal)]->GetBinContent(ist); 
				Nexp_systmp+=hName[name_sys]->GetBinContent(ist); 
			}
			cout << "nJ: " << nJ << " mu " << mu << " Nobssys: " << Nexp_systmp << " obs: "<< Nexp[ibin]<< endl; 
			Nexp_sys[ibin]=1+Nexp_systmp/Nexp[ibin]; 
			ibin++;
		}//nJ loop
	}
	//fill observation
	datacard << "bin b1 b2 b3 b4 b5 b6 b7 " << endl;
	datacard << "observation "; 
	
	for(int bin=0; bin<7; bin++){
		datacard << Nobs[bin] << " ";
	}
	datacard << endl; 
	datacard << "------------" << endl; 
	
	datacard << "bin b1 b1 b2 b2 b3 b3 b4 b4 b5 b5 b6 b6 b7 b7 " << endl;
	datacard << "process sig bkg sig bkg sig bkg sig bkg sig bkg sig bkg sig bkg " << endl; 
	datacard << "process ";
	for(int bin=0; bin<7; bin++) datacard << " 0 1 ";
	datacard << endl; 
	
	//fill signal expected and background expected 
	datacard << "rate     "; 
	for(int bin=0; bin<7; bin++){
		datacard << Nsig[bin] << " " <<  Nexp[bin] << " ";
	}
	
	datacard <<endl; 
	datacard << "------------" << endl; 
	datacard << " bkg   lnN  ";
	
	for(int bin=0; bin<7; bin++){
		datacard << "- " << Nexp_sys[bin] << " "; 
	}
	datacard <<endl; 
}

void compute_expected_obs(int nJ, int mu, int iRPV){
	TString signal=Form("_RPV%d",iRPV); 
	ofstream datacard; 
	string file_name=Form("datacard_nJets%d_%dMu",nJ,mu); 
	file_name=file_name+signal;
	file_name=file_name+".txt";
	datacard.open(file_name.c_str());
	datacard << "imax 1 number of channels " << endl; 
	datacard << "jmax 1 number of backgrounds " << endl; 
	datacard << "kmax 1 number of nuisance parameters (sources of systematical uncertainties) " << endl; 
	datacard << "bin b1 " << endl;
	datacard << "------------" << endl; 
	bool print=false; 
	if(print)cout << "nJ: " << nJ << " mu: " << mu << " signal: " << signal << endl; 
	TString name_data=st_name(nJ,mu,"","_singleMu");
	TString name_MC=st_name(nJ,mu,"","_correctedMC");
	TString name_sys=Form("st_nJets%d_%dMu_systematic",nJ,mu);
	
	int Nbins=hName[name_data]->GetNbinsX(); //number of bins
	int isto=hName[sig_bkg(nJ,mu,signal)]->GetMaximumBin(); //st threshold bin
	//if(isto=Nbins)isto=Nbins-1;

	cout << "nJ: " << nJ << " mu: " << mu << " stThreshold: "<< hName[name_data]->GetBinLowEdge(isto) << endl; 

	double Nexp=0; 
	double Nobs=0;
	double Nsig=0; 
	double Nexp_sys=0; 
	for(int ist=isto; ist<=Nbins; ist++){
		Nexp+=hName[name_MC]->GetBinContent(ist); 
		Nobs+=hName[name_data]->GetBinContent(ist); 
		Nexp_sys+=hName[name_sys]->GetBinContent(ist);
		Nsig+=hName[st_name(nJ,mu,"",signal)]->GetBinContent(ist); 
	}
	
	hName[Form("expected_sys_mass_nJets%d_%dMu",nJ,mu)]->Fill(iRPV,Nexp_sys); 
	hName[Form("expected_mass_nJets%d_%dMu",nJ,mu)]->Fill(iRPV,Nexp); 
	hName[Form("observed_mass_nJets%d_%dMu",nJ,mu)]->Fill(iRPV,Nobs); 
	hName[Form("signal_mass_nJets%d_%dMu",nJ,mu)]->Fill(iRPV,Nsig); 
	
	Nexp_sys=1+Nexp_sys/Nexp;//systematic is supposed to be deltaX/x
	
	if(print)cout << "Nexp: " << Nexp << " Nobs: " << Nobs << " Nsig: " << Nsig << endl; 
	


	
	datacard << "observation " << Nobs << endl; 
	datacard << "------------" << endl; 
	datacard << "bin     b1     b1" << endl;
	datacard << "process sig    bkg" << endl; 
	datacard << "process 0      1" << endl;
	datacard << "rate   " << Nsig << "      " << Nexp << endl; 
	datacard << "-------------" << endl; 
	datacard << "bkg  lnN  -      " << Nexp_sys << endl; 
	datacard.close();
	
	
}

void compute_sensitivity(int nJ, int mu, TString signal){
	//cout << "sensitivity (mu): " << mu << " nJ: " << nJ << " file: " <<signal << endl; 
	TString file = "_singleMu"; 
	int bkg_nJ=3; 
	float st_bins[]={300,400,500,600,700,800,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	
	CreateHistogram(sig_bkg(nJ,mu,signal),"","S_{T}^{thresh} [GeV]","S/#sqrt(B)",N,st_bins); 
	
	int Nbins=hName[st_name(nJ,mu,"",signal)]->GetNbinsX(); 
	for(int ist=1; ist<=Nbins; ist++){
		
		double bkg=hName[st_name(nJ,mu,"","_allMC")]->Integral(ist,Nbins); 
		double sig=hName[st_name(nJ,mu,"",signal)]->Integral(ist,Nbins);
		
		if(isnan(sig/TMath::Sqrt(bkg))==1 || isinf(sig/TMath::Sqrt(bkg))==1) continue; 
		//cout << "ist: " << ist << sig/TMath::Sqrt(bkg) << endl; 
		hName[sig_bkg(nJ,mu,signal)]->SetBinContent(ist,sig/TMath::Sqrt(bkg)); 
		hName[sig_bkg(nJ,mu,signal)]->SetBinError(ist,0); 
	}
	
}

void jet_meanpt(){

	for(int nJ=1; nJ<=5; nJ++){
		int ijet=0;
		TString variable=Form("W_jetpt_jet%d_st1_nJets%d",ijet,nJ);
		TString file="_singleMu"; 
		double data_mean=hName[variable+file]->GetMean();
		//cout << "data name: "<< hName[variable+file]->GetName() << endl; 
		file="_allMC";
		double mc_mean=hName[variable+file]->GetMean(); 
		//cout << "mc name: " << hName[variable+file]->GetName() << endl; 
		
		//cout << "nJ: "<< nJ << " <data> " << data_mean << " <mc> " << mc_mean << " delta(data-mc): " << data_mean-mc_mean<< endl; 

	}
	
}

TString acc_name(int nJ, int mu){
	TString name=Form("acc_nJets%d_%dMu",nJ,mu);
	return name; 
}

void table_obs(int nJ, int mu){
	ofstream output; 
	output.open(Form("observed_expected_nJets%d_%dMu.tex",nJ,mu)); 
	output << "\\begin{tabular}{";
	int N=5; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}" << endl; 
	
	output << "$M_{\\tilde{t}}$ [GeV] & $S_{T}^{Thresh}$ & $N_{sig}$ & $N_{bkg} \\pm N_{sys}$ & $N_{obs}$ \\\\" << endl; 
	for(int irpv=300; irpv<=1000; irpv+=100){
		TString signal=Form("_RPV%d",irpv); 
		int isto=hName[sig_bkg(nJ,mu,signal)]->GetMaximumBin(); //st threshold bin
		double st_threshold=hName[sig_bkg(nJ,mu,signal)]->GetBinLowEdge(isto);
		int iMass=hName[Form("expected_sys_mass_nJets%d_%dMu",nJ,mu)]->FindBin(irpv); 
		double Nexp_sys=hName[Form("expected_sys_mass_nJets%d_%dMu",nJ,mu)]->GetBinContent(iMass);
		double Nobs=hName[Form("observed_mass_nJets%d_%dMu",nJ,mu)]->GetBinContent(iMass);
		double Nexp=hName[Form("expected_mass_nJets%d_%dMu",nJ,mu)]->GetBinContent(iMass);
		double Nsig=hName[Form("signal_mass_nJets%d_%dMu",nJ,mu)]->GetBinContent(iMass);
		
		output << irpv << " & " << st_threshold << " & " << Nsig << " & " << Nexp << " $\\pm $ " << Nexp_sys; 
		output << " & " << Nobs << "\\\\" << endl; 
			
	}
	output << "\\end{tabular}"; 
	output.close();
	
}

void table_contamination(int mu){
//compute_contamination(int mu, int nJ, string signal, float st1, float st2)	
	float st1=300;
	float st2=500; 
	ofstream output; 
	output.open(Form("contamination_%dMu.tex",mu)); 
	output << "\\begin{tabular}{";
	int N=6; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}" << endl;  
	output << "M [GeV] & 3-jets & 4-jets & 5-jets & 6-jets & 7-jets \\\\ \\hline" << endl; 
	for(int irpv=300; irpv<=1000; irpv+=100){
		output << irpv << " & ";
		TString signal=Form("_RPV%d",irpv); 
	
		for(int nJ=3; nJ<=nJetsmax; nJ++) output << compute_contamination(mu, nJ, signal, st1, st2) << " & "; 
		output << " \\\\ " << endl; 
	}

	output << "\\end{tabular}"; 
}

double compute_contamination(int mu, int nJ, TString signal, float st1, float st2){

	double bkg=0; 
	double sig=0; 
	for(int ist=1; ist<=hName[st_name(nJ,mu,"","_allMC")]->GetNbinsX();ist++){
		double st=hName[st_name(nJ,mu,"","_allMC")]->GetBinCenter(ist); 
		if(st>st1 && st<=st2) {
			bkg+=hName[st_name(nJ,mu,"","_allMC")]->GetBinContent(ist); 
			sig+=hName[st_name(nJ,mu,"",signal)]->GetBinContent(ist);
		}
	}
	return sig/(sig+bkg); 
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

void compute_ratios(){
	for (std::map<TString,TFile*>::iterator it=fName.begin(); it!=fName.end(); it++){
		TString file=it->first; 
		if(file.Contains("RPV") || file.Contains("singleTop") || file.Contains("stealth") ) continue;
		if(file.Contains("diboson") ) continue;

		for(int mu=1; mu<=2; mu++){
			for(int num=3;num<=nJetsmax; num++){

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
	hName[variable+"_ttbar_frac"]->SetLineColor(419);
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

double get_NE(double st1, double st2, TH1F *h){
	//cout << "name: " << h->GetName() << endl; 
	if(st2<0) st2=h->GetBinCenter(h->GetNbinsX()); 
	
	double NE=0;
	for(int ipt=1; ipt<=h->GetNbinsX(); ipt++){
		double st=h->GetBinCenter(ipt); 
		if(st>st1 && st<=st2) NE+=TMath::Power(h->GetBinError(ipt),2); 
	}
	return TMath::Sqrt(NE); 
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
	
	double NdataE[3];
	double NwE[3]; 
	double NdyE[3]; 
	double NttbarE[3];
	double NsingleTopE[3]; 
	double NdibosonE[3];
	double NQCDE[3]; 
	double NotherE[3];
	
	get_column(name,file,st1,st2,Ndata,NdataE); 
	file="_wJets"; 
	get_column(name,file,st1,st2,Nw,NwE); 
	file="_dy";
	get_column(name,file,st1,st2,Ndy,NdyE); 
	file="_ttbar";
	get_column(name,file,st1,st2,Nttbar,NttbarE); 
	file="_singleTop";
	get_column(name,file,st1,st2,NsingleTop,NsingleTopE); 
	file="_diboson";
	get_column(name,file,st1,st2,Ndiboson,NdibosonE); 
	file="_QCD";
	get_column(name,file,st1,st2,NQCD,NQCDE); 
	Nother[0]=NsingleTop[0]+Ndiboson[0]+NQCD[0];//contribution from "other" in each region 
	Nother[1]=NsingleTop[1]+Ndiboson[1]+NQCD[1];
	Nother[2]=NsingleTop[2]+Ndiboson[2]+NQCD[2];
	
	for(int i=0; i<3; i++) {
		NotherE[i]=TMath::Power(NsingleTopE[i],2)+TMath::Power(NdibosonE[i],2)+TMath::Power(NQCDE[i],2);
		NotherE[i]=TMath::Sqrt(NotherE[i]); 
	}
	
	const int N=3; 
	
	TH1F *W     = new TH1F("W"    ,"W"    , 3, 0.5, 3.5);
	TH1F *Z     = new TH1F("Z"    ,"Z"    , 3, 0.5, 3.5);
	TH1F *tt    = new TH1F("tt"   ,"tt"   , 3, 0.5, 3.5);
	TH1F *data  = new TH1F("data" ,"data" , 3, 0.5, 3.5);
	
	TH1F *weightW  = new TH1F("W_weight" , "W_weight" , 3, 0.5, 3.5);	
	TH1F *weightZ  = new TH1F("Z_weight" , "Z_weight" , 3, 0.5, 3.5);	
	TH1F *weighttt = new TH1F("tt_weight", "tt_weight", 3, 0.5, 3.5);	
	
	for(int i=0; i<3; i++){
		Nevent[Form("N_%d_data_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Ndata[i];
		Nevent[Form("N_%d_W_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Nw[i];
		Nevent[Form("N_%d_DY_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Ndy[i];
		Nevent[Form("N_%d_ttbar_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Nttbar[i];
		Nevent[Form("N_%d_other_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=Nother[i];
		
		NeventE[Form("N_%d_data_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=NdataE[i];
		NeventE[Form("N_%d_W_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=NwE[i];
		NeventE[Form("N_%d_DY_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=NdyE[i];
		NeventE[Form("N_%d_ttbar_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=NttbarE[i];
		NeventE[Form("N_%d_other_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)]=NotherE[i];
		cout << Form("Nw[%d]",i) << Nw[i] << " +/- " <<  NwE[i] <<  endl; 
		double Nw_prime     = TMath::Power(Nw[i]/NwE[i]        , 2); 
		double NZ_prime     = TMath::Power(Ndy[i]/NdyE[i]        , 2); 
		double Ntt_prime    = TMath::Power(Nttbar[i]/NttbarE[i]      , 2); 
		
		if(isnan(Nw_prime)==1) Nw_prime=0;
		if(isnan(NZ_prime)==1) NZ_prime=0;
		if(isnan(Ntt_prime)==1) Ntt_prime=0;

		double W_weight=Nw[i]/Nw_prime;
		double Z_weight=Ndy[i]/NZ_prime;
		double tt_weight=Nttbar[i]/Ntt_prime; 
		
		if(isnan(W_weight)==1) W_weight=1; 
		if(isnan(Z_weight)==1) Z_weight=1; 
		if(isnan(tt_weight)==1) tt_weight=1; 

		W->SetBinContent(i+1,Nw_prime);
		Z->SetBinContent(i+1,NZ_prime);
		tt->SetBinContent(i+1,Ntt_prime);
		data->SetBinContent(i+1,Ndata[i]-Nother[i]); 
		weightW->SetBinContent(i+1,W_weight); 
		weightZ->SetBinContent(i+1,Z_weight); 
		weighttt->SetBinContent(i+1,tt_weight); 

	}
	
	TObjArray *mc= new TObjArray(3);
	mc->Add(W);
	mc->Add(Z);
	mc->Add(tt);
	TFractionFitter *fit = new TFractionFitter(data,mc,"Q");
	fit->SetWeight(0,weightW);
	fit->SetWeight(1,weightZ);
	fit->SetWeight(2,weighttt);
	/*
	for(int i=1; i<=3; i++){
		cout << "data: " << data->GetBinContent(i) << " W: " << W->GetBinContent(i) << " Z: " << Z->GetBinContent(i) << " tt " << tt->GetBinContent(i) << endl; 
		cout << "weights, W: " << weightW->GetBinContent(i) << " Z: " << weightZ->GetBinContent(i) << " tt " << weighttt->GetBinContent(i) << endl; 
	}
	*/
	for(int i=0; i<3; i++){
		cout << Form("Nw[%d]=%.1f",i,Nw[i]) << Form(" Nz[%d]=%.1f",i,Ndy[i]) << Form(" Ntt[%d]=%.1f",i,Nttbar[i])  << Form(" Nother[%d]=%.1f",i,Nother[i]);
		cout << Form(" Ndata[%d]=%.1f",i,Ndata[i]) << endl;
	}
	fit->Fit();
		
	double R[3]; 
	double RE[3];
	
	double Fo[3];
	double Nwtot=0; 
	double Nztot=0; 
	double Nttbartot=0; 
	double Ndatatot=0; 
	for(int i=0; i<3; i++){
		Nwtot+=Nw[i]; 
		Nztot+=Ndy[i]; 
		Nttbartot+=Nttbar[i]; 
		Ndatatot+=Ndata[i]-Nother[i];
	}
	Fo[0]=Nwtot/Ndatatot; 
	Fo[1]=Nztot/Ndatatot;
	Fo[2]=Nttbartot/Ndatatot; 
	
	for(int i=0; i<3; i++){
		double x=0;
		double xE=0;
		fit->GetResult(i,x,xE);
		
		R[i]=x/Fo[i];
		RE[i]=(xE/x)*R[i] ;
		cout << "R: " << R[i] << " +/- " << (xE/x)*R[i] << endl; 
		
	}
	
	delete W;
	delete Z;
	delete tt;
	delete data; 
	delete weightW;
	delete weightZ;
	delete weighttt;
	
	TMatrixD NMC(N,N); 
	TVectorD aOther; aOther.Use(N,Nother);

	TVectorD aW; aW.Use(N,Nw);
	TVectorD aDY; aDY.Use(N,Ndy);
	TVectorD aTTBAR; aTTBAR.Use(N,Nttbar);
	TVectorD aDATA; aDATA.Use(N,Ndata);
	//cout << "print before subtracting other: " << endl; 
	//cout << "Data[0]: " << Ndata[0] << endl; 

	//aDATA.Print();
	//cout << "print after subtracting other: " << endl; 
	aDATA-=aOther;
	//aDATA.Print();
	//cout << "Data[0]: " << Ndata[0] << endl; 
	
	TVectorD Corr=aDATA;
	
	TMatrixDColumn(NMC,0)=aW;
	TMatrixDColumn(NMC,1)=aDY;
	TMatrixDColumn(NMC,2)=aTTBAR;

	TMatrixD invNMC=NMC.Invert();
	Corr*=invNMC;
	//Corr.Print();
	
	cout << "R(W): " << Corr(0) << " R(Z) " << Corr(1) << " R(ttbar) " << Corr(2) << endl; 
	
	int ist=hName["h_Rw_nJ2"]->FindBin((st1+st2)/2)-1; 
	cout << "ist: " << ist << endl; 
	//fill nJ binned
	hName[Form("h_Rw_st%d",ist)]->SetBinContent(nJ-1,R[0]) ; 
	hName[Form("h_Rw_st%d",ist)]->SetBinError(nJ-1,RE[0]) ; 

	hName[Form("h_Rz_st%d",ist)]->SetBinContent(nJ-1,R[1]) ; 
	hName[Form("h_Rz_st%d",ist)]->SetBinError(nJ-1,RE[1]) ; 
	
	hName[Form("h_Rtt_st%d",ist)]->SetBinContent(nJ-1,R[2]) ; 
	hName[Form("h_Rtt_st%d",ist)]->SetBinError(nJ-1,RE[2]) ; 
	
	//fill st binned 
	
	
	hName[Form("h_Rw_nJ%d",nJ)]->SetBinContent(ist+1,R[0]) ; 
	hName[Form("h_Rw_nJ%d",nJ)]->SetBinError(ist+1,RE[0]) ; 
	
	hName[Form("h_Rz_nJ%d",nJ)]->SetBinContent(ist+1,R[1]) ; 
	hName[Form("h_Rz_nJ%d",nJ)]->SetBinError(ist+1,RE[1]) ; 
	
	hName[Form("h_Rtt_nJ%d",nJ)]->SetBinContent(ist+1,R[2]) ; 
	hName[Form("h_Rtt_nJ%d",nJ)]->SetBinError(ist+1,RE[2]) ; 

	
}

void print_numbers(int nJ, float st1, float st2){

//	int nJ=2;
//	float st1=300;
//	float st2=600; 
	
	ofstream output; 
	output.open(Form("numbers_table_nJets%d_st%.0f-%.0f.tex",nJ,st1,st2));
	output << "\\begin{tabular}{";
	int N=5; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}" << endl;  
	output << "Control & Data & W+Jets & DY & TTBar & Other \\\\ \\hline" << endl; 
	for(int i=0; i<3; i++){
		output << i+1 << " & "; 
		output << Nevent[Form("N_%d_data_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " $\\pm$ "; 
		output << NeventE[Form("N_%d_data_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 
		
		output << Nevent[Form("N_%d_W_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " $\\pm$ ";
		output << NeventE[Form("N_%d_W_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 

		output << Nevent[Form("N_%d_DY_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " $\\pm$ "; 
		output << NeventE[Form("N_%d_DY_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 
		
		output << Nevent[Form("N_%d_ttbar_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " $\\pm$ "; 
		output << NeventE[Form("N_%d_ttbar_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " & "; 

		output << Nevent[Form("N_%d_other_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " $\\pm$ "; 
		output << NeventE[Form("N_%d_other_nJ%d_st%.0f-%.0f",i+1,nJ,st1,st2)] << " \\\\ " << endl; 

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
	/*
	for(int nJ=2; nJ<=4; nJ++){
		double st1=300;
		double st2=500;

		output << nJ << " & ";
		if(st2>0)output << st1 << "-" << st2 << " & "; 
		else output << "$\\ge$ " << st1 << " & "; 
		output << MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " \\\\ ";
		output << endl; 
		
		st1=500;
		st2=800;
		
		output << nJ << " & ";
		if(st2>0)output << st1 << "-" << st2 << " & "; 
		else output << "$\\ge$ " << st1 << " & ";  
		output << MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " \\\\ ";
		output << endl; 
		
		st1=800;
		st2=-1;
		
		output << nJ << " & ";
		if(st2>0)output << st1 << "-" << st2 << " & "; 
		else output << "$\\ge$ " << st1 << " & "; 
		
		output << MC_corr[Form("R_W_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_Z_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " & ";
		output << MC_corr[Form("R_ttbar_nJets%d_st%.0f-st%.0f",nJ,st1,st2)] << " \\\\ \\hline";
		output << endl; 
		
	}
	 */
	output << "\\hline" << endl; 
	output.close(); 
}

void get_column(TString name, TString file,double st1,double st2, double *N, double *NE){
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
	
	hist_name=name+"_W"+file; 
	double N1E=get_NE(st1,st2,hName[hist_name]); 
	hist_name=name+"_Z"+file;
	double N2E=get_NE(st1,st2,hName[hist_name]); 
	hist_name=name+"_1Mu_1El"+file; 
	double N3E=get_NE(st1,st2,hName[hist_name]); 
	NE[0]=N1E;//W enriched
	NE[1]=N2E; // Z enriched
	NE[2]=N3E; //ttbar enriched
}

double scaleStBins(TString hist_name,double st1, double st2, double Corr){
	double N=0;
	for(int ipt=1; ipt<=hName[hist_name]->GetNbinsX();ipt++){
		double st=hName[hist_name]->GetBinCenter(ipt); 

		//cout << "Corr: " << Corr << endl; 
		if(st>st1 && st<=st2)hName[hist_name]->SetBinContent(ipt, hName[hist_name]->GetBinContent(ipt)*Corr); 
		if(st>st1 && st<=st2)N+=hName[hist_name]->GetBinContent(ipt); 
	}
	return N; 
}

void bookHisto(){
	
	
	for(int ist=0; ist<NbinsCor; ist++){
		TString title=Form("%.0f < S_{T} < %.0f GeV",correction_bins[ist],correction_bins[ist+1]);
		CreateHistogram(Form("h_Rw_st%d",ist),title,"n-jets","R",3,1.5,4.5); 
		CreateHistogram(Form("h_Rz_st%d",ist),title,"n-jets","R",3,1.5,4.5);
		CreateHistogram(Form("h_Rtt_st%d",ist),title,"n-jets","R",3,1.5,4.5);
	}
	CreateHistogram("h_Rw_avg", " S_{T} averaged","n-jets","R",3,1.5,4.5); 
	CreateHistogram("h_Rz_avg", " S_{T} averaged","n-jets","R",3,1.5,4.5); 
	CreateHistogram("h_Rtt_avg", " S_{T} averaged","n-jets","R",3,1.5,4.5); 

	for(int nJ=2; nJ<=4; nJ++){
		CreateHistogram(Form("h_Rw_nJ%d",nJ),"","S_{T}","R",NbinsCor,correction_bins); 
		CreateHistogram(Form("h_Rz_nJ%d",nJ),"","S_{T}","R",NbinsCor,correction_bins);
		CreateHistogram(Form("h_Rtt_nJ%d",nJ),"","S_{T}","R",NbinsCor,correction_bins);
	}
	float st_bins[]={300,400,500,600,700,800,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	for(int mu=1; mu<=2; mu++){
		for(int nJ=2; nJ<=nJetsmax; nJ++){
			CreateHistogram(Form("st_nJets%d_%dMu_systematic",nJ,mu),"","S_{T}","systematic",N ,st_bins); 
			CreateHistogram(Form("st_nJets%d_%dMu_ratio_systematic",nJ,mu),"","S_{T}","systematic",N ,st_bins); 
			CreateHistogram(Form("expected_mass_nJets%d_%dMu",nJ,mu),"","Mass","Expected",8,250,1050); 
			CreateHistogram(Form("expected_sys_mass_nJets%d_%dMu",nJ,mu),"","Mass","Expected",8,250,1050); 
			CreateHistogram(Form("observed_mass_nJets%d_%dMu",nJ,mu),"","Mass","Expected",8,250,1050); 
			CreateHistogram(Form("signal_mass_nJets%d_%dMu",nJ,mu),"","Mass","Signal",8,250,1050); 

		}
		
	}
	
	CreateHistogram("Rw_pull","Rw_pull","Pull","Bins",25,-6,6);
	CreateHistogram("Rz_pull","Rz_pull","Pull","Bins",25,-6,6);
	CreateHistogram("Rtt_pull","Rtt_pull","Pull","Bins",25,-6,6);

	
	CreateHistogram("fit_ttbar_prob_pol1","Pol1", "Prob","Fits", 10,0,1); 
	CreateHistogram("fit_ttbar_prob_pol0","Pol0", "Prob","Fits",10,0,1); 

}

void fit_corrections(){

	double Rwavg[3];
	double Wwavg[3];
	double Swavg[3];
	
	double Rzavg[3];
	double Wzavg[3];
	double Szavg[3];

	
	double Rttavg[3];
	double Wttavg[3];
	double Sttavg[3];

	for(int nJ=0; nJ<3; nJ++){
		Rwavg[nJ]=0;
		Wwavg[nJ]=0;
		Swavg[nJ]=0;
		
		Rzavg[nJ]=0;
		Wzavg[nJ]=0;
		Szavg[nJ]=0;	
		
		Rttavg[nJ]=0;
		Wttavg[nJ]=0;
		Sttavg[nJ]=0;
	
	}
	
	
	for(int ist=0; ist<NbinsCor; ist++){
		TF1 *f = new TF1(Form("W_pol1_st%d",ist),"pol1",0.5,4.5); 
		f->SetLineColor(kRed);
		hName[Form("h_Rw_st%d",ist)]->Fit(Form("W_pol1_st%d",ist),"Q"); 
		f1Name[f->GetName()]=f; 
		
		
		TF1 *fz = new TF1(Form("Z_pol1_st%d",ist),"pol1",0.5,4.5);
		fz->SetLineColor(kBlue);
		hName[Form("h_Rz_st%d",ist)]->Fit(Form("Z_pol1_st%d",ist),"Q"); 
		f1Name[fz->GetName()]=fz; 
		
		TF1 *ftt0 = new TF1(Form("tt_pol0_st%d",ist),"pol0",0.5,4.5);
		ftt0->SetLineColor(419);
		ftt0->SetLineStyle(kDashed);
		hName[Form("h_Rtt_st%d",ist)]->Fit(Form("tt_pol0_st%d",ist),"Q0"); 
		f1Name[ftt0->GetName()]=ftt0;
		
		TF1 *ftt = new TF1(Form("tt_pol1_st%d",ist),"pol1",0.5,4.5);
		ftt->SetLineColor(419);
		hName[Form("h_Rtt_st%d",ist)]->Fit(Form("tt_pol1_st%d",ist),"Q"); 
		f1Name[ftt->GetName()]=ftt; 
		
	
		for(int nJ=0; nJ<3; nJ++){
			Rwavg[nJ]+=hName[Form("h_Rw_st%d",ist)]->GetBinContent(nJ+1)/hName[Form("h_Rw_st%d",ist)]->GetBinError(nJ+1);
			Wwavg[nJ]+=1/hName[Form("h_Rw_st%d",ist)]->GetBinError(nJ+1);
			Swavg[nJ]+=TMath::Power(hName[Form("h_Rw_st%d",ist)]->GetBinError(nJ+1)/hName[Form("h_Rw_st%d",ist)]->GetBinContent(nJ+1),2);

			Rzavg[nJ]+=hName[Form("h_Rz_st%d",ist)]->GetBinContent(nJ+1)/hName[Form("h_Rz_st%d",ist)]->GetBinError(nJ+1);
			Wzavg[nJ]+=1/hName[Form("h_Rz_st%d",ist)]->GetBinError(nJ+1);
			Szavg[nJ]+=TMath::Power(hName[Form("h_Rz_st%d",ist)]->GetBinError(nJ+1)/hName[Form("h_Rz_st%d",ist)]->GetBinContent(nJ+1),2);
			
			Rttavg[nJ]+=hName[Form("h_Rtt_st%d",ist)]->GetBinContent(nJ+1)/hName[Form("h_Rtt_st%d",ist)]->GetBinError(nJ+1);
			Wttavg[nJ]+=1/hName[Form("h_Rtt_st%d",ist)]->GetBinError(nJ+1);
			Sttavg[nJ]+=TMath::Power(hName[Form("h_Rtt_st%d",ist)]->GetBinError(nJ+1)/hName[Form("h_Rtt_st%d",ist)]->GetBinContent(nJ+1),2);

			
			hName["Rw_pull"]->Fill((1-hName[Form("h_Rw_st%d",ist)]->GetBinContent(nJ+1))/hName[Form("h_Rw_st%d",ist)]->GetBinError(nJ+1)); 
			hName["Rz_pull"]->Fill((1-hName[Form("h_Rz_st%d",ist)]->GetBinContent(nJ+1))/hName[Form("h_Rz_st%d",ist)]->GetBinError(nJ+1)); 
			hName["Rtt_pull"]->Fill((1-hName[Form("h_Rtt_st%d",ist)]->GetBinContent(nJ+1))/hName[Form("h_Rtt_st%d",ist)]->GetBinError(nJ+1)); 

			
		}

		
		hName["fit_ttbar_prob_pol1"]->Fill(TMath::Prob(ftt->GetChisquare(),ftt->GetNDF())); 
		hName["fit_ttbar_prob_pol0"]->Fill(TMath::Prob(ftt0->GetChisquare(),ftt0->GetNDF())); 

		cout << "pol1 Prob: " << TMath::Prob(ftt->GetChisquare(),ftt->GetNDF()) << endl; 
		cout << "pol0 Prob: " << TMath::Prob(ftt0->GetChisquare(),ftt0->GetNDF()) << endl;
		
	}
	
	for(int nJ=0; nJ<3;nJ++){
		cout << "nJ: " << nJ+1 << endl; 
		hName["h_Rw_avg"]->SetBinContent(nJ+1,Rwavg[nJ]/Wwavg[nJ]); 
		hName["h_Rz_avg"]->SetBinContent(nJ+1,Rzavg[nJ]/Wzavg[nJ]); 
		hName["h_Rtt_avg"]->SetBinContent(nJ+1,Rttavg[nJ]/Wttavg[nJ]); 
		
		Swavg[nJ]=TMath::Sqrt(Swavg[nJ])*hName["h_Rw_avg"]->GetBinContent(nJ+1);
		Szavg[nJ]=TMath::Sqrt(Szavg[nJ])*hName["h_Rz_avg"]->GetBinContent(nJ+1);
		Sttavg[nJ]=TMath::Sqrt(Sttavg[nJ])*hName["h_Rtt_avg"]->GetBinContent(nJ+1);
		
		cout << "W unc " << Swavg[nJ] << " z unc: " << Szavg[nJ] << " tt unc " << Sttavg[nJ] << endl; 
		
		hName["h_Rw_avg"]->SetBinError(nJ+1,Swavg[nJ]);
		hName["h_Rz_avg"]->SetBinError(nJ+1,Szavg[nJ]);
		hName["h_Rtt_avg"]->SetBinError(nJ+1,Sttavg[nJ]);
		
	}
}	

void draw_corrections(){
	
	CreateCanvas("Correction_Pull","",600,600);
	CName["Correction_Pull"]->cd();
	hName["Rw_pull"]->SetFillColor(kRed);
	hName["Rz_pull"]->SetFillColor(kBlue);
	hName["Rtt_pull"]->SetFillColor(419);
	
	hName["Rw_pull"]->Draw("histo");
	hName["Rz_pull"]->Draw("histo same");
	hName["Rtt_pull"]->Draw("histo same");

	
	
	CreateCanvas("Correction_Plot", "",900,600); 
	CName["Correction_Plot"]->Divide(3,2,0.01,0.01); 
	
	for(int ist=0; ist<NbinsCor; ist++){
		hName[Form("h_Rw_st%d",ist)]->SetMarkerStyle(8); 
		hName[Form("h_Rz_st%d",ist)]->SetMarkerStyle(22);
		hName[Form("h_Rtt_st%d",ist)]->SetMarkerStyle(23); 
		
		hName[Form("h_Rw_st%d",ist)]->SetMarkerSize(1.5);
		hName[Form("h_Rz_st%d",ist)]->SetMarkerSize(1.5); 
		hName[Form("h_Rtt_st%d",ist)]->SetMarkerSize(1.5); 
		
		hName[Form("h_Rw_st%d",ist)]->SetMarkerColor(kRed); 
		hName[Form("h_Rz_st%d",ist)]->SetMarkerColor(kBlue); 
		hName[Form("h_Rtt_st%d",ist)]->SetMarkerColor(419); 
	
		hName[Form("h_Rw_st%d",ist)]->SetLineColor(kRed); 
		hName[Form("h_Rz_st%d",ist)]->SetLineColor(kBlue); 
		hName[Form("h_Rtt_st%d",ist)]->SetLineColor(419); 
	}	

	
	for(int nJ=2; nJ<=4; nJ++){
		hName[Form("h_Rw_nJ%d",nJ)]->SetMarkerStyle(8); 
		hName[Form("h_Rz_nJ%d",nJ)]->SetMarkerStyle(22);
		hName[Form("h_Rtt_nJ%d",nJ)]->SetMarkerStyle(23); 
		
		hName[Form("h_Rw_nJ%d",nJ)]->SetMarkerSize(1.5);
		hName[Form("h_Rz_nJ%d",nJ)]->SetMarkerSize(1.5); 
		hName[Form("h_Rtt_nJ%d",nJ)]->SetMarkerSize(1.5); 
		
		hName[Form("h_Rw_nJ%d",nJ)]->SetMarkerColor(kRed); 
		hName[Form("h_Rz_nJ%d",nJ)]->SetMarkerColor(kBlue); 
		hName[Form("h_Rtt_nJ%d",nJ)]->SetMarkerColor(419); 
		
		hName[Form("h_Rw_nJ%d",nJ)]->SetLineColor(kRed); 
		hName[Form("h_Rz_nJ%d",nJ)]->SetLineColor(kBlue); 
		hName[Form("h_Rtt_nJ%d",nJ)]->SetLineColor(419);
	}
	
	TLatex txt;
	
	TLegend Leg(0.44,0.16,0.84,0.36);
	Leg.SetFillColor(10);
	Leg.SetLineColor(10);
	Leg.SetLineWidth(0);
	Leg.SetBorderSize(0);
	Leg.SetTextSize(0.05); 
	
	Leg.AddEntry(hName[Form("h_Rw_st%d",2)],"W+Jets","P"); 
	Leg.AddEntry(hName[Form("h_Rz_st%d",2)],"Z+Jets","P"); 
	Leg.AddEntry(hName[Form("h_Rtt_st%d",2)],"tt+Jets","P"); 

	
	for(int ist=0; ist<NbinsCor; ist++){
		CName["Correction_Plot"]->cd(ist+1); 
		txt.SetNDC(kTRUE); 
		hName[Form("h_Rw_st%d",ist)]->SetMinimum(0.5);
		hName[Form("h_Rw_st%d",ist)]->SetMaximum(1.5); 
		txt.DrawLatex(0.3,0.8,Form("%.0f < S_{T} < %.0f",correction_bins[ist],correction_bins[ist+1])); 
		//hName[Form("h_Rw_st%d",ist)]->SetTitle(""); 
		hName[Form("h_Rw_st%d",ist)]->Draw(" E1X0"); 
		hName[Form("h_Rz_st%d",ist)]->Draw("  E1X0 same"); 
		hName[Form("h_Rtt_st%d",ist)]->Draw("  E1X0 same"); 
		
		
	}
	CName["Correction_Plot"]->cd(6);
	hName["h_Rw_avg"]->SetMarkerStyle(8);
	hName["h_Rz_avg"]->SetMarkerStyle(22);
	hName["h_Rtt_avg"]->SetMarkerStyle(23);
	
	hName["h_Rw_avg"]->SetMarkerSize(1.5);
	hName["h_Rz_avg"]->SetMarkerSize(1.5);
	hName["h_Rtt_avg"]->SetMarkerSize(1.5);
	
	hName["h_Rw_avg"]->SetMarkerColor(kRed);
	hName["h_Rz_avg"]->SetMarkerColor(kBlue);
	hName["h_Rtt_avg"]->SetMarkerColor(419);
	
	hName["h_Rw_avg"]->SetLineColor(kRed);
	hName["h_Rz_avg"]->SetLineColor(kBlue);
	hName["h_Rtt_avg"]->SetLineColor(419);
	
	hName["h_Rw_avg"]->SetMinimum(0.5);
	hName["h_Rw_avg"]->SetMaximum(1.5); 
	hName["h_Rw_avg"]->Draw(" E1X0"); 
	hName["h_Rz_avg"]->Draw("  E1X0 same"); 
	hName["h_Rtt_avg"]->Draw("  E1X0 same");
	Leg.DrawClone(); 
	
	
	CreateCanvas("Correction_Plot_st", "",900,600); 
	CName["Correction_Plot_st"]->Divide(3,1,0.01,0.01); 

	Leg.SetX1NDC(0.4);
	Leg.SetX2NDC(0.87);
	Leg.SetY1NDC(0.7);
	Leg.SetY2NDC(0.9);
	
	for(int nJ=2; nJ<=4; nJ++){
		CName["Correction_Plot_st"]->cd(nJ-1); 
		txt.SetNDC(kTRUE); 
		hName[Form("h_Rw_nJ%d",nJ)]->SetMinimum(0.8);
		hName[Form("h_Rw_nJ%d",nJ)]->SetMaximum(1.2); 
		
		hName[Form("h_Rw_nJ%d",nJ)]->Draw(" E1X0"); 
		hName[Form("h_Rz_nJ%d",nJ)]->Draw("  E1X0 same"); 
		hName[Form("h_Rtt_nJ%d",nJ)]->Draw("  E1X0 same"); 
		
		if(nJ==2) txt.DrawLatex(0.3,0.8,"2-jets"); 
		if(nJ==3) txt.DrawLatex(0.3,0.8,"3-jets"); 
		if(nJ==4) txt.DrawLatex(0.3,0.8,"4-jets"); 
		
	}
	CName["Correction_Plot_st"]->cd(1); 
	Leg.DrawClone(); 
								 
}

void rescaleMC(int nJ, int mu){

	
	TString variable=Form("st_nJets%d_%dMu",nJ,mu);
	TString file="_wJets"; 
//	int bkgnJ=nJ-1; 
//	if(mu==2) bkgnJ=nJ;
	
	double CW[NbinsCor]; 
	double CZ[NbinsCor]; 
	double Ctt[NbinsCor];

	
	int bkgnJ=nJ; 
	double sys=0; 
	for(int ist=0; ist<NbinsCor; ist++){
		float st1=correction_bins[ist];
		float st2=correction_bins[ist+1];
		file="_wJets"; 
		double Corr=0; 
				
		if(nJ<4) Corr=hName[Form("h_Rw_st%d",ist)]->GetBinContent(nJ-1); 
		else Corr=f1Name[Form("W_pol1_st%d",ist)]->Eval(nJ);
		scaleStBins(variable+file,st1,st2,Corr);
		CW[ist]=Corr; 
		
		file="_dy";
		if(nJ<4) Corr=hName[Form("h_Rz_st%d",ist)]->GetBinContent(nJ-1); 
		else Corr=f1Name[Form("Z_pol1_st%d",ist)]->Eval(nJ);
		scaleStBins(variable+file,st1,st2,Corr);

		CZ[ist]=Corr; 
		
		file="_ttbar";
		if(nJ<4) Corr=hName[Form("h_Rtt_st%d",ist)]->GetBinContent(nJ-1); 
		else Corr=f1Name[Form("tt_pol0_st%d",ist)]->Eval(nJ);
		scaleStBins(variable+file,st1,st2,Corr);

		Ctt[ist]=Corr; 

		
	}
	
	for(int ist=1; ist<=hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->GetNbinsX();ist++){
		double sys=0; 
		file="_wJets"; 
		double c=0;
		
		if(ist<NbinsCor)c=CW[ist-1];
		else c=CW[NbinsCor-1];
		double Nw=hName[variable+file]->GetBinContent(ist); 
		sys+=TMath::Power(Nw*(1-c),2);
		
		file="_dy"; 
		if(ist<NbinsCor)c=CZ[ist-1];
		else c=CZ[NbinsCor-1];
		double Nz=hName[variable+file]->GetBinContent(ist); 
		sys+=TMath::Power(Nz*(1-c),2);
		
		file="_ttbar";
		if(ist<NbinsCor)c=Ctt[ist-1];
		else c=Ctt[NbinsCor-1];
		double Nttbar=hName[variable+file]->GetBinContent(ist); 
		sys+=TMath::Power(Nttbar*(1-c),2);
		sys=TMath::Sqrt(sys); 
		
		hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->SetBinContent(ist,sys);
	}
	
}

void divide_hist(TH1F *data, TH1F *mc, double *x, int N){

	//cout << "N: " << N << endl; 
	TString name=data->GetName();
	name=name+"_dataMC";
	TH1F *rebindata=(TH1F*)data->Rebin(N,name,x);
	TH1F *rebinMC=(TH1F*)mc->Rebin(N,"RebinMC",x);
	
	rebindata->Divide(rebinMC); 
	rebindata->GetYaxis()->SetTitle("Data/MC"); 
	rebindata->SetMinimum(0.5);
	rebindata->SetMaximum(1.5);
	
	hName[rebindata->GetName()]=rebindata; 
	delete rebinMC; 
	
}

void add_corrected_bkg(int nJ, int mu){
	TString name=Form("st_nJets%d_%dMu",nJ,mu);
	TString file="_QCD"; 
	
	TH1F *h=(TH1F*)hName[name+file]->Clone(Form("st_nJets%d_%dMu_correctedMC",nJ,mu));
	hName[h->GetName()]=h; 
	
	file="_diboson"; 
	hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]->Add(hName[name+file]); 
	
	file="_singleTop"; 
	hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]->Add(hName[name+file]); 
	
	file="_ttbar"; 
	hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]->Add(hName[name+file]); 
	
	file="_wJets"; 
	hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]->Add(hName[name+file]); 
	
	file="_dy"; 
	hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]->Add(hName[name+file]); 
	
}

void relative_systematic(int nJ, int mu){
	cout << "systematic: " << endl; 
	for(int ist=1; ist<=hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->GetNbinsX();ist++){
		hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->SetBinError(ist,0); 
	}
	cout << "divide: " << endl; 
	hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->Divide(hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]); 
	cout << "scale to ratio: " << endl; 
	for(int ist=1; ist<=hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->GetNbinsX();ist++){
		double R=hName[Form("st_nJets%d_%dMu_singleMu",nJ,mu)]->GetBinContent(ist)/hName[Form("st_nJets%d_%dMu_correctedMC",nJ,mu)]->GetBinContent(ist);
		double Rsys=hName[Form("st_nJets%d_%dMu_systematic",nJ,mu)]->GetBinContent(ist)*R;
		
		cout << "R: " << R << " Rsys: "<< Rsys << endl;
		
		hName[Form("st_nJets%d_%dMu_ratio_systematic",nJ,mu)]->SetBinContent(ist, 1); 
		hName[Form("st_nJets%d_%dMu_ratio_systematic",nJ,mu)]->SetBinError(ist, Rsys); 

	}
	TGraphAsymmErrors *gr=new TGraphAsymmErrors(hName[Form("st_nJets%d_%dMu_ratio_systematic",nJ,mu)]);
	gr->SetName(Form("st_nJets%d_%dMu_ratio_systematic",nJ,mu)); 
	for(int i=0; i<gr->GetN(); i++){
		double ex=hName[Form("st_nJets%d_%dMu_ratio_systematic",nJ,mu)]->GetBinWidth(i+1)/2; 

		gr->SetPointEXhigh(i,ex);
		gr->SetPointEXlow(i,ex); 
		
	}
	grName[gr->GetName()]=gr; 
				
	
}

void fill_stack(TString variable, bool blind){
	//modify this code to be given a string with the name
	//create just one stack per call
	//fix y axis min thing
	//plot h_mumu in bins
//	cout << "variable: " << variable << endl; 

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
	hName[name+file]->SetMarkerColor(419); 
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
	
	shape_pad->SetTopMargin(0.1);
	shape_pad->SetBottomMargin(0.05);
	ratio_pad->SetTopMargin(0.05);
	ratio_pad->SetBottomMargin(0.3); 
	
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
		hName[name+file]->GetYaxis()->SetTitleOffset(0.65);
		hName[name+file]->GetYaxis()->SetTitleSize(0.08);
		hName[name+file]->GetYaxis()->SetLabelSize(0.05);
		hName[name+file]->GetXaxis()->SetLabelSize(0);
		hName[name+file]->GetXaxis()->SetTitleSize(0);		
		
		hName[name+file]->DrawCopy(" E1 ");
		stackName[jet_stack]->DrawClone("hist same");
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
	double xBins[]={0,30,50,100,150,200,300,400,1500};
	double stBins[]={0,300,400,500,600,800,1000,3000};
	
	int Nx=sizeof(xBins)/sizeof(double)-1;
	int Nst=sizeof(stBins)/sizeof(double)-1;
	
	if(name.Contains("jetpT") || name.Contains("met") || name.Contains("muonpT")){
	
		divide_hist( hName[name+file], hName[tot_MC],xBins, Nx);
	
		TString dataMC=name+file;
		dataMC=dataMC+"_dataMC"; 
		hName[dataMC]->GetYaxis()->SetTitleOffset(0.45);
		hName[dataMC]->GetYaxis()->SetTitleSize(0.15);
		hName[dataMC]->GetYaxis()->SetTitleOffset(0.45);
		hName[dataMC]->GetYaxis()->SetLabelSize(0.08);
		
		hName[dataMC]->GetXaxis()->SetTitleSize(0.15); 
		hName[dataMC]->GetXaxis()->SetLabelSize(0.12); 
		hName[dataMC]->GetXaxis()->SetTitleSize(0.15); 
		
		hName[dataMC]->DrawCopy("histo E1"); 
	
	}
	else{
		if(name.Contains("st_nJets2_1Mu")==1){
		cout << "bin center: " << hName[name+file]->GetBinCenter(10) << endl; 
		cout << "Ndata: " << hName[name+file]->GetBinContent(10) << " +/- " << hName[name+file]->GetBinError(10) << endl; 
		cout << "Nmc: " << hName[tot_MC]->GetBinContent(10) << " +/- "  << hName[tot_MC]->GetBinError(10) << endl; 
		}
		TString dataMC=name+file;
		dataMC=dataMC+"_dataMC"; 
		TH1F *R=(TH1F*)hName[name+file]->Clone(dataMC);
		R->Divide(hName[tot_MC]); 
		hName[dataMC]=R; 
			
		hName[dataMC]->GetYaxis()->SetTitle("Data/MC"); 
		hName[dataMC]->SetMinimum(0.5);
		hName[dataMC]->SetMaximum(1.5);
		
		if(name.Contains("st_nJets2_1Mu")==1)cout	<< "R: " << hName[name+file]->GetBinContent(10) << " +/- " <<  hName[name+file]->GetBinError(10) << endl; 

		hName[dataMC]->GetYaxis()->SetTitleOffset(0.45);
		hName[dataMC]->GetYaxis()->SetTitleSize(0.15);
		hName[dataMC]->GetYaxis()->SetTitleOffset(0.45);
		hName[dataMC]->GetYaxis()->SetLabelSize(0.08);

		
		hName[dataMC]->GetXaxis()->SetTitleSize(0.15); 
		hName[dataMC]->GetXaxis()->SetLabelSize(0.12); 
		hName[dataMC]->GetXaxis()->SetTitleSize(0.15); 
		
		if(!blind)hName[name+file]->DrawCopy("histo E1"); 
		if(name.Contains("st_nJets5_1Mu") || name.Contains("st_nJets5_2Mu") || name.Contains("st_nJets6_1Mu") || name.Contains("st_nJets6_2Mu") || name.Contains("st_nJets7_1Mu") || name.Contains("st_nJets7_2Mu")){
			TString sysType="_ratio_systematic";
			if(name.Contains("1El")==0)grName[name+sysType]->SetFillColor(38);
			if(name.Contains("1El")==0)grName[name+sysType]->SetFillStyle(3244);

			if(name.Contains("1El")==0)grName[name+sysType]->Draw("E2 same");
		}
		
	}

	
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