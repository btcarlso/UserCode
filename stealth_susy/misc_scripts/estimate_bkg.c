/*
 *  estimate_bkg.c
 *  
 *
 *  Created by Benjamin Carlson on 10/30/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "estimate_bkg.h"

void estimate_bkg(){
	gROOT->SetBatch();
	open_files(); 
	get_hists();
	get_ratios();
	CreateHistogram("nJet_exp_1Mu","Expected","n-jets","Events",8,0.5,8.5); 
	CreateHistogram("nJet_exp_2Mu","Expected","n-jets","Events",8,0.5,8.5); 

	bkg_est(4,1,"singleMu");
	bkg_est(5,1,"singleMu");
	
	bkg_est(4,2,"singleMu");
	bkg_est(5,2,"singleMu");
	bkg_est(6,2,"singleMu");
	bkg_est(7,2,"singleMu");
	
	for(int nJ=4; nJ<=7; nJ++)bkg_est(nJ,2,"allMC");
	
	
	//fit(1);
	//fit(2); 
	
	TFile *output_file = new TFile("plots_bkg_exp.root","RECREATE"); 
	output_file->cd();
	bool print=1; 
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		if(print)CName[it->first]->Print("/uscms_data/d3/btcarlso/Figures_SUSYAN/"+it->first+".png");
		CName[it->first]->Write();
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

double power_law(double *x, double *par){
	double x_=x[0]/8000; 
	return  par[0]/(TMath::Power(x_,par[1])); 
}

double pol_scaled(double *x, double *par){
	return par[0]+par[1]*x[0]/8000; 
}

double expo_quad(double *x, double *par){
	double x_=x[0]/8000;
	return par[0]*TMath::Exp(par[1]*x_+par[2]*x_*x_);
}


double expo_normalized(double *x, double *par){
	double x_=x[0]/8000;
	return par[0]*TMath::Exp(par[1]*x_);
}

double pareto(double *x, double *par){
	double alpha=par[0]; 
	double beta=par[1];
	double mu=par[2];
	
	double power=-(1+beta)/beta;
	
	return (1/alpha)*TMath::Power(1+beta*(x[0]-mu)/alpha,power); 
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
	for(int num=3; num<=7; num++){
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
		for(int num=3; num<=7; num++){
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

void fit(int mu){
	cout << "mu: " << mu << endl; 
	double st1=500; 
	double st2=3000; 
	TGraphErrors *gr = new TGraphErrors(7); 
	cout << "fit exponential: " << endl; 
	TF1 *expo=new TF1("expo",expo_normalized,st1,st2,2);

	int ip=0; 
	
	for(int nJ=2; nJ<=7; nJ++){
		CreateCanvas(Form("Fits_nJets%d_%dMu",nJ,mu),"Fit",Cx,Cy); 
		CName[Form("Fits_nJets%d_%dMu",nJ,mu)]->cd(); 
		hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->SetAxisRange(500,3000); 
		gPad->SetLogy();
		hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Draw("E1"); 
		double A0=hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Integral();
		cout << "initial A: " << A0 << endl; 
		expo->SetParameter(0,A0); 
		expo->SetParameter(1,-10); 
		
		hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Fit("expo","QR0"); 
		expo->DrawClone("same"); 
		cout << "nJets: " << nJ << " slope: " << expo->GetParameter(1) << " +/- " << expo->GetParError(1) << endl; 
		gr->SetPoint(ip,nJ,expo->GetParameter(1)); 
		gr->SetPointError(ip,0.,expo->GetParError(1)); 
		ip++;
	}
	cout << "fit exponential with quadratic term: " << endl; 
	TF1 *expo_q=new TF1("expo_quad",expo_quad,st1,st2,3); 
	
	TF1 *pareto_fcn = new TF1("pareto",pareto,st1,st2,3); 
	pareto_fcn->SetParNames("alpha","beta","mu"); 
	pareto_fcn->SetLineColor(kRed);
	pareto_fcn->SetLineStyle(kDashed);

	
	expo_q->SetLineColor(kBlue); 
	for(int nJ=2; nJ<=7; nJ++){
		expo_q->SetParameter(0,hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Integral()); 
		expo_q->SetParameter(1,-100);
		expo_q->SetParameter(2,70);
		CName[Form("Fits_nJets%d_%dMu",nJ,mu)]->cd(); 
		hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Fit("expo_quad","QR0"); 
		hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Fit("expo_quad","QR0"); 

		expo_q->DrawClone("same"); 
		
		pareto_fcn->SetParLimits(0,1,100000);
		pareto_fcn->SetParameter(0,400);
		pareto_fcn->FixParameter(2,500);
		pareto_fcn->SetParameter(1,0.05);
		hName[Form("st_nJets%d_%dMu_50GeVallMC",nJ,mu)]->Fit("pareto","WR0");		

		pareto_fcn->DrawClone("same"); 
		
		cout << "nJets: " << nJ << " A: " << expo_q->GetParameter(0) << " +/- " << expo_q->GetParError(0) << endl; 
		cout << "nJets: " << nJ << " slope: " << expo_q->GetParameter(1) << " +/- " << expo_q->GetParError(1) << endl; 
		cout << "nJets: " << nJ << " curvature: " << expo_q->GetParameter(2) << " +/- " << expo_q->GetParError(2) << endl; 
	}
	
	CreateCanvas(Form("expo_slope_%dMu",mu),"",Cx,Cy); 
	CName[Form("expo_slope_%dMu",mu)]->cd(); 
	gr->DrawClone("ap"); 
	delete gr; 

}

void bkg_est(int nJ, int mu, string file){
	cout << "nJ: " << nJ << " mu: "<< mu  << " " << file << endl; 
	int bkg_nJ=3; 
	int bkg_nJ_sys=4;
	int ratio_num=nJ; 
	int ratio_den=nJ-1;
	if(nJ>=7) {
		ratio_num=7; 
		ratio_den=5; 
	}
	cout << "Get ratios and create exp. histograms: " << endl; 
	float st_bins[]={0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1200,1500,2000,3000}; 
	int N=sizeof(st_bins)/sizeof(float)-1;
	
	TString ratio_name=Form("Ratio%d-%d_%dMuallMC",ratio_num,bkg_nJ,mu); //Get ratio correction from MC 
	TString ratio_name_sys=Form("Ratio%d-%d_%dMuallMC",ratio_num,bkg_nJ_sys,mu); //Get ratio correction from MC 

	CreateHistogram(Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str()),"","S_{T}","Expected Events",N,st_bins); 
	CreateHistogram(Form("stexpEp_nJets%d_%dMu%s",nJ,mu,file.c_str()),"","S_{T}","Expected Events",N,st_bins);
	CreateHistogram(Form("stexpEm_nJets%d_%dMu%s",nJ,mu,file.c_str()),"","S_{T}","Expected Events",N,st_bins);

//	double st_threshold[]={500,500,875,1100,1300,1300,2000};
	
	cout << "Fill expected histogram: " << endl; 
	
	for(int ist=1; ist<=hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetNbinsX(); ist++){
		double st=hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinCenter(ist);
		double N=hName[Form("st_nJets%d_%dMu%s",bkg_nJ,mu,file.c_str())]->GetBinContent(ist); 
		double NE=hName[Form("st_nJets%d_%dMu%s",bkg_nJ,mu,file.c_str())]->GetBinError(ist); 
		int iB=hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->FindBin(st)-1; 
		double Ncor=grName[ratio_name]->GetY()[iB]*N; 
		//cout << "iB: " << iB << endl; 
		
		double fracEp=0; 
		double fracEm=0; 
		
		fracEp+=TMath::Power((grName[ratio_name]->GetErrorYhigh(iB))/grName[ratio_name]->GetY()[iB],2); 
		fracEm+=TMath::Power((grName[ratio_name]->GetErrorYlow(iB))/grName[ratio_name]->GetY()[iB],2); 

		if(N>0){
			fracEp+=TMath::Power(NE/N,2); 
			fracEm+=TMath::Power(NE/N,2); 
		}
		
		fracEp=TMath::Sqrt(fracEp); 
		fracEm=TMath::Sqrt(fracEm);
		
		double NcorEp=(1+fracEp)*Ncor; 
		double NcorEm=(1-fracEm)*Ncor; 
		
		
		//cout << "stat: " << (grName[ratio_name]->GetErrorYhigh(iB))/grName[ratio_name]->Eval(st)<< " sys: " << (grName[ratio_name]->Eval(st)-grName[ratio_name_sys]->Eval(st))/grName[ratio_name]->Eval(st) << endl;
		//cout << "st: " << st << " frac+: " << fracEp << " frac-: " << fracEm << endl; 
		if(st>=500){
			//cout << "st: " << st << " NE/N: " << NE/N << endl; 
			//cout << "st: " << st << " Ncor : " << Ncor << " Ncor+: " << NcorEp << " Ncor-: " << NcorEm << endl; 
			
			hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetBinContent(ist,Ncor); 
			hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetBinError(ist,0);
			hName[Form("stexpEp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetBinContent(ist, NcorEp); 
			hName[Form("stexpEm_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetBinContent(ist, NcorEm);
			
			hName[Form("stexpEp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetBinError(ist,0); 
			hName[Form("stexpEm_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetBinError(ist,0);
		}

	}
	
	cout << "Normalized 3 jet exp. " << endl; 
	
	int B1=hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->FindBin(600);
	int B2=hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->FindBin(800);
	
	double SFnocor=hName[Form("st_nJets%d_%dMu_50GeV%s",nJ,mu,file.c_str())]->Integral(B1,B2)/hName[Form("st_nJets%d_%dMu_50GeV%s",bkg_nJ,mu,file.c_str())]->Integral(B1,B2);
	
	cout << "Draw: " << endl; 
	
	CreateCanvas(Form("bkg_expected_nJets%d_%dMu%s",nJ,mu,file.c_str()),"",Cx,Cy); 
	
	TPad *shape_pad = new TPad(Form("stpad_%d",nJ),"st-shape",0.0,0.3,1,1);
	TPad *ratio_pad = new TPad(Form("ratiopad_%d",nJ),"Ratio",0,0.,1,0.3);
	CName[Form("bkg_expected_nJets%d_%dMu%s",nJ,mu,file.c_str())]->cd();

	shape_pad->Draw();
	ratio_pad->Draw();
	
	shape_pad->cd(); 
	gPad->SetLogy();
	hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetLineColor(kRed); 
	hName[Form("stexpEm_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetLineColor(kOrange); 
	hName[Form("stexpEp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetLineColor(kOrange); 

	
	hName[Form("st_nJets%d_%dMu_50GeV%s",bkg_nJ,mu,file.c_str())]->SetLineColor(kGreen); 
	TH1F *h=(TH1F*)hName[Form("st_nJets%d_%dMu_50GeV%s",bkg_nJ,mu,file.c_str())]->Clone(Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())); 
	h->Scale(SFnocor); 
	
	TH1F *ratio=(TH1F*)hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->Clone(Form("Obs/exp_nJets%d_%dMu%s",nJ,mu,file.c_str())); 
	
	ratio->Divide(hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]); 
	ratio->SetAxisRange(500,3000);
	
	TGraphAsymmErrors *grRatio=new TGraphAsymmErrors(ratio); 
	TGraphAsymmErrors *grRatio_stat = new TGraphAsymmErrors(ratio); 
	
	for(int i=1; i<=ratio->GetNbinsX(); i++){
		double st=hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinCenter(i);
		
		double sigma_data=hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinError(i); 
		double data=hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinContent(i); 

		double sigma_predEp=hName[Form("stexpEp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinContent(i); 
		double sigma_predEm=hName[Form("stexpEm_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinContent(i); 

		double pred=hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->GetBinContent(i); 
		
		
		double sigmaEp=0; 
		double sigmaEm=0; 
		
		sigmaEp+=TMath::Power(sigma_data/data,2); 
		sigmaEp+=TMath::Power(sigma_predEp/pred,2); 
		
		sigmaEm+=TMath::Power(sigma_data/data,2); 
		sigmaEm+=TMath::Power(sigma_predEm/pred,2); 
	
		sigmaEp=TMath::Sqrt(sigmaEp); 
		sigmaEm=TMath::Sqrt(sigmaEm); 
		
		cout << "sigmaEp: " << sigmaEp << " stat: " << sigma_data/data << endl; 
		
		double R=ratio->GetBinContent(i);
		grRatio->SetPoint(i-1,st,R); 
		grRatio_stat->SetPoint(i-1,st,R); 
		grRatio->SetPointError(i-1,0,0,R*sigmaEp,R*sigmaEm); 
		grRatio_stat->SetPointError(i-1,0,0,R*sigma_data/data,R*sigma_data/data);

	}


	hName[Form("st_nJets%d_%dMuallMC",nJ,mu,file.c_str())]->SetFillColor(kBlue); 
	hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetMarkerStyle(22); 
	hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->SetLineColor(kBlack); 
	
	hName[Form("st_nJets%d_%dMuallMC",nJ,mu,file.c_str())]->SetAxisRange(500,2000);
	hName[Form("st_nJets%d_%dMuallMC",nJ,mu,file.c_str())]->DrawCopy("histo"); 
	hName[Form("st_nJets%d_%dMu%s",nJ,mu,file.c_str())]->DrawCopy("E1 same"); 
	hName[Form("stexp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->Draw("same");
	hName[Form("stexpEp_nJets%d_%dMu%s",nJ,mu,file.c_str())]->DrawCopy("same");
	hName[Form("stexpEm_nJets%d_%dMu%s",nJ,mu,file.c_str())]->DrawCopy("same");

	//h->DrawCopy("same");
	
	TF1 *line= new TF1("line","pol0",500,3000); 
	line->SetParameter(0,1); 
	line->SetLineColor(kRed); 
	
	ratio_pad->cd();
	grRatio->SetMarkerStyle(7);
	grRatio->SetMarkerSize(1); 
	grRatio_stat->SetLineColor(kRed); 
	grRatio->Draw("ap"); 
	grRatio->GetXaxis()->SetRangeUser(500,2000); 
	grRatio->DrawClone("ap"); 
	grRatio_stat->DrawClone("|| same"); 
	line->Draw("same"); 
	
	delete h; 
}

void get_hists(){
	TString file="allMC"; 
	for(int mu=1; mu<=2; mu++){
		for(int nJ=2; nJ<=7; nJ++){
			TString name=Form("st_nJets%d_%dMu_uW_50GeV",nJ,mu); 
			GetHistogram(name,file);
			name=Form("st_nJets%d_%dMu_uW",nJ,mu); 
			GetHistogram(name,file);
			name=Form("st_nJets%d_%dMu",nJ,mu); 
			GetHistogram(name,file);
			name=Form("st_nJets%d_%dMu_50GeV",nJ,mu); 
			GetHistogram(name,file);
		}
	}
	
	 file="singleMu"; 
	for(int mu=1; mu<=2; mu++){
		for(int nJ=2; nJ<=7; nJ++){
			TString name=Form("st_nJets%d_%dMu_uW_50GeV",nJ,mu); 
			GetHistogram(name,file);
			name=Form("st_nJets%d_%dMu_uW",nJ,mu); 
			GetHistogram(name,file);
			name=Form("st_nJets%d_%dMu",nJ,mu); 
			GetHistogram(name,file);
			name=Form("st_nJets%d_%dMu_50GeV",nJ,mu); 
			GetHistogram(name,file);
		}
	}
	
	
}

void GetGraph(TString graph_name, TString file){
	TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fName[file]->FindObjectAny(graph_name);
	grName[graph_name+file]=gr; 
}

void GetHistogram(TString hist_name, TString file){
	//cout << "file: " << file << " name: " << hist_name << endl; 
	TH1F *h=(TH1F*)fName[file]->FindObjectAny(hist_name);
	h->SetStats(kFALSE);
	hName[hist_name+file]=h; 
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
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	hName[name] = h;
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	createC->SetFillColor(10); 
	CName[Name]=createC;
}

