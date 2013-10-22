/*
 *  Generate_LS.c
 *  
 *
 *  Created by Benjamin Carlson on 7/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "Generate_LS.h"

using namespace std;

void Generate_LS(int job){
	gROOT->SetBatch();
	int iy=int(job)/fNpt;
	int ipt=int(job)%fNpt;
	cout << "JOB: " << job << endl; 
	cout << "iy: " << iy << " ipt: " << ipt << endl; 
	cout << "y: " << fYbin[iy] << "-" << fYbin[iy+1] << endl; 
	cout << "pT: " << fPTbin[ipt] << "-" << fPTbin[ipt+1] << endl; 
	TStopwatch t; 
	t.Start();
	loop(iy, ipt); 
	t.Stop(); 
	cout << "Real Time: " << t.RealTime() << endl;
	cout << "CPU Time: " << t.CpuTime() << endl; 
	
}

void open_file(TString name){
	input_file = new TFile(name, "READ"); 
}

void initialize_tree(){
	tree = (TTree*)input_file->FindObjectAny("data");
	/*
	if(data){
		tree->SetBranchAddressBranch("UM", &UM, "UM/F");
		tree->Branch("UPt", &UPt, "UPt/F");
		tree->Branch("URapidity", &URapidity, "URapidity/F");
	}
	 */
	
}

double PDF_shape(double *x, double *par){
	//PDF shape for Nevents, which is a sum of gaussians. Fill the parameters for each Gaussian 
	double f=0; 
	double neventd=static_cast<double> (N_event);
	
	for(int i=0; i<N_event; i++){
		int im=2*i;//index for mean
		int is=2*i+1;//index for sigma
		double meanprime=par[im]; // +x[1]
		double sigmaprime=par[is]*x[1];

		double t2=(11.3-meanprime)/(TMath::Sqrt(2)*sigmaprime); 
		double t1=(8.7-meanprime)/(TMath::Sqrt(2)*sigmaprime); 
		
		double Kappa = (1./2)*(TMath::Erf(t2)-TMath::Erf(t1)); 
		
		double norm=1./( TMath::Sqrt(2*TMath::Pi())*sigmaprime*Kappa); // Normalization
		f+=norm*TMath::Gaus(x[0],meanprime,sigmaprime);// Gaussian at point x, with position mean, with width sigma, normalized by 1/2pi
		
	}
	
	return f/neventd; 
}


double get_dm(TH1F *dm_distribution){
	//Function randomly samples a dm distribution, and returns a dm value in units of GeV 
	
	double dm=999; 
	
	TH1F *dm_clone = (TH1F*)dm_distribution->Clone("dm_clone"); 
	dm_clone->SetAxisRange(40,dm_max); 
	
	while (dm*1000>dm_max) {
		dm=dm_clone->GetRandom()/1000.;
	}
	
	delete dm_clone;
	
	return dm; 
	
}

void LS(int ups, int iy, int ipt){
	TString LS_name = Form("LS_y%d_pt%d_ups%d",iy,ipt,ups);
	double dm_scale_max;
	double dm_scale_min;
	if(ipt==fNpt){
		dm_scale_max=1.0+dm_scale_width_last[0];
		dm_scale_min=1.0-dm_scale_width_last[1];
	}
	else {
		dm_scale_min=1.0-dm_scale_width; 
		dm_scale_max=1.0+dm_scale_width; 
	}
	
	TF2 *LS = new TF2(LS_name, PDF_shape, 8.7,11.3,dm_scale_min,dm_scale_max,N_event*2); 

	//Main loop
	//R was 11.6
	int Ntail=static_cast<int>(static_cast<double>(N_event)/(R+1) );

	for(int iSample=0; iSample<N_event; iSample++){
		int im=2*iSample; 
		int is=2*iSample+1; 
		
		double zeta=get_dm(dm);//Get mass uncertainty 
		
		double m=0; 	
		if(iSample<Ntail)m=mass_function->GetRandom(8.7,fill_cutoff[ups-1]); //should be pm-0.002
		else m=PDG_mass[ups-1]; // Get Mass 
			
		LS->SetParameter(im,m);
		LS->SetParameter(is,zeta); 	
	}
	//Write output
	LS->SetNpx(Npx); 
	LS->SetNpy(Npy); 
	output_file->cd();
	LS->Write(); 
	//delete pointers 
	delete LS; 
	
}

double QED_Mass(double *x, double *par){
	double mmuon=0.105658;
	double mmuon2=mmuon*mmuon;
	//Par[0] - arbitrary
	//Par[1] - Upsilon Mass
	//par[2]
	double mups=par[1];
	double mups2=par[1]*par[1];
	double muu=x[0];
	double muu2=x[0]*x[0];
	double beta=TMath::Sqrt(1-4*mmuon2*mups2/TMath::Power(mmuon2+mups2,2));
	double xprime=TMath::Abs(mups-muu);
	double Corr=1/(1-beta*(mups/xprime)*(1-muu2/mups2));
	
	return par[0]*Corr/xprime;
}

void fit_QED_mass(int ups){
	bool print = false; 
	if(print) cout << "fit QED mass: " << ups << endl; 
	if(print) cout << "cutoff: " << cutoff[ups-1] << " fill cutoff: " << fill_cutoff[ups-1] << endl; 
	mass_function = new TF1("mass_function",QED_Mass,8.7,11.5,2);
	mass_function->FixParameter(1,PDG_mass[ups-1]); 
	mass_function->SetRange(8.7,cutoff[ups-1]); 
	TH1F *corrected_peak = new TH1F("corrected_peak",";M_{#mu#mu}^{gen} [GeV];Events",hName[Form("M_mumu%dS",ups)]->GetNbinsX(),8.7,11.3); 
	corrected_peak->SetMarkerColor(kBlue);
	corrected_peak->SetMarkerStyle(8); 
	hName[Form("M_mumu%dS",ups)]->Fit("mass_function","LR");
	
	cout << "chi2: " << mass_function->GetChisquare() << " NDOF: " << mass_function->GetNDF() << endl; 
	
	int peakbin=hName[Form("M_mumu%dS",ups)]->FindBin(PDG_mass[ups-1]);
	
	double P=hName[Form("M_mumu%dS",ups)]->GetBinContent(peakbin);
	double BW=hName[Form("M_mumu%dS",ups)]->GetBinWidth(1);
	double Na=hName[Form("M_mumu%dS",ups)]->Integral(hName[Form("M_mumu%dS",ups)]->FindBin(mlo[ups-1]),hName[Form("M_mumu%dS",ups)]->FindBin(cutoff[ups-1]));
	R=P/Na; 
	cout << "Actual R: " << R << endl; 
	
	double Pba=mass_function->Integral(cutoff[ups-1],fill_cutoff[ups-1])/mass_function->Integral(mlo[ups-1],cutoff[ups-1]); 
	double Nb=Pba*Na;
	R=R-Pba; 
	double Npeak_corr=P-Nb; 

	cout << "Extrapolated R: " << R << endl; 
	/*
	double T=mass_function->Integral(mlo[ups-1],fill_cutoff[ups-1])/BW; // usuall L-pm-0.002... 9.405,9.97,10.3
	double extrap=mass_function->Integral(cutoff[ups-1],fill_cutoff[ups-1])/BW; 
	if(print) cout << "Total number of events: " << T << endl; 
	
	double a=0; 
	double b=0; 
	double NFi=0; 
	double Pi=0;
	double Nhisti=0; 
	double Npeak_corr=0; 
	
	double total_events=hName[Form("M_mumu%dS",ups)]->Integral(hName[Form("M_mumu%dS",ups)]->FindBin(mlo[ups-1]),hName[Form("M_mumu%dS",ups)]->FindBin(PDG_mass[ups-1]+0.05)); 
	
	for(int i=hName[Form("M_mumu%dS",ups)]->FindBin(cutoff[ups-1]); i<hName[Form("M_mumu%dS",ups)]->FindBin(fill_cutoff[ups-1]);i++){
		a=hName[Form("M_mumu%dS",ups)]->GetBinLowEdge(i);
		b=a+BW;
		if(print) cout << "integral from a-b: " << a <<"-"<<b << endl; 
		NFi=mass_function->Integral(a,b)/BW;
//		Pi=NFi/T;
		total_events+=NFi; 
		Pi=NFi/total_events; 
		Nhisti=P*Pi; 
		Npeak_corr+=Nhisti;
	}
	
	double Npeak_corr=P-extrap; 
	
	cout << "Peak: " << P << endl; 
//	cout << "Corrected Peak: " << P+Npeak_corr << endl; 
	cout << "Corrected Peak: " << Npeak_corr << endl; 
	
//	R=(P+Npeak_corr)/T;
	R=Npeak_corr/T;
	 */	
	
	corrected_peak->SetBinContent(corrected_peak->FindBin(PDG_mass[ups-1]),Npeak_corr); 
	
	hName[Form("M_mumu%dS",ups)]->Rebin(4); 
	TF1 *mass_draw = new TF1("mass_draw",QED_Mass,mlo[ups-1],cutoff[ups-1],2); 
	mass_draw->SetParameters(mass_function->GetParameter(0)*4,mass_function->GetParameter(1)); 
	mass_draw->SetLineStyle(1); 
	CreateCanvas(Form("M_mumu_fit_%dS",ups),"",800,600); 
	
	
	TLegend *Leg = new TLegend(0.55,0.55,.85,0.85,"","NDC"); 
	Leg->SetLineColor(10);
	Leg->SetFillColor(10);
	Leg->AddEntry((TF1*)mass_draw->Clone(),"QED Fit","L");
	Leg->AddEntry(corrected_peak,"Extrapolated number of events","P"); 
	Leg->AddEntry(hName[Form("M_mumu%dS",ups)],"PHOTOS Monte Carlo","LEP");
	
	
	hName[Form("M_mumu%dS",ups)]->SetStats(kFALSE);
	hName[Form("M_mumu%dS",ups)]->SetTitle(""); 
	
	
	CName[Form("M_mumu_fit_%dS",ups)]->cd(); 
	gPad->SetLogy();
	hName[Form("M_mumu%dS",ups)]->Draw("histo"); 
	corrected_peak->Draw("E1 same"); 
	mass_draw->DrawClone("same"); 
	mass_draw->SetRange(cutoff[ups-1],PDG_mass[ups-1]-Emin); 
	mass_draw->SetLineStyle(kDashed); 
	Leg->AddEntry((TF1*)mass_draw->Clone(),"Extrapolated QED Fit","L"); 
	mass_draw->DrawClone("same"); 
	Leg->Draw("same"); 
	
	output_file->cd();
	CName[Form("M_mumu_fit_%dS",ups)]->Write(); 
	
	

}

void loop(int iy, int ipt){
	if(data){
		fill_cutoff[0]=PDG_mass[0]-Emin; 
		fill_cutoff[1]=PDG_mass[1]-Emin;
		fill_cutoff[2]=PDG_mass[2]-Emin;
	}
	output_file=new TFile(Form("output_file_y%d_pt%d.root",iy,ipt),"RECREATE"); 
	gen_file = new TFile("gen_mass.root","READ"); 
	if(gen_file->IsOpen()==0) return;
	for(int ups=1; ups<=3; ups++){
		//if(ups!=1) continue; 
		input_file=new TFile("hist_file.root","READ");
		if(input_file->IsOpen()==0)continue;
		get_histogram(ups,iy,ipt);
		fit_QED_mass(ups);
		LS(ups,iy,ipt); 		
	}
	gen_file->Close();
	output_file->Close();
}


void get_histogram(int ups, int iy, int ipt){
	
	dm_m=(TH2F*)input_file->FindObjectAny(Form("zeta_m_y%d_pt%d",iy,ipt)); 
	m=(TH1F*)input_file->FindObjectAny(Form("m_y%d_pt%d",iy,ipt)); 
	for(int ups=1; ups<=3; ups++){
		cout << "Open gen histograms: " << Form("M_mumu%dS",ups) << endl; 
		TH1F *h=(TH1F*)gen_file->FindObjectAny(Form("M_mumu%dS",ups)); 
		hName[Form("M_mumu%dS",ups)]=h; 
	}

	if(ups==1) dm=(TH1F*)dm_m->ProjectionY("dm",m->FindBin(y1m[0]),m->FindBin(y1m[1])); 
	if(ups==2) dm=(TH1F*)dm_m->ProjectionY("dm",m->FindBin(y2m[0]),m->FindBin(y2m[1])); 
	if(ups==3) dm=(TH1F*)dm_m->ProjectionY("dm",m->FindBin(y3m[0]),m->FindBin(y3m[1])); 
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas( Name,  Title,x,y);
	CName[Name]=createC;
}
