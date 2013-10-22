/*
 *  plot_efficiencies.c
 *  
 *
 *  Created by Benjamin Carlson on 7/12/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "plot_efficiencies.h"

void plot_efficiencies(){
	int UPS=1; 
	output_canvas = new TFile(Form("output_file_efficiencies%dS.root",UPS),"RECREATE");
	 cms_pre= Form("CMS Preliminary Simulation, #Upsilon(%dS)",UPS); 
	gROOT->SetBatch();
	load_histograms();
	plot_vertex_efficiency();
	plot_polarization(); 
	draw_histograms("DeltaPhi_DeltaEta_gen",2); 
	draw_histograms("rho_dR_pt_den",2);
	draw_histograms("rho_pt_distM1_den",2);
	draw_histograms("eff_mass_den",1); 
	
	
	projectionY("rho_y_pt",0,0.6); 
	projectionY("rho_y_pt",0.6,1.2); 
	
	divide_histograms("eff_mass",1);
	divide_histograms("rho_cosTheta_pt",2);//2D rho, cosTheta vs pT
	divide_histograms("rho_pt",1); 
	divide_histograms("rho_dR",1); 
	divide_histograms("vertex_y0",1); 


	divide_histograms("sg_y0",1); 
	divide_histograms("sg_y1",1); 

	divide_histograms("rho_dRPtE",1);
	
	divide_histograms("rho_dRPtE_y0",1);
	divide_histograms("rho_dRPtE_y1",1);

	
	divide_histograms("rho_cosTheta",1);//rho as a function of cosTheta
	divide_histograms("rho_cosTheta_dR_LT1p7",1); //rho as a function of cosTheta dR<1.7
	divide_histograms("rho_cosTheta_dR_ge1p7",1); //-- dR>1.7

	divide_histograms("rho_dRPtE_pt",2);
	divide_histograms("rho_DM1",1);
	divide_histograms("rho_dR_pt",2);
	divide_histograms("rho_Dphi_Deta",2);
	divide_histograms("rho_pt_distM1",2);
	divide_histograms("rho_y_pt",2);
	
	//divide_histograms("rho_y_dRPtE",2);
	
	plot_2Deta();
	divide_histograms("rho_cosTheta_DeltaRPtE",2); // 2D rho cosTheta vs dR
	
	for(int iy=0; iy<2; iy++){
		
		compute_acceptance("",iy);
		compute_acceptance("Em",iy);
		compute_acceptance("Ep",iy);
		compute_acceptance("_longitudinal",iy);
		compute_acceptance("_unPol",iy);
		compute_acceptance("_transverse",iy);
		
		compute_acceptance("",iy);
		
		plot_acceptance(iy);
	}
	compute_acceptance("");

	
	output_canvas->cd();
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		bool print = 0; 
		if(print){
			CName[it->first]->Print("/uscms_data/d3/btcarlso/UpsAN/Figures/MC/"+it->first+".pdf"); 
			CName[it->first]->Print("/uscms_data/d3/btcarlso/UpsAN/FigJPG/MC/"+it->first+".jpg"); 
		}
	}
}

void draw_histograms(TString name, int dim){
	CreateCanvas(name,name,600,600);
	CName[name]->cd();
	if(dim==1) hName[name]->DrawCopy(); 
	if(dim==2) hName2D[name]->DrawCopy("COLZ"); 
}

void projectionY(TString name, float y1, float y2){
	TString N=name+Form("_num_%.1f_%.1f",y1,y2); 
	TString D=name+Form("_den_%.1f_%.1f",y1,y2); 


	TH1F *y=(TH1F*)hName2D[name+"_den"]->ProjectionY();
	int BinY1=y->FindBin(y1); 
	int BinY2=y->FindBin(y2);
	
	cout << "Projection from: " << BinY1 << "-" << BinY2 << endl; 
	
	TH1F *h_den=(TH1F*)hName2D[name+"_den"]->ProjectionX(D,BinY1,BinY2);
	TH1F *h_num=(TH1F*)hName2D[name+"_num"]->ProjectionX(N,BinY1,BinY2);
	
	h_num->SetStats(kFALSE); 
	h_num->Divide(h_den); 
	
	
	CreateCanvas(N,N,600,600); 
	CName[N]->cd(); 
	h_num->DrawCopy(); 
	hName[h_num->GetName()]=h_num; 

	
}

void draw_header(){
		L1.SetNDC(kTRUE); 
		L1.DrawLatex(0.15,0.92, cms_pre); 
		
}

void remove_initial_points(TGraph *gr){
	//for some reason polarization plots have non-sense numbers for the first few bins. This procedure removes thos points. 
	//Otherwise the interpolation feature can use the point at say 8,1000
    double x=0;
	
	while (x<10){
		x=gr->GetX()[0]; 
		if(x<10) gr->RemovePoint(0); 
	}
	
}
void plot_polarization(){
	cout << "plot polarization: " << endl; 
	TFile *pol_file = new TFile("/uscms/home/btcarlso/polarization/TGraphResults_1SUps_1sigma.root","READ"); 
	if(pol_file->IsOpen()!=1) 
	{
		cout << "File failed to open " << endl; 
		return; 
	}
	TGraphAsymmErrors *grpol = (TGraphAsymmErrors*)pol_file->FindObjectAny("lth_HX_rap1"); 
	remove_initial_points(grpol); 
	CreateCanvas("polarization_comparison","Polarization", 600,600); 
	CName["polarization_comparison"]->cd(); 
	grpol->Draw("ap"); 
	grpol->SetMarkerStyle(8); 
	grpol->SetMarkerSize(0.5); 
	grpol->GetYaxis()->SetRangeUser(-1,1);
	grpol->GetXaxis()->SetRangeUser(10,100);  
	
	prName["lambda_theta_y0"]->SetAxisRange(-1,1,"Y"); 
	
	prName["lambda_theta_y0"]->Draw("histo"); 
	prName["lambda_thetaEm_y0"]->SetLineColor(kRed);
	prName["lambda_thetaEm_y0"]->Draw("histo same"); 
	prName["lambda_thetaEp_y0"]->Draw("histo same");
	prName["lambda_thetaEp_y0"]->SetLineColor(kGreen);

	TLegend *L = new TLegend(0.2,0.2,0.4,0.4); 
	L->AddEntry(grpol,"Measured Points","LEP"); 
	L->AddEntry(prName["lambda_theta_y0"],"Interpolated","L"); 
	L->AddEntry(prName["lambda_thetaEm_y0"],"Interpolated +68%CL","L"); 
	L->AddEntry(prName["lambda_thetaEp_y0"],"Interpolated -68%CL","L"); 

	grpol->Draw("p same");
	L->Draw("same"); 
	
}

void plot_2Deta(){

	CreateCanvas("eta1_eta2_y0","",700,600); 
	CName["eta1_eta2_y0"]->cd();
	hName2D["eta1_eta2_y0"]->Draw("COLZ");
	draw_header();
	
	TLatex L0(0.65,0.8,"|y|<0.6"); 
	TLatex L1(0.65,0.8,"0.6<|y|<1.2"); 
	L0.SetNDC(kTRUE);
	L1.SetNDC(kTRUE);
	L0.DrawClone();
	
	CreateCanvas("eta1_eta2_y1","",700,600); 
	CName["eta1_eta2_y1"]->cd();
	hName2D["eta1_eta2_y1"]->Draw("COLZ"); 
	draw_header(); 
	L1.DrawClone();
}

void plot_vertex_efficiency(){

	float pT[]={10,12,15,20,30,100};
	float eff[]={0.990,0.990,0.985,0.981,0.998}; 
	float effE[]={0.006,0.007,0.008,0.01,0.025};
	CreateHistogram("vertex_Jpsi","","p_{T} [GeV]","#epsilon_{vp}",5,pT); 
	
	
	for(int i=1; i<=hName["vertex_Jpsi"]->GetNbinsX(); i++){
		hName["vertex_Jpsi"]->SetBinContent(i,eff[i-1]);
		hName["vertex_Jpsi"]->SetBinError(i,effE[i-1]); 
	}
	hName["vertex_Jpsi"]->SetLineColor(kRed); 
	hName["vertex_Jpsi"]->SetStats(kFALSE);
}


void plot_acceptance(int iy){
	TString name=Form("Acceptance_y%d",iy); 
	CreateCanvas(name,name,600,600); 
	CName[name]->cd();
	string method=""; 
	
	TLegend *Lacc= new TLegend(0.35,0.25,0.85,0.5); 
	
	hName[Form("acceptance%s_y%d",method.c_str(),iy)]->Draw("histo p"); 
	TString labeltxt="measured polarization"; 
	Lacc->AddEntry(hName[Form("acceptance%s_y%d",method.c_str(),iy)],labeltxt,"L"); 

	method="Em";
	hName[Form("acceptance%s_y%d",method.c_str(),iy)]->Draw("histo same c"); 
	labeltxt="measured polarization -68%CL"; 
	Lacc->AddEntry(hName[Form("acceptance%s_y%d",method.c_str(),iy)],labeltxt,"L"); 
	
	
	method="Ep";
	hName[Form("acceptance%s_y%d",method.c_str(),iy)]->Draw("histo same c"); 
	
	labeltxt="measured polarization +68%CL"; 
	Lacc->AddEntry(hName[Form("acceptance%s_y%d",method.c_str(),iy)],labeltxt,"L"); 

	
	method="_transverse";
	hName[Form("acceptance%s_y%d",method.c_str(),iy)]->Draw("histo same c"); 
	
	labeltxt="Transverse polarization-HX"; 
	Lacc->AddEntry(hName[Form("acceptance%s_y%d",method.c_str(),iy)],labeltxt,"L"); 

	
	method="_longitudinal";
	hName[Form("acceptance%s_y%d",method.c_str(),iy)]->Draw("histo same c"); 
	
	labeltxt="Longitudinal polarization-HX"; 
	Lacc->AddEntry(hName[Form("acceptance%s_y%d",method.c_str(),iy)],labeltxt,"L"); 
	
	Lacc->Draw("same"); 
	draw_header();


}

void compute_acceptance(string method){
	//Compute 2D acceptance map and plot 
	
	TH2F *N=(TH2F*)hName[Form("acceptance_passed%s",method.c_str())]->Clone(Form("acceptance%s",method.c_str()));
	N->Divide(hName[Form("acceptance_total%s",method.c_str())]); 
	N->SetMinimum(0);
	N->SetMaximum(1); 
	hName2D[N->GetName()]=N; 
	
	CreateCanvas(N->GetName(),N->GetName(),600,600); 
	CName[N->GetName()]->cd();
	hName2D[N->GetName()]->Draw("COLZ"); 
	draw_header();
}

void Rebin(TString name){

	double fYbin[]={0,0.6,1.2};
	double fPTbin[]={10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,43,46,50,55,60,70,100};
	
	const int fNy=sizeof(fYbin)/sizeof(double)-1;
	const int fNpt=sizeof(fPTbin)/sizeof(double)-1;

	for(int i=1; i<=hName["bin_width"]->GetNbinsX(); i++){
		hName["bin_width"]->SetBinContent(i,hName["bin_width"]->GetBinWidth(i)); 
	}
	cout << "binwidth: " << hName["bin_width"]->GetBinContent(13) << endl; 
	
	TH1F *h=(TH1F*)hName[name]->Rebin(fNpt,name+"_rebin",fPTbin); 
	h->Divide(hName["bin_width"]); 
	hName[h->GetName()]=h; 
}

void rho_table(TGraphAsymmErrors *gr){
	ofstream output("rho_table.txt"); 
	int N=7;
	output <<"\\begin{tabular}{"; 
	for(int i=0; i<N; i++) output << "c"; 
	output << "}\\hline" << endl; 
	output << "$p_{T}$ & $\\rho$ & $\\delta \\rho_{pol}$ & $\\delta \\rho_{pol+}$ & $\\delta \\rho_{pol-}$ & $\\delta \\rho_{transverse}$ & $\\delta \\rho_{longidudinal}$" << endl;
	output << "\\\\ \\hline " << endl; 

	float pt_b[]={50,55,60,70,100}; 
	
	double rho[4]; 
	for(int i=0; i<4; i++){
		double avg=0; 
		int n=0; 
		for(int j=0; j<gr->GetN(); j++){
			float x=gr->GetX()[j]; 
			if(x>pt_b[i] && x<=pt_b[i+1]){
				avg+=gr->GetY()[j]; 
				n++; 
			}
		}
		avg=avg/static_cast<float>(n); 
		rho[i]=avg; 
	}
	
	for(int i=0; i<4; i++){
	
		output << pt_b[i] << "-" << pt_b[i+1] << " & "; 
		output << setprecision(2) << rho[i] <<  " & "; 
		string method[]={"","Ep","Em","_transverse","_longitudinal"}; 
		for(int k=0;k<5;k++){
			int Bin=hName[Form("rho_polarizationweighted_pt%s",method[k].c_str())]->FindBin((pt_b[i]+pt_b[i+1])/2); 
			double rho_prime=hName[Form("rho_polarizationweighted_pt%s",method[k].c_str())]->GetBinContent(Bin); 
			cout << "rho: " << rho[i]; 
			cout << " rho_ " <<method[k] << " " << rho_prime << endl; 
			double sys=100*TMath::Abs(rho_prime-rho[i])/rho[i];
			if(k<4)output << setprecision(2) << sys << " & "; 
			else output << setprecision(2) << sys; 
		}
		output << "\\\\ \\hline " << endl; 
	}
	output << "\\hline" << endl << "\\end{tabular}" << endl; 
	
}

void reweight_rho(TH2F *rho, string method){
	cout << "re-weight rho: " << endl; 
	
//	CreateHistogram(Form("rho_polarizationweighted_pt_y%d",iy),"","p_{T} [GeV]","#rho",90,10,100); 
	
	TH1F *Num = new TH1F("Num",";p_{T} [GeV];#rho",rho->GetNbinsX(),10,100); 
	TH1F *Den = new TH1F("Den",";p_{T} [GeV];#rho",rho->GetNbinsX(),10,100); 

	
	for(int ipt=1; ipt<=rho->GetNbinsX(); ipt++){
		double WR=0; 
		double W=0; 
		for(int icos=1; icos<=rho->GetNbinsY(); icos++){
			double weight=hName2D[Form("weight_cosTheta_pt%s",method.c_str())]->GetBinContent(ipt,icos);
	//		cout <<"weight: " << weight << endl; 
			double R=rho->GetBinContent(ipt,icos); 
	//		cout << "rho: " << R << endl; 
			WR+=weight*R;
			W+=weight; 
		}//icos loop	
		Num->SetBinContent(ipt,WR);
		Den->SetBinContent(ipt,W); 
	}//ipt loop 
	
	Num->Divide(Den); 
	TH1F *rho_pol = (TH1F*)Num->Clone(Form("rho_polarizationweighted_pt%s",method.c_str())); 
	hName[Form("rho_polarizationweighted_pt%s",method.c_str())]=rho_pol; 
	delete Num;
	delete Den; 
	
	if(method=="")hName[Form("rho_polarizationweighted_pt%s",method.c_str())]->SetLineColor(kBlack);
	if(method=="_longitudinal")hName[Form("rho_polarizationweighted_pt%s",method.c_str())]->SetLineColor(kOrange);
	if(method=="_transverse")hName[Form("rho_polarizationweighted_pt%s",method.c_str())]->SetLineColor(kBlue);
	if(method=="Ep")hName[Form("rho_polarizationweighted_pt%s",method.c_str())]->SetLineColor(kGreen);
	if(method=="Em")hName[Form("rho_polarizationweighted_pt%s",method.c_str())]->SetLineColor(kRed);
	cout << "polarization weighted rho: " <<  Form("rho_polarizationweighted_pt%s",method.c_str()) << endl; 
	
	output_canvas->cd();
	hName[Form("rho_polarizationweighted_pt%s",method.c_str())]->Write();

	Rebin(Form("rho_polarizationweighted_pt%s",method.c_str())); 
	
	
}

void compute_acceptance(string method, int iy){

	
	TH1F *N=(TH1F*)hName[Form("acceptance_passed%s_y%d",method.c_str(),iy)]->Clone(Form("acceptance%s_y%d",method.c_str(),iy));
	N->Divide(hName[Form("acceptance_total%s_y%d",method.c_str(),iy)]); 
	if(method=="Em") N->SetLineColor(kRed); 
	if(method=="Ep") N->SetLineColor(kGreen);
	if(method=="_transverse")N->SetLineColor(kBlue);
	if(method=="_longitudinal")N->SetLineColor(kOrange); 
	N->SetMinimum(0);
	N->SetMaximum(1); 
	hName[N->GetName()]=N; 

	
}

void divide_histograms(TString name, int dim){
	cout << "name: " << name  << " dim: " << dim << endl; 
	CreateCanvas(name,name,600,600); 
	CName[name]->cd();
	//if(dim==2) cout <<"Dividing" << name << "_num by " << name <<"_den" << endl; 

	
	if(name=="vertex_y0"){ 
		TLegend *L = new TLegend(0.7,0.7,0.85,0.85);
		L->SetLineColor(10);
		L->SetFillColor(10);
		TEfficiency *effy0=new TEfficiency(*hName[name+"_num"],*hName[name+"_den"]);
		name="vertex_y1";
		TEfficiency *effy1=new TEfficiency(*hName[name+"_num"],*hName[name+"_den"]);
		effy1->SetMarkerColor(kBlue);
		
		L->AddEntry(effy0,"MC, |y|<0.6","LEP");
		L->AddEntry(effy1,"MC, 0.6<|y|<1.2","LEP");

		L->AddEntry(hName["vertex_Jpsi"],"data","LEP");
		hName["vertex_Jpsi"]->SetMinimum(0.95);
		hName["vertex_Jpsi"]->SetMaximum(1.05);
		hName["vertex_Jpsi"]->DrawCopy("E1"); 
		effy0->Draw("same");
		effy1->Draw("same");
		L->Draw("same");
		draw_header();
		return; 
	}

	if(dim==1) {
		hName[name+"_num"]->Divide(hName[name+"_den"]); 
		hName[name+"_num"]->DrawCopy(); 
	}
	if(dim==2){
		hName2D[name+"_num"]->Divide(hName2D[name+"_den"]); 
		hName2D[name+"_num"]->SetMinimum(0.65);
		hName2D[name+"_num"]->SetMaximum(1.05);
		hName2D[name+"_num"]->SetTitle("");
		hName2D[name+"_num"]->DrawCopy("COLZ"); 
		draw_header();
		if(name=="rho_cosTheta_pt"){
			reweight_rho(hName2D[name+"_num"],""); 
			reweight_rho(hName2D[name+"_num"],"Ep");
			reweight_rho(hName2D[name+"_num"],"Em");
			reweight_rho(hName2D[name+"_num"],"_longitudinal"); 
			reweight_rho(hName2D[name+"_num"],"_transverse"); 
			
		}
		
	}
	
	output_canvas->cd();

	if(name=="rho_dRPtE" || name =="rho_pt"){
		
		TEfficiency *rho_unc=new TEfficiency(*hName[name+"_num_uW"],*hName[name+"_den"]); 
		TFile *data = new TFile("rho_data.root","READ"); 
		TGraphAsymmErrors *gr_data = (TGraphAsymmErrors*)data->FindObjectAny(name); 
		TGraphAsymmErrors *gr_rho = new TGraphAsymmErrors(hName[name+"_num"]); 

		for(int i=0; i<=gr_rho->GetN(); i++) {
			float yEH=gr_rho->GetY()[i]*rho_unc->GetEfficiencyErrorUp(i)/rho_unc->GetEfficiency(i); 
			float yEL=gr_rho->GetY()[i]*rho_unc->GetEfficiencyErrorLow(i)/rho_unc->GetEfficiency(i); 

			if(!isnan(yEH) && !isnan(yEL)){
				//cout << "effU: " << i << " " <<  yEH << endl; 
				gr_rho->SetPointEYhigh(i,yEH); 
				gr_rho->SetPointEYlow(i,yEL); 
			}
		}
		
		TLegend *L = new TLegend(0.5,0.2,0.8,0.6);
		L->SetLineColor(10);
		L->SetFillColor(10);
		L->AddEntry(gr_rho,"#rho MC","LE");
		L->AddEntry(gr_data,"#rho data","LE");
		gr_rho->SetTitle(""); 
		gr_rho->SetMarkerStyle(7); 
		gr_data->SetMarkerStyle(9);
		gr_data->RemovePoint(0);
		gr_rho->Draw("ap"); 
		gr_rho->GetYaxis()->SetTitle("#rho"); 
		if(name=="rho_dRPtE")gr_rho->GetXaxis()->SetTitle("#DeltaR_{p_{T}}^{elliptic}"); 
		if(name=="rho_pt")gr_rho->GetXaxis()->SetTitle("p_{T} [GeV]"); 
		if(name=="rho_pt")gr_rho->GetXaxis()->SetRangeUser(10,100);
		gr_rho->DrawClone("ap");
		if(name=="rho_dRPtE")gr_data->DrawClone("p same"); 
		if(name=="rho_dRPtE")L->DrawClone("same"); 
		
		if(name=="rho_pt"){
			TString namepol=name+"_polarization_weighted"; 
			CreateCanvas(namepol,namepol,600,600); 
			CName[namepol]->cd();
			gr_rho->DrawClone("ap");
			gr_rho->GetYaxis()->SetRangeUser(0.75,1.05); 
			gr_rho->DrawClone("ap"); 
			hName["rho_polarizationweighted_pt"]->DrawCopy("same histo"); 
			hName["rho_polarizationweighted_pt_longitudinal"]->DrawCopy("same histo"); 
			hName["rho_polarizationweighted_pt_transverse"]->DrawCopy("same histo"); 
			hName["rho_polarizationweighted_ptEp"]->DrawCopy("same histo"); 
			hName["rho_polarizationweighted_ptEm"]->DrawCopy("same histo"); 

			TLegend *Lpol = new TLegend(0.2,0.2,0.6,0.6); 
			Lpol->AddEntry(hName["rho_polarizationweighted_pt"],"#rho-measured pol","L"); 
			Lpol->AddEntry(hName["rho_polarizationweighted_ptEp"],"#rho-measured pol +68%CL","L"); 
			Lpol->AddEntry(hName["rho_polarizationweighted_ptEm"],"#rho-measured pol -68%CL","L"); 

			Lpol->AddEntry(hName["rho_polarizationweighted_pt_longitudinal"],"#rho-longitudinal","L"); 
			Lpol->AddEntry(hName["rho_polarizationweighted_pt_transverse"],"#rho-transverse","L"); 
			Lpol->Draw("same"); 
			draw_header();
			rho_table(gr_rho);

		}
		
		output_canvas->cd();
		gr_rho->Write();
	}
	
	if(dim==1) {
		TH1F *h=(TH1F*)hName[name+"_num"]->Clone(name+"_TH1");
		h->Write();
	}
	
	if(dim==2) {
		TH2F *h=(TH2F*)hName2D[name+"_num"]->Clone(name+"_TH2");
		h->Write();
	}
	
}


void GetProfile(TString name){
	TProfile *pr = (TProfile*)input_file->FindObjectAny(name); 
	pr->SetTitle("");
	prName[name]=pr; 
}

void GetHistogram(TString name){
	TH1F *h=(TH1F*)input_file->FindObjectAny(name);
	h->SetStats(kFALSE);
	h->SetTitle(""); 
	hName[name]=h; 
}

void GetHistogram2D(TString name){
	cout << " name: " << name << endl;
	TH2F *h=(TH2F*)input_file->FindObjectAny(name);
	h->SetStats(kFALSE); 
	h->SetTitle(""); 
	hName2D[name]=h; 
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
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



void load_histograms(){

	GetHistogram("bin_width"); 
	GetHistogram("rho_pt_num");
	GetHistogram("rho_pt_num_uW");
	GetHistogram("rho_pt_den");
	
	GetHistogram("rho_dR_num");
	GetHistogram("rho_dR_den");
	
	GetHistogram("rho_dRPtE_den");
	GetHistogram("rho_dRPtE_num_uW");
	GetHistogram("rho_dRPtE_num");
	
	GetHistogram("rho_dRPtE_y0_den");
	GetHistogram("rho_dRPtE_y0_num_uW");
	GetHistogram("rho_dRPtE_y0_num");
	
	GetHistogram("rho_dRPtE_y1_den");
	GetHistogram("rho_dRPtE_y1_num_uW");
	GetHistogram("rho_dRPtE_y1_num");
	
	GetHistogram("rho_cosTheta_den");
	GetHistogram("rho_cosTheta_num_uW");
	GetHistogram("rho_cosTheta_num");
	
	GetHistogram("rho_cosTheta_dR_LT1p7_den");
	GetHistogram("rho_cosTheta_dR_LT1p7_num_uW");
	GetHistogram("rho_cosTheta_dR_LT1p7_num");
	
	GetHistogram("rho_cosTheta_dR_ge1p7_den");
	GetHistogram("rho_cosTheta_dR_ge1p7_num_uW");
	GetHistogram("rho_cosTheta_dR_ge1p7_num");
	
	GetHistogram("rho_DM1_num");
	GetHistogram("rho_DM1_den");
	
	GetHistogram("eff_mass_den");
	GetHistogram("eff_mass_num");

	GetHistogram2D("rho_cosTheta_pt_num"); 
	GetHistogram2D("rho_cosTheta_pt_num_uW"); 
	GetHistogram2D("rho_cosTheta_pt_den"); 

	GetHistogram2D("rho_cosTheta_DeltaRPtE_num"); 
	GetHistogram2D("rho_cosTheta_DeltaRPtE_num_uW");
	GetHistogram2D("rho_cosTheta_DeltaRPtE_den"); 
	
	
	GetHistogram2D("rho_dR_pt_num");
	GetHistogram2D("rho_dR_pt_den");
	
	GetHistogram2D("eta1_eta2_y0");
	GetHistogram2D("eta1_eta2_y1");
	
//	GetHistogram2D("rho_y_dRPtE_num");
//	GetHistogram2D("rho_y_dRPtE_den");

	
	GetHistogram2D("rho_dRPtE_pt_num");
	GetHistogram2D("rho_dRPtE_pt_den");

	
	GetHistogram2D("rho_y_pt_num");
	GetHistogram2D("rho_y_pt_den");
	
	GetHistogram2D("rho_Dphi_Deta_num");
	GetHistogram2D("rho_Dphi_Deta_den");
	
	GetHistogram2D("rho_pt_distM1_num");
	GetHistogram2D("rho_pt_distM1_den");
	
	GetHistogram("vertex_y0_num");
	GetHistogram("vertex_y0_den");
	GetHistogram("vertex_y1_num");
	GetHistogram("vertex_y1_den");
	
	GetHistogram("sg_y0_num");
	GetHistogram("sg_y0_den");
	
	GetHistogram("sg_y1_num");
	GetHistogram("sg_y1_den");
	
	GetHistogram2D("DeltaPhi_DeltaEta_gen");
	

	GetHistogram2D("weight_cosTheta_pt");
	GetHistogram2D("weight_cosTheta_pt_transverse");
	GetHistogram2D("weight_cosTheta_pt_longitudinal");

	GetHistogram2D("weight_cosTheta_ptEp");
	GetHistogram2D("weight_cosTheta_ptEm");
	
	
	string list[]={"","Em","Ep","_unPol","_transverse","_longitudinal"};
	for (int iy=0; iy<2; iy++) {
		for(int j=0; j<=5;j++){
			GetHistogram(Form("acceptance_passed%s_y%d",list[j].c_str(),iy)); 
			GetHistogram(Form("acceptance_total%s_y%d",list[j].c_str(),iy)); 
		}
	}
	string list2[]={"","_unPol","_transverse","_longitudinal"};

	for(int j=0; j<=3;j++){
		cout << "Method: " <<  list2[j] << endl; 
		GetHistogram(Form("acceptance_passed%s",list2[j].c_str())); 
		GetHistogram(Form("acceptance_total%s",list2[j].c_str())); 

	}
	
	string list3[]={"","Em","Ep"}; 
	for(int j=0; j<3;j++){
	GetProfile(Form("lambda_phi%s_y0",list3[j].c_str())); 
	GetProfile(Form("lambda_phi%s_y1",list3[j].c_str())); 
	
	if(list3[j]=="Em")prName[Form("lambda_phi%s_y0",list3[j].c_str())]->SetLineColor(kRed); 
	if(list3[j]=="Ep")prName[Form("lambda_phi%s_y1",list3[j].c_str())]->SetLineColor(kGreen); 
	
	GetProfile(Form("lambda_theta%s_y0",list3[j].c_str())); 
	GetProfile(Form("lambda_theta%s_y1",list3[j].c_str())); 
	
	if(list3[j]=="Em")prName[Form("lambda_theta%s_y0",list3[j].c_str())]->SetLineColor(kRed); 
	if(list3[j]=="Ep")prName[Form("lambda_theta%s_y1",list3[j].c_str())]->SetLineColor(kGreen); 
	
	GetProfile(Form("lambda_theta_phi%s_y0",list3[j].c_str())); 
	GetProfile(Form("lambda_theta_phi%s_y1",list3[j].c_str())); 
	
	if(list3[j]=="Em")prName[Form("lambda_theta_phi%s_y0",list3[j].c_str())]->SetLineColor(kRed); 
	if(list3[j]=="Ep")prName[Form("lambda_theta_phi%s_y1",list3[j].c_str())]->SetLineColor(kGreen); 
	}
}