/*
 *  compute_xs.c
 *  
 *
 *  Created by Benjamin Carlson on 1/23/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "compute_xs.h"

void initialize(){
	gROOT->SetStyle("Plain"); 
	gStyle->SetPalette(1);
	gErrorIgnoreLevel = kWarning;
	gROOT->SetBatch();
	
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetFillColor(0);
	
	// set the paper & margin sizes
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.12);
	
	// use large Times-Roman fonts
	gStyle->SetTextFont(132);
	gStyle->SetTextSize(0.08);
	gStyle->SetLabelFont(132,"x");
	gStyle->SetLabelFont(132,"y");
	gStyle->SetLabelFont(132,"z");
	gStyle->SetLabelSize(0.05,"x");
	gStyle->SetTitleSize(0.06,"x");
	gStyle->SetLabelSize(0.05,"y");
	gStyle->SetTitleSize(0.06,"y");
	gStyle->SetLabelSize(0.05,"z");
	gStyle->SetTitleSize(0.06,"z");
	
	// use bold lines and markers
//	gStyle->SetMarkerStyle(20);
	gStyle->SetHistLineWidth(1.85);
	//gStyle->SetHistLineWidth(3.85);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	
	// get rid of X error bars and y error bar caps
	gStyle->SetErrorX(0.001);
	
	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	
	gROOT->SetStyle("Plain");
	//gStyle->SetOptStat(1100);
//	gStyle->SetOptStat(1);
//	gStyle->SetOptFit(1111);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	
	
	prelimText.SetTextSize(label_text_size);
	lumiText.SetTextSize(label_text_size);
	prelimTextL.SetTextSize(label_text_size);
	lumiTextL.SetTextSize(label_text_size);
	
	cout << "Compute Cross Section. " << endl; 
	cout << "Number of Pt bins: " << fNpt2 << endl; 
	cout << "Pt bins:"; 
	for(int ipt=0; ipt<fNpt2; ipt++)
		cout << fPTbin2[ipt] << ", "; 
	cout << endl;  
	
	output_files->mkdir("acceptance"); 
	output_files->mkdir("Ratio");
	output_files->mkdir("xs_pulls"); 
	output_files->mkdir("various_otherxs");
	output_files->mkdir("xs_fit"); 
	output_files->mkdir("summary_plots"); 
	output_files->mkdir("mass_shape_pull"); 
	output_files->mkdir("fit_parameters"); 
	output_files->mkdir("yield_fit"); 
	output_files->mkdir("yield_sys"); 
	output_files->mkdir("yield_histograms");
	output_files->mkdir("comparisons"); 
	output_files->mkdir("range_sys"); 
	
	
}

void CreateCanvas(string Name,string Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name.c_str(), Title.c_str(),x,y);
	CName[Name]=createC;
}

void book_canvas(){
	CreateCanvas("fits","fits",600,600);
	CreateCanvas("fits2","fits2",600,600);
	CreateCanvas("xs_canvas","xs_canvas",600,600);
	CreateCanvas("atlas_canvas","atlas_canvas",600,600);
	CreateCanvas("error_summary","error_summary",600,600); 
	CName["error_summary"]->SetFillColor(10); 
	CName["error_summary"]->Divide(2,2); 
}

void compute_xs(){
	
	initialize();
	book_canvas();
	cout << "Compute Yields: " << endl; 
	
	setBW();
	
	compute_yields();
	
	cout << "Compute Acceptance: " << endl; 
	
	acceptance(1); 
	acceptance(2);
	acceptance(3);
	
	cout << "Make Acceptance Plot: " << endl; 
	acceptance_summary_plot();
		
	cout << "Total Systematic: " << endl; 
	stat_error(0,1); 
	stat_error(0,2); 
	stat_error(0,3); 
	
	total_systematic(0,1); 
	total_systematic(0,2); 
	total_systematic(0,3);
	
	plot_total_error(0,1); 
	plot_total_error(0,2); 
	plot_total_error(0,3); 
	
	CName["error_summary"]->cd(4);
	error_legend->Draw();
	
	cout << "Divide Acceptance: " << endl; 
	
	divide_accept(0,1,"");
	divide_accept(0,2,"");
	divide_accept(0,3,"");
	
	cout << "Make xs graph: " << endl; 
	make_xs_stat_graph(0,1);
	make_xs_stat_graph(0,2);
	make_xs_stat_graph(0,3);
	
	CName["xs_canvas"]->SetFillColor(10);
	CName["xs_canvas"]->Divide(1,3);
	
	
	xs_fit(1);
	xs_fit(2);
	xs_fit(3);
	
	reweight_acceptance(1,0); 
	reweight_acceptance(2,0); 
	reweight_acceptance(3,0); 

	CName["xs_canvas"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_canvas.pdf");
	CName["xs_canvas"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_canvas.jpg");

	
	print_tables(0,"");
	cout << "Plot Atlas: " << endl; 
	plot_atlas(1);
	plot_atlas(2);
	plot_atlas(3);
	cout << "Plot CMS: " << endl; 
	plot_cms2010(1);
	plot_cms2010(2);
	plot_cms2010(3);
	
	CName["atlas_canvas"]->SetFillColor(10);
	CName["atlas_canvas"]->Divide(1,3);
	cout << "Compare ATLAS: " << endl; 
	compare_atlas(1);
	compare_atlas(2);
	compare_atlas(3);
	
	output_files->cd("comparisons");
	CName["atlas_canvas"]->Write();
	CName["atlas_canvas"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/atlas_canvas.pdf");
	CName["atlas_canvas"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/atlas_canvas.jpg");

	cout << "Compare CMS: " << endl; 
	compare_cms2010(1);
	compare_cms2010(2);
	compare_cms2010(3);
	
	compare_xs(1);
	compare_xs(2);
	compare_xs(3);

	total_ratio_uncertainty(2);
	total_ratio_uncertainty(3);


	plot_ratios(2,1,"Plot_Pull",false);
	plot_ratios(3,1,"Plot_Pull",false);
	
	plot_ratios(2,1,"Plot_Pull",true);
	plot_ratios(3,1,"Plot_Pull",true);
	
	overlay_xs("");
	overlay_xs("CMS");
	overlay_xs("ATLAS");

	
	CName["error_summary"]->Write();
	CName["error_summary"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/error_summary.pdf"); 
	CName["error_summary"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/error_summary.jpg"); 

		
}


void reweight_acceptance(int ups, int iy){
	cout << "Reweight Acceptance. " << endl; 
	TH1D *acceptance_reweighted = new TH1D(Form("acceptance_reweighted_%dS",ups),Form("acceptance_reweighted_%dS",ups), fNpt2, fPTbin2); 
	TH1D *acceptance1GeV = (TH1D*)output_files->FindObjectAny(Form("acceptance_1GeV%dS",ups));
	TF1 *xs_fit; 
	
									
	for (int i=1; i<=acceptance_reweighted->GetNbinsX(); i++) {
		if(fPTbin2[i]<=20)
			xs_fit= (TF1*)output_files->FindObjectAny(Form("xs_fit_expo_%dS",ups));
		else 
			xs_fit=(TF1*)output_files->FindObjectAny(Form("xs_fit_fcn_%dS",ups));
		

		int nSteps = static_cast<int>((fPTbin2[i]-fPTbin2[i-1])/2);
		int start_bin = acceptance1GeV->FindBin(fPTbin2[i-1]);
		double weighted_acceptanceN=0;//numerator
		double weighted_acceptanceD=0; //denominator
		for (int j=start_bin;j<=start_bin+nSteps; j++) {
			double pt_center = acceptance1GeV->GetBinCenter(j); 
			double weight=xs_fit->Eval(pt_center); 
			double Ao=acceptance1GeV->GetBinContent(j);
			weighted_acceptanceN+=weight*Ao; 
			weighted_acceptanceD+=weight; 
		}
		double Aweighted=weighted_acceptanceN/weighted_acceptanceD; 
		acceptance_reweighted->SetBinContent(i,Aweighted); 
	}
	output_files->cd("acceptance"); 
	acceptance_reweighted->Write(); 
	output_files->cd();
	delete acceptance1GeV; 
	delete acceptance_reweighted; 
	delete xs_fit; 
}

void divide_accept(int iy,int ups,string mode){
	
	TH1D *yield=(TH1D*)output_files->FindObjectAny(yield_histogram_GeV(iy,ups,mode).c_str());
	cout << "Yield before dividing by acceptance: " << yield->GetBinContent(1) << endl;  
	yield->SetName(Form("xs_%dS%s",ups,mode.c_str()));
	TH1D *accept = (TH1D*)output_files->FindObjectAny(Form("acceptance_%dS",ups)); 

	yield->Divide(accept); 
	yield->Scale(1./4900); //Lumi
	output_files->cd();
	yield->Write();
	
	delete yield;
	delete accept; 
}

void make_xs_stat_graph(int iy, int ups){
	TH1D *xs=(TH1D*)output_files->FindObjectAny(Form("xs_%dS",ups));
	TH1D *stat = (TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups)); 
	TH1D *sysP = (TH1D*)output_files->FindObjectAny(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups));
	TH1D *sysM = (TH1D*)output_files->FindObjectAny(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups));

	const int N=xs->GetNbinsX();
	
	double y[N];
	double pt[N];
	double y_eP[N];//y error, P 
	double y_eM[N];//y error, M 

	double xl_e[N];//bin edge low
	double xh_e[N];//bin edge high
	
	double y_stat[N];
	double x_0[N];
	
	for (int i=0; i<N; i++) {
		y[i]=xs->GetBinContent(i+1);
		double SE=stat->GetBinContent(i+1)/100; // stat error
		double SysEP=sysP->GetBinContent(i+1)/100; // sys error +
		double SysEM=sysM->GetBinContent(i+1)/100; //sys error -

		double lE=0.04; //lumi 
		y_eP[i]=TMath::Sqrt(SE*SE+SysEP*SysEP)*y[i]; 
		y_eM[i]=TMath::Sqrt(SE*SE+SysEM*SysEM)*y[i]; 

		y_stat[i]=SE*y[i];
		
		pt[i]=get_bin_center(i,ups);
		//if(pt[i]>=30){
			xl_e[i]=pt[i]-(xs->GetBinCenter(i+1)-xs->GetBinWidth(i+1)/2);
			xh_e[i]=(xs->GetBinWidth(i+1)/2+xs->GetBinCenter(i+1))-pt[i];
		/*
		}
		else{
			xl_e[i]=0; 
			xh_e[i]=0; 
		}
		*/
		x_0[i]=0; 
		
		//cout << "ptbin: " << fPTbin2[i] << "-" << fPTbin2[i+1] << " pt center: " << pt[i] << " +/- " << xl_e[i] << " " << xh_e[i] << " sigma: " << y[i] << " +/- " << y_e[i] << endl;
	}
	
	TGraphAsymmErrors *gr= new TGraphAsymmErrors(N,pt,y,xl_e,xh_e,y_eM,y_eP); 
	gr->SetName(Form("xs_graph_%dS",ups)); 
	output_files->cd();
	gr->Write();
	
	TGraphAsymmErrors *grstat= new TGraphAsymmErrors(N,pt,y,x_0,x_0,y_stat,y_stat); 
	grstat->SetName(Form("xs_stat_graph_%dS",ups)); 
	
	output_files->cd();
	grstat->Write();
	
	delete xs; 
	delete stat;
	delete sysP;
	delete sysM;
	delete gr; 
	delete grstat;
	
}

double ratio_function(double *x, double *par){

	//Par[0] = Anum, Par[1]=Cnum, Par[2]=Alphanum
	//Par[3] = Aden, Par[4]=Cden, Par[5]=Alphaden

	double SN=par[0]/(par[1]+TMath::Power((x[0]/20),par[2])); 
	double SD=par[3]/(par[4]+TMath::Power((x[0]/20),par[5])); 

	double RR=SN/SD;
	
	return RR;
}

void compute_ratio_function(int num,int den){
	TF1 *Fnum=(TF1*)output_files->FindObjectAny(Form("xs_fit_fcn_%dS",num));
	TF1 *Fden=(TF1*)output_files->FindObjectAny(Form("xs_fit_fcn_%dS",den));

	TF1 *RatioF=new TF1(Form("RatioF_%dS_%dS",num,den),ratio_function, 10,100,6); 
	RatioF->SetParameters(Fnum->GetParameter(0),Fnum->GetParameter(2),Fnum->GetParameter(1),Fden->GetParameter(0),Fden->GetParameter(2),Fden->GetParameter(1));
	RatioF->SetName(Form("xs_fit_Ratio_%dS_%dS",num,den));
	output_files->cd("xs_fit");
	RatioF->Write(); 
	
	TF1 *Fexpo_Num=(TF1*)output_files->FindObjectAny(Form("xs_fit_expo_%dS",num)); 
	TF1 *Fexpo_Den=(TF1*)output_files->FindObjectAny(Form("xs_fit_expo_%dS",den)); 
	
	TF1 *expo_R = new TF1(Form("expo_R_%dS_%dS",num,den),"expo",0,30); 
	
	expo_R->SetParameters(Fexpo_Num->GetParameter(0)-Fexpo_Den->GetParameter(0),Fexpo_Num->GetParameter(1)-Fexpo_Den->GetParameter(1));
	expo_R->SetName(Form("expo_R_%dS_%dS",num,den));
	expo_R->Write(); 
	
	delete Fnum; 
	delete Fden; 
	
	delete RatioF; 
	delete Fexpo_Num;
	delete Fexpo_Den; 
	}

void total_ratio_uncertainty(int num){
	cout << "Ratio uncertainty calculations: " << num << "-1" << endl; 

	TFile *output_filesEp=new TFile("output_feb11Ep.root","READ");
	TFile *output_filesEm=new TFile("output_feb11Em.root","READ");
	
	TH1D *r_num=(TH1D*)output_files->FindObjectAny(ratio_histogram(0,num,"").c_str());

	TH1D *rEp_num=(TH1D*)output_filesEp->FindObjectAny(ratio_histogram(0,num,"").c_str());
	TH1D *rEm_num=(TH1D*)output_filesEm->FindObjectAny(ratio_histogram(0,num,"").c_str());

	TH1D *r_sys=(TH1D*)output_files->FindObjectAny(ratio_histogram(0,num,"fit_sys").c_str());
	
	TH1D *An=(TH1D*)output_files->FindObjectAny(Form("acceptance_%dS",num));//acceptance N
	TH1D *Ad=(TH1D*)output_files->FindObjectAny(Form("acceptance_%dS",1));//acceptance D
	
	TH1D *AnEp=(TH1D*)output_files->FindObjectAny(Form("acceptance_error_Ep_%dS",num));//acceptance N
	TH1D *AnEm=(TH1D*)output_files->FindObjectAny(Form("acceptance_error_Em_%dS",num));//acceptance N
	
	TH1D *AdEp=(TH1D*)output_files->FindObjectAny(Form("acceptance_error_Ep_%dS",1));//acceptance N
	TH1D *AdEm=(TH1D*)output_files->FindObjectAny(Form("acceptance_error_Em_%dS",1));//acceptance N
	
	ofstream output_table; 
	output_table.open(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/ratio_%d-1_table.txt",num));
	output_table << "pT " << Form("R(%dS)/R(1S)",num) << " & " << " Efficiency p(m) & "<< " Acceptance p(m) & "  << " fit & " << " systematic & " << " stat & "  << endl;
	
	
	TH1D *RT_uP = new TH1D(ratio_histogram(0,num,"ratio_total_sys_uncertaintyP").c_str(),"Total systematic uncertainty + [%]", fNpt2, fPTbin2); 
	TH1D *RT_uM = new TH1D(ratio_histogram(0,num,"ratio_total_sys_uncertaintyM").c_str(),"Total systematic uncertainty - [%]", fNpt2, fPTbin2); 

	TH1D *Rcorr = new TH1D(ratio_histogram(0,num,"Rcor").c_str(),"Acceptance corrected ratio", fNpt2, fPTbin2); 
	TH1D *R_stat= new TH1D(ratio_histogram(0,num,"Rstat").c_str(),"R statistical uncertainty [%]", fNpt2, fPTbin2); 

	for(int ipt=1; ipt<=r_num->GetNbinsX(); ipt++){
		
		double REp=0; 
		double REm=0; 
		
		double AR=Ad->GetBinContent(ipt)/An->GetBinContent(ipt); 
		double R_fit=r_sys->GetBinContent(ipt); //systematic is in percent
		double R=r_num->GetBinContent(ipt)*AR; 
		output_table << Form("%.0f &",Rcorr->GetBinCenter(ipt)); 
		output_table << Form("%.2f & ",R); 
		
		double stat=100*r_num->GetBinError(ipt)/r_num->GetBinContent(ipt); 
		
		double EEp=100*TMath::Abs(rEp_num->GetBinContent(ipt)-r_num->GetBinContent(ipt));
		double EEm=100*TMath::Abs(rEm_num->GetBinContent(ipt)-r_num->GetBinContent(ipt));

		REp+=TMath::Power(EEp,2);
		REm+=TMath::Power(EEm,2);
		
		output_table << Form("%.1f(%.1f) & ", EEp, EEm); 
		
		double AEp=TMath::Sqrt(TMath::Power(AnEp->GetBinContent(ipt),2)+TMath::Power(AdEp->GetBinContent(ipt),2));
		double AEm=TMath::Sqrt(TMath::Power(AnEm->GetBinContent(ipt),2)+TMath::Power(AdEm->GetBinContent(ipt),2)); 
		
		REp+=TMath::Power(AEp,2);
		REm+=TMath::Power(AEm,2);
		
		output_table << Form("%.1f(%.1f) & ", AEp, AEm); 
		
		output_table << Form("%.1f & ",R_fit); 
		
		REp+=TMath::Power(R_fit,2); 
		REm+=TMath::Power(R_fit,2);
		
		REp=TMath::Sqrt(REp); 
		REm=TMath::Sqrt(REm); 
		
		RT_uP->SetBinContent(ipt,REp); 
		RT_uP->SetBinError(ipt,0);
		RT_uM->SetBinContent(ipt,REm); 
		RT_uM->SetBinError(ipt,0); 
		
		R_stat->SetBinContent(ipt,stat);
		R_stat->SetBinError(ipt,0); 
		Rcorr->SetBinContent(ipt,R); 
		Rcorr->SetBinError(ipt,0); 
		
		output_table << Form("%.1f(%.1f) & ",REp, REm); 
		output_table << Form("%.0f \\\\", stat) << endl; 

	}
	output_table.close();
	output_files->cd("Ratio");
	RT_uP->Write();
	RT_uM->Write();
	R_stat->Write();
	Rcorr->Write();
	
	delete RT_uP;
	delete RT_uM;
	delete R_stat; 
	delete Rcorr; 
	
	
	delete An;
	delete Ad; 
	delete AnEm;
	delete AnEp;
	delete AdEm;
	delete AdEp;
	
	delete r_num;
	delete rEp_num;
	delete rEm_num;
	delete r_sys; 
	
	output_filesEp->Close();
	output_filesEm->Close();
	delete output_filesEp;
	delete output_filesEm;
	
}

void plot_ratios(int num, int den, string PP, bool show_fits){
	cout << Form("Plot Ratios: %d-%d",num,den) << endl; 
	compute_ratio_function(num,den); 
	
	TGraphAsymmErrors *gr_num=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",num));
	TGraphAsymmErrors *gr_den=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",den));
	
	TH1D *r_num=(TH1D*)output_files->FindObjectAny(ratio_histogram(0,num,"Rcor").c_str());
	TH1D *r_stat=(TH1D*)output_files->FindObjectAny(ratio_histogram(0,num,"Rstat").c_str());
	TH1D *r_sysP=(TH1D*)output_files->FindObjectAny(ratio_histogram(0,num,"ratio_total_sys_uncertaintyP").c_str());
	TH1D *r_sysM=(TH1D*)output_files->FindObjectAny(ratio_histogram(0,num,"ratio_total_sys_uncertaintyP").c_str());

	TGraphAsymmErrors *gr2010=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("Ratio2010_%d_%d",num,den));
	gr2010->SetLineColor(kRed); 

	TF1 *xs_fit=(TF1*)output_files->FindObjectAny(Form("xs_fit_fcn_%dS",num)); 
	TF1 *fit_ratio = (TF1*)output_files->FindObjectAny(Form(Form("xs_fit_Ratio_%dS_%dS",num,den))); 
	TF1 *expo_ratios = (TF1*)output_files->FindObjectAny(Form("expo_R_%dS_%dS",num,den)); 
	
	const int N=r_num->GetNbinsX();

	//cout << "Number of Points: " << N << endl; 
	
	double *pt_den=gr_den->GetX();
	double *pt_num=gr_num->GetX();

	double y[N];
	double pt[N];
	double y_le[N];//y error, low 
	double y_he[N];//y error, high 

	double ystat_le[N];//y error, low 
	double ystat_he[N];//y error, high 
	
	double xl_e[N];//bin edge low
	double xh_e[N];//bin edge high
	
	double zero[N];
	cout << "Start filling graph: " << endl; 
	for (int i=0; i<N; i++) {

	
		double R=r_num->GetBinContent(i+1); // Acceptance corrected ratio  
		double RP=R*r_sysP->GetBinContent(i+1)/100; 
		double RM=R*r_sysM->GetBinContent(i+1)/100; 
		
		double stat=R*r_stat->GetBinContent(i+1)/100;

		double TEP=TMath::Sqrt(RP*RP+stat*stat);
		double TEM=TMath::Sqrt(RM*RM+stat*stat);

		pt[i]=pt_den[i];
		y[i]=R;
		ystat_le[i]=stat;  
		ystat_he[i]=stat;
		
		y_le[i]=TEM; 
		y_he[i]=TEP; 
		
		xl_e[i]=gr_den->GetEXlow()[i];
		xh_e[i]=gr_den->GetEXhigh()[i];
		
		zero[i]=0; 
		//cout << "pt: " << pt[i] << " y: " << y[i] << " stat: " << ystat_he[i] << " tot: " << y_he[i] << endl; 
	}

	TGraphAsymmErrors *ratio=new TGraphAsymmErrors(N,pt,y,xl_e,xh_e,y_le,y_he); 
	ratio->SetName(Form("Ratio_2011_%dS_%dS",num,den));
	
	TGraphAsymmErrors *ratio_stat=new TGraphAsymmErrors(N,pt,y,zero,zero,ystat_le,ystat_he); 
	ratio_stat->SetName(Form("Ratio_stat_2011_%dS_%dS",num,den));
	
	TCanvas *C = new TCanvas(Form("Ratio%dS-%dS",num,den), "",800,800);
	ratio->SetTitle("");
	ratio->SetMinimum(0);
	double max=1.3;
	ratio->SetMaximum(max);
	set_CMS2011(ratio,"ratio",1);
	set_CMS2010(gr2010,"ratio"); 
	set_CMS2011(ratio_stat,"stat",1); 
	
	ratio->GetXaxis()->SetRangeUser(0,100);
	gr2010->GetXaxis()->SetRangeUser(0,100);
	ratio->GetXaxis()->SetTitle("p_{T}(#mu#mu) [GeV]"); 
	ratio->GetYaxis()->SetTitle(Form("#sigma#timesBr(%dS)/#sigma#timesBr(%dS)",num,den));

	
	TLegend *L = new TLegend(0.3,0.2,0.65,0.35,"","NDC");//0.15,0.73,0.45,0.93
	L->SetTextSize(label_text_size*0.9);
	L->SetFillColor(10); 
	L->SetLineColor(10);
	L->AddEntry(ratio,"CMS 2011 |y(#mu#mu)|<0.6","LEP"); 
	L->AddEntry(gr2010, "CMS 2010 |y(#mu#mu)|<2.4, 36 pb^{-1}","LEP"); 
	
	TMultiGraph * grM = new TMultiGraph(); 
	grM->Add(gr2010); 
	grM->Add(ratio); 
	grM->Add(ratio_stat,"||"); 
	grM->Draw("AP"); 
	grM->GetHistogram()->GetXaxis()->SetRangeUser(0,100);
	grM->GetHistogram()->GetXaxis()->SetTitle("p_{T}(#mu#mu) [GeV]");
	grM->GetHistogram()->GetYaxis()->SetTitle(Form("#sigma#timesBr(%dS)/#sigma#timesBr(%dS)",num,den));
	//grM->GetHistogram()->GetYaxis()->SetLabelSize(0.02);

	grM->SetMinimum(0);
	grM->SetMaximum(1); 
	if(num==2 && den==1) grM->SetMaximum(1);
	if(num==3 && den==1) grM->SetMaximum(1);
	if(num==3 && den==2) grM->SetMaximum(1);
	grM->Draw("AP");
	//Draw ratios
	fit_ratio->SetRange(20,100); 
	fit_ratio->SetLineWidth(1);
	expo_ratios->SetLineColor(kRed);
	expo_ratios->SetLineWidth(1);
	if(show_fits)L->AddEntry(expo_ratios,"#scale[0.9]{Exponential fit: 10<p_{T}<20 GeV}","L"); 
	if(show_fits)L->AddEntry(fit_ratio,"#scale[0.9]{Power-law fit: 20<p_{T}<100 GeV}","L"); 
	
	L->DrawClone("same"); 

	
	if(show_fits)fit_ratio->DrawClone("same"); 
	
	fit_ratio->SetRange(10,20); 
	fit_ratio->SetLineStyle(kDashed); 
	if(show_fits)fit_ratio->DrawClone("same"); 
	
	expo_ratios->SetRange(10,20); 
	if(show_fits)expo_ratios->DrawClone("same"); 
	
	expo_ratios->SetRange(0,10); 
	expo_ratios->SetLineStyle(kDashed); 
	if(show_fits)expo_ratios->DrawClone("same");
	
	expo_ratios->SetRange(20,30); 
	expo_ratios->SetLineStyle(kDashed); 
	if(show_fits)expo_ratios->DrawClone("same"); 
	
	prelimTextL.SetNDC(kTRUE); 
	lumiTextL.SetNDC(kTRUE);
	prelimTextL.DrawLatex(0.15,0.85,cms_pre);
	lumiTextL.DrawLatex(0.15,0.85-label_text_size,lumi_string);
	TLatex lumi_exc_txt;
	lumi_exc_txt.SetNDC(kTRUE);
	lumi_exc_txt.SetTextSize(0.025);
	//lumi_exc_txt.DrawLatex(0.15,0.85-label_text_size*2.1,lumi_exc);
	output_files->cd("Ratio"); 
	ratio->Write(); 
	ratio_stat->Write();
	
	C->Write();
	if(show_fits==false)C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/Ratio%dS-%dS_nofits.pdf",num,den));
	else C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/Ratio%dS-%dS.pdf",num,den));

	if(show_fits==false)C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/Ratio%dS-%dS_nofits.jpg",num,den));
	else C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/Ratio%dS-%dS.jpg",num,den));
	
	delete C;

	output_files->cd();
	cout <<"Free pointers: " << endl; 
	delete L; 
	
	/*delete gr_den;
	delete gr_num;
	delete r_num; 
	delete r_sysP;
	delete r_sysM;
	delete r_stat;
	 */
	delete expo_ratios;
	delete fit_ratio;
	delete xs_fit; 
	delete gr2010;
	
	delete ratio;
	delete ratio_stat; 

	
}

void scale_graph(TGraphAsymmErrors *gr, double SF){
	//Scale TGraphAsymmErrors 
	
	for (int i=0;i<gr->GetN();i++){
		gr->GetY()[i] *= SF;
		gr->GetEYhigh()[i]*=SF; 
		gr->GetEYlow()[i]*=SF;
	}
}

void clear_points(TGraphAsymmErrors *gr, double Xmin, double Xmax){
	//clear points between Xmin-Xmax
	//cout << "Clearing points from: " << Xmin << "-" << Xmax << endl;
	
	int k=0; 
	
	for(int i=0; i<gr->GetN(); i++){
		double x = gr->GetX()[i];
		if(x>Xmin && x<=Xmax) {
			//cout << "Removing Point: " << x << " i " << i << endl; 
			gr->RemovePoint(k); 
		}
		else {
			k++;
		}

	}
	
	
}

void set_titles(TMultiGraph *GR){
	//cout << "Set Titles: " << endl; 
	GR->GetHistogram()->GetXaxis()->SetTitle(x_label.c_str()); 
	GR->GetHistogram()->GetYaxis()->SetTitle(xs_y.c_str()); 	
	
}

void set_CMS2010(TGraphAsymmErrors *gr,string MODE){

	gr->SetMarkerStyle(5); 
	gr->SetMarkerSize(1); //0.2
	gr->SetLineWidth(0.025); 
	gr->SetMarkerColor(kBlue); 
	gr->SetFillStyle(3005);
	gr->SetFillColor(9);
//	if(MODE=="ratio") gr->SetMarkerSize(0.75); 
	gr->SetLineColor(kBlue); 
}

void set_CMS2011(TGraphAsymmErrors *gr,string MODE, int ups){
	
	if(ups==1) gr->SetMarkerStyle(20); 
	if(ups==2) gr->SetMarkerStyle(21);
	if(ups==3) gr->SetMarkerStyle(22); 
	
	gr->SetMarkerSize(0.5); 
	gr->SetLineWidth(0.1); 
	if(ups==3 && MODE!="stat") gr->SetMarkerSize(0.35); 
	if(MODE=="ratio") gr->SetMarkerSize(0.75); 
	
	if(MODE=="stat"){
		gr->SetLineColor(kRed); 
		//clear_points(gr,10,40);
	}
}

void set_ATLAS(TGraphAsymmErrors *gr){
	gr->SetMarkerStyle(24); 
	gr->SetMarkerSize(0.5); //makersize(0.5), width(0.025), marker 24
	gr->SetFillColor(9); 
	gr->SetFillStyle(3005); 
	gr->SetLineWidth(0.025); 
	gr->SetLineColor(kBlue); 
}

void overlay_xs(string experiment){
	int iy=0;
	TCanvas *C = new TCanvas("","",800,800); 

	TMultiGraph *gRM= new TMultiGraph(); 
	
	TLegend *L = new TLegend(0.17,.2,0.4,0.37,"","NDC");
	L->SetFillColor(10); 
	L->SetLineColor(10);
	L->SetTextSize(label_text_size);
	
	TGraphAsymmErrors *xs_graph1 = (TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",1));
	TGraphAsymmErrors *xs_graph_stat1=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_stat_graph_%dS",1)); 
		
	TGraphAsymmErrors *xs_graph2 = (TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",2));
	TGraphAsymmErrors *xs_graph_stat2=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_stat_graph_%dS",2)); 
	
	TGraphAsymmErrors *xs_graph3 = (TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",3));
	TGraphAsymmErrors *xs_graph_stat3=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_stat_graph_%dS",3)); 
	
	
	TGraphAsymmErrors *cms2010Y1S=(TGraphAsymmErrors*)output_files->FindObjectAny(Form(Form("cms2010_%dS",1)));
	TGraphAsymmErrors *cms2010Y2S=(TGraphAsymmErrors*)output_files->FindObjectAny(Form(Form("cms2010_%dS",2)));
	TGraphAsymmErrors *cms2010Y3S=(TGraphAsymmErrors*)output_files->FindObjectAny(Form(Form("cms2010_%dS",3)));

	TGraphAsymmErrors *gr_atlas1 =(TGraphAsymmErrors*)output_files->FindObjectAny(Form("Graph_Atlas_%dS",1));
	TGraphAsymmErrors *gr_atlas2 =(TGraphAsymmErrors*)output_files->FindObjectAny(Form("Graph_Atlas_%dS",2));
	TGraphAsymmErrors *gr_atlas3 =(TGraphAsymmErrors*)output_files->FindObjectAny(Form("Graph_Atlas_%dS",3));

	
	double SY2=0.1;
	double SY3=0.01;
	
	scale_graph(xs_graph2,SY2); 
	scale_graph(xs_graph3,SY3);
	
	scale_graph(xs_graph_stat2,SY2); 
	scale_graph(xs_graph_stat3,SY3);
	
	scale_graph(cms2010Y2S, SY2); 
	scale_graph(cms2010Y3S, SY3); 
	
	scale_graph(gr_atlas2, SY2); 
	scale_graph(gr_atlas3, SY3); 
	
	set_CMS2011(xs_graph1,"",1);
	set_CMS2011(xs_graph2,"",2);
	set_CMS2011(xs_graph3,"",3);
	
	set_CMS2011(xs_graph_stat1,"stat",1);
	set_CMS2011(xs_graph_stat2,"stat",2);
	set_CMS2011(xs_graph_stat3,"stat",3);

	set_CMS2010(cms2010Y1S,""); 
	set_CMS2010(cms2010Y2S,""); 
	set_CMS2010(cms2010Y3S,""); 
	set_ATLAS(gr_atlas1);
	set_ATLAS(gr_atlas2);
	set_ATLAS(gr_atlas3);
 
	gRM->Add(xs_graph1);
	gRM->Add(xs_graph_stat1,"||");
		
	gRM->Add(xs_graph2);
	gRM->Add(xs_graph_stat2,"||");
	
	gRM->Add(xs_graph3);
	gRM->Add(xs_graph_stat3,"||");
	//set_titles(gRM); 
	
	gPad->SetLogy(); 
	//gRMstat->Draw("||");
	TLatex UpsPeak;
	UpsPeak.SetNDC(kTRUE);
	UpsPeak.SetTextSize(label_text_size*0.7);

	L->AddEntry(xs_graph1,"CMS 2011 |y(#mu#mu)|<0.6","LP");

	if(experiment==""){
		gRM->Draw("ap");
		gRM->GetXaxis()->SetRangeUser(10,100);
		set_titles(gRM);
		gRM->Draw("ap"); 
		
		//	L->Draw("same");
		prelimText.SetNDC(kTRUE); 
		lumiText.SetNDC(kTRUE);
		prelimText.Draw();
		lumiText.Draw();
		
		UpsPeak.DrawLatex(0.25,0.83,Form("#Upsilon(%dS)",1)); 
		UpsPeak.DrawLatex(0.25,0.65,Form("#Upsilon(%dS)#times%.1f",2,SY2)); 
		UpsPeak.DrawLatex(0.25,0.4,Form("#Upsilon(%dS)#times%.2f",3,SY3)); 
		UpsPeak.SetTextSize(label_text_size);
		UpsPeak.DrawLatex(label_x+0.1, label_y-3*label_text_size,rapidity_text(iy).c_str());
		UpsPeak.DrawLatex(0.2,0.2,lumi_exc);
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_points_overlay.pdf"); 
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_points_overlay.jpg"); 
	}	
	
	if(experiment=="CMS"){
		double max=1.4*(cms2010Y1S->Eval(4.0)+cms2010Y1S->GetErrorYhigh(cms2010Y1S->GetXaxis()->FindBin(4))); 
		gRM->SetMaximum(max);
		gRM->Add(cms2010Y1S);
		gRM->Add(cms2010Y2S);
		gRM->Add(cms2010Y3S); 
		L->AddEntry(cms2010Y1S,"CMS 2010 |y(#mu#mu)|<2.4#times0.25, 36 pb^{-1}", "LAP");
		gRM->Draw("ap"); 
		set_titles(gRM);
		gRM->GetXaxis()->SetRangeUser(0,100);
		gRM->Draw("ap"); 
		L->Draw("same");
		prelimText.SetNDC(kTRUE); 
		lumiText.SetNDC(kTRUE);
		prelimText.Draw();
		lumiText.Draw();
		
		UpsPeak.DrawLatex(0.25,0.85,Form("#Upsilon(%dS)",1)); 
		UpsPeak.DrawLatex(0.15,0.8,Form("#Upsilon(%dS)#times%.1f",2,SY2)); 
		UpsPeak.DrawLatex(0.17,0.5,Form("#Upsilon(%dS)#times%.2f",3,SY3)); 
		UpsPeak.DrawLatex(label_x, label_y-label_text_size*3,lumi_exc);
		
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_points_cms2010_overlay.pdf"); 
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_points_cms2010_overlay.jpg"); 

	}
	
	if(experiment=="ATLAS"){
		double max=1.4*(gr_atlas1->Eval(4.0)+gr_atlas1->GetErrorYhigh(gr_atlas1->GetXaxis()->FindBin(4))); 
		gRM->SetMaximum(max);
		gRM->Add(gr_atlas1, "p2");
		gRM->Add(gr_atlas2,"p2");
		gRM->Add(gr_atlas3,"p2"); 
		L->AddEntry(gr_atlas1,"ATLAS 2011 |y(#mu#mu)|<1.2, 1.8 fb^{-1}", "f");
		gRM->Draw("ap");
		set_titles(gRM);
		gRM->GetXaxis()->SetRangeUser(0,100);
		gRM->Draw("ap"); 
		L->Draw("same");
		prelimText.SetNDC(kTRUE); 
		lumiText.SetNDC(kTRUE);
		prelimText.Draw();
		lumiText.Draw();
		
		UpsPeak.DrawLatex(0.35,0.7,Form("#Upsilon(%dS)",1)); 
		UpsPeak.DrawLatex(0.35,0.57,Form("#Upsilon(%dS)#times%.1f",2,SY2)); 
		UpsPeak.DrawLatex(0.17,0.5,Form("#Upsilon(%dS)#times%.2f",3,SY3)); 
		
		UpsPeak.DrawLatex(0.45,label_y-3.5*label_text_size,"ATLAS points normalized to ");
		UpsPeak.DrawLatex(0.45,label_y-4.5*label_text_size,"CMS 40<p_{T}(#mu#mu)<50 GeV");
		UpsPeak.DrawLatex(0.45,label_y-5.5*label_text_size,lumi_exc);

		
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_points_atlas_overlay.pdf"); 
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_points_atlas_overlay.jpg"); 

	}
	
	delete xs_graph1;
	delete xs_graph_stat1; 
	delete xs_graph2;
	delete xs_graph_stat2; 
	delete xs_graph3;
	delete xs_graph_stat3; 
	
	delete cms2010Y1S;
	delete cms2010Y2S;
	delete cms2010Y3S;

	delete gr_atlas1;
	delete gr_atlas2;
	delete gr_atlas3;
	
	delete gRM;
	
	delete L; 
	delete C; 
}

void xs_fit_hist(TF1 *xs_fit,int ups){
	cout <<"XS fit option I: " << endl; 
	int iy=0; 
	TH1D *xsP=(TH1D*)((TH1D*)output_files->FindObjectAny(Form("xs_%dS",ups)))->Clone("xsP");
	TH1D *xsM=(TH1D*)((TH1D*)output_files->FindObjectAny(Form("xs_%dS",ups)))->Clone("xsM");

	TH1D *stat = (TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups)); 
	TH1D *sysP = (TH1D*)output_files->FindObjectAny(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups));
	TH1D *sysM = (TH1D*)output_files->FindObjectAny(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups));
	
	for (int ipt=1; ipt<=xsP->GetNbinsX(); ipt++) {
		double errP=0; 
		double errM=0; 
		double xs=xsP->GetBinContent(ipt);
		double S=stat->GetBinContent(ipt);
		double P=sysP->GetBinContent(ipt); 
		double M=sysM->GetBinContent(ipt);
		
		errP=xs*TMath::Sqrt(S*S+P*P)/100;
		errM=xs*TMath::Sqrt(S*S+M*M)/100; 
		
		xsP->SetBinError(ipt,errP); 
		xsM->SetBinError(ipt,errM); 
	}
	double A=0; 
	double alpha=0; 
	double C=0; 
	
	double AE=0; 
	double alphaE=0; 
	double CE=0; 
	
	xsP->Fit(xs_fit,"QRI"); 
	
	A+=xs_fit->GetParameter(0); 
	alpha+=xs_fit->GetParameter(1);
	C+=xs_fit->GetParameter(2); 
	
	AE+=xs_fit->GetParError(0);
	alphaE+=xs_fit->GetParError(1);
	CE+=xs_fit->GetParError(2); 
	
	xsM->Fit(xs_fit,"QRI"); 
	
	A+=xs_fit->GetParameter(0); 
	alpha+=xs_fit->GetParameter(1);
	C+=xs_fit->GetParameter(2); 
	
	AE+=xs_fit->GetParError(0);
	alphaE+=xs_fit->GetParError(1);
	CE+=xs_fit->GetParError(2); 
	
	A=A/2;
	alpha=alpha/2;
	C=C/2;
	
	AE=AE/2;
	alphaE=alphaE/2;
	CE=CE/2;
	//average both fit parameters and uncertainties
	
	xs_fit->SetParameters(A,alpha,C); 
	xs_fit->SetParError(0,AE);
	xs_fit->SetParError(1,alphaE);
	xs_fit->SetParError(2,CE); 
		
	delete xsP; 
	delete xsM; 
	
	delete stat; 
	delete sysP; 
	delete sysM; 

	
}

void xs_fit(int ups){
	int iy=0;
	TGraphAsymmErrors *xs_graph = (TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",ups));
	TGraphAsymmErrors *xs_graph_stat=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_stat_graph_%dS",ups)); 
	TF1 *xs_fit = new TF1("xs_fit", power_law,20,100,3); 
	xs_fit->SetParNames("A","alpha","C");
	
	xs_fit->SetParameters(1,5,1);
	
	TF1 *xs_low = new TF1("xs_low",power_law,10,20,3);	
	TF1 *expo = new TF1("expo", "expo",10,20);

	TLegend *L = new TLegend(0.4,.55,0.85,0.7,"","NDC");
	L->SetFillColor(10); 
	L->SetTextSize(label_text_size*0.9);
	expo->SetParameters(15,.5);
	
	TCanvas *C = new TCanvas(Form("xs_%d_fit",ups),"xs fit",800,800);
	TPaveText *txt = new TPaveText(.1,.10,.95,.95);

	set_CMS2011(xs_graph,"",ups);
	set_CMS2011(xs_graph_stat,"stat",ups);
	
	xs_graph->Fit(xs_fit,"QR0");
	xs_fit_hist(xs_fit,ups); 
	//xs_fit->SetRange(20,100);
	cout << "Fit expo: " << endl;
	xs_graph->Fit(expo,"QR0"); 
	
	cout << "After expo fit. " << endl; 
	xs_low->SetParameters(xs_fit->GetParameter(0),xs_fit->GetParameter(1),xs_fit->GetParameter(2));
	xs_low->SetLineStyle(kDashed);
	expo->SetLineColor(kRed); 
	
	txt->AddText(Form("A=%.2f +/- %.3f",xs_fit->GetParameter(0),xs_fit->GetParError(0))); 
	txt->AddText(Form("#alpha=%.2f +/- %.3f",xs_fit->GetParameter(1),xs_fit->GetParError(1))); 
	txt->AddText(Form("C=%.2f +/- %.3f",xs_fit->GetParameter(2),xs_fit->GetParError(2))); 
	txt->AddText(Form("#chi^{2}=%.2f",xs_fit->GetChisquare()));
	txt->AddText(Form("NDOF=%d",xs_fit->GetNDF()));
	//txt->AddText(Form("#chi^{2}/NDOF=%.2f",xs_fit->GetChisquare()/static_cast<double>(xs_fit->GetNDF())));
	
	txt->AddText("Exponential Fit: ");
	txt->AddText(Form("A=%.1f +/- %.2f",expo->GetParameter(0),expo->GetParError(0))); 
	txt->AddText(Form("B=%.1f +/- %.2f",expo->GetParameter(1),expo->GetParError(1))); 
	txt->AddText(Form("#chi^{2}=%.2f",expo->GetChisquare()));
	txt->AddText(Form("NDOF=%d",expo->GetNDF()));
	//txt->AddText(Form("#chi^{2}/NDOF=%.2f",expo->GetChisquare()/static_cast<double>(expo->GetNDF())));
	
	xs_graph->GetXaxis()->SetTitle(x_label.c_str()); 
	xs_graph->GetYaxis()->SetTitle(xs_y.c_str()); 
	
	xs_graph->SetTitle(""); 
	double max=xs_graph->Eval(10.);
	double min=xs_fit->Eval(100);
	xs_graph->SetMaximum(max*1.1);
	xs_graph->SetMinimum(min); 
	gPad->SetLogy();
	xs_graph->GetXaxis()->SetRangeUser(10,100);
	xs_graph_stat->GetXaxis()->SetRangeUser(10,100);
	xs_graph_stat->SetLineColor(kRed);
	xs_graph->Draw("ap"); 
	xs_graph_stat->Draw("|| same");
	xs_low->SetLineWidth(1);
	xs_fit->SetLineWidth(1);
	xs_fit->Draw("same");
	xs_low->Draw("same"); 
	expo->SetLineWidth(1);
	
	L->SetLineColor(10);
	L->AddEntry(expo,"Exponential fit: 10<p_{T}<20 GeV","L"); 
	L->AddEntry(xs_fit,"Power-law fit: 20<p_{T}<100 GeV","L");  

	L->DrawClone("same");
	expo->DrawClone("same"); 
	expo->SetRange(20,35);
	expo->SetLineStyle(kDashed); 
	expo->Draw("same");
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	TLatex UpsPeak(0.45,label_y-3.*label_text_size,Form("#Upsilon(%dS), %s",ups,rapidity_text(iy).c_str()));
	UpsPeak.SetTextSize(label_text_size);
	UpsPeak.SetNDC(kTRUE);
	UpsPeak.Draw();
	UpsPeak.DrawLatex(0.2,0.2, lumi_exc); 
	//UpsPeak.DrawLatex(0.5,label_y-4.*label_text_size,rapidity_text(iy).c_str()); 
	
	TCanvas *CC = new TCanvas("","",800,800);
	txt->DrawClone();		  
				  
	output_files->cd("xs_fit");
	xs_fit->SetName(Form("xs_fit_fcn_%dS",ups));
	xs_fit->Write();
	expo->SetName(Form("xs_fit_expo_%dS",ups)); 
	expo->Write();
	
	CName["xs_canvas"]->cd(ups); 
	gPad->SetLogy();
	xs_graph->SetMarkerStyle(7); 
	xs_graph->DrawClone("ap"); 
	xs_fit->DrawClone("same");
	xs_low->DrawClone("same");
	expo->SetRange(10,20); 
	expo->SetLineStyle(1); 
	expo->DrawClone("same"); 
	expo->SetRange(20,35); 
	expo->SetLineStyle(kDashed);
	expo->DrawClone("same"); 
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	UpsPeak.DrawLatex(0.45,label_y-3.1*label_text_size,Form("#Upsilon(%dS), %s",ups, rapidity_text(iy).c_str()));
	UpsPeak.DrawLatex(0.2,0.2, lumi_exc); 
	//UpsPeak.DrawLatex(0.47,label_y-4.1*label_text_size,rapidity_text(iy).c_str()); 

	
	C->Write();
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_fit_%dS.pdf",ups));
	CC->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_fit_params_%dS.pdf",ups));
	
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_fit_%dS.jpg",ups));
	CC->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_fit_params_%dS.jpg",ups));
	
	C->cd();
	xs_graph->Draw("AP");
	xs_graph_stat->Draw("|| same");
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	UpsPeak.DrawLatex(0.5,label_y-3.1*label_text_size,Form("#Upsilon(%dS), %s",ups, rapidity_text(iy).c_str()));
	UpsPeak.DrawLatex(0.2,0.2, lumi_exc); 

	
	//UpsPeak.Draw();
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_points_only_%dS.pdf",ups));
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_points_only_%dS.jpg",ups));

	delete txt;
	delete C;
	delete CC;
	delete expo; 
	delete xs_fit;
	delete xs_low;
	delete xs_graph; 
	delete xs_graph_stat;
}

void plot_cms2010(int ups){
		//Table 38, BPH-11-001
	TGraphAsymmErrors *gr; 
	TGraphAsymmErrors *grRatio;
	
	double pTbins[]={0,2,5,8,10,13,16,18,22,38};

	int nR=9; 
	double ptR[9];
	double ptRE[9];
	
	//cout << "pT bins: " << endl; 
	
	for (int i=0; i<nR; i++) {
		ptR[i]=(pTbins[i+1]-pTbins[i])/2.+pTbins[i];
		ptRE[i]=pTbins[i+1]-ptR[i];
		//cout << pTbins[i+1] << "-" << pTbins[i] << "ptR: " << ptR[i] << " +/- " << ptRE[i] << endl;  
	}
	
	if(ups==1){

		//Table 58, page 97 from AN2011_107_v6
		//xs from http://arxiv.org/pdf/1303.5900v1.pdf
		// and exact numbers can be downloaded from the twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsBPH11001
		// though it's in a terrible format
		
		//xs measured for |y|<2.4, then scaled 
		double ptB[]={0,0.5,1,1.5,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,25,30,50};
		double pT[]={0.327,0.772,1.257,1.751,2.489,3.481,4.478,5.485,6.487,7.493,8.484,9.478,10.479,11.485,12.487,13.475,14.489,15.477,16.909,18.978,20.937,23.284,27.024,35.968};

		const int N=sizeof(ptB)/sizeof(double)-1;
		double ptP[N];
		double ptM[N];
	
		double xs[]={0.0859,	0.263,	0.374,	0.505,	1.16,	1.212,	1.084,	0.879,	0.68,	0.556,	0.419,	0.331,	0.238,	0.1793,	0.1451,	0.099,	0.075,	0.0595,	0.0732,0.05,0.0302,	0.0237,	0.0205,	0.0123};
		double xsEp[]={0.0094894678,0.0219317122,0.0297321375,0.0436806593,0.0921954446,0.0905207159,0.0680073525,0.0633245608,0.0481040539,0.0338378486,0.0256320112,0.0187882942,0.0134164079,0.0101355809,0.0090354856,0.0058830264,0.0054203321,0.0042720019,0.0048877398,0.0034669872,0.0024207437,0.0020615528,0.0019209373,0.001421267};

		double xsEm[]={0.0084599054,0.021023796,0.0278567766,0.0323109888,0.1118033989,0.0696419414,0.0596154342,0.071063352,0.0362353419,0.0273130006,0.0201246118,0.017,0.0116619038,0.0101355809,0.0078230429,0.0059682493,0.004580393,0.0038600518,0.0043011626,0.0032202484,0.0024207437,0.0020615528,0.0019209373,0.0013453624};
	
		for (int i=0; i<N; i++) {
			ptP[i]=(ptB[i+1]-pT[i]); 
			ptM[i]=(pT[i]-ptB[i]); 
			double BW=ptB[i+1]-ptB[i]; 
			xs[i]=xs[i]/BW; 
			//subtract lumi uncertainty
			double sysTot=xsEm[i];
			xsEm[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			sysTot=xsEp[i];
			xsEp[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			
			xsEm[i]=xsEm[i]/BW; 
			xsEp[i]=xsEp[i]/BW;
		}
		
		gr = new TGraphAsymmErrors(N,pT,xs,ptM,ptP,xsEp,xsEm);
		scale_graph(gr,1000*(0.25)); // scale from nb, then correction for rapidity scaling. 
		
		//|y|<2.4 acceptance corrected 
		
		double R[]={0.11,0.10,0.12,0.15,0.20,0.26,0.31,0.28,0.32};
		double Rp[]={0.02,0.02,0.02,0.02,0.03,0.030,0.04,0.04,0.07};
		double Rm[]={0.02,0.01,0.02,0.02,0.02,0.03,0.04,0.04,0.09};
	
		grRatio=  new TGraphAsymmErrors(nR,ptR,R,ptRE,ptRE,Rp,Rm); 
		
		grRatio->SetName(Form("Ratio2010_%d_%d",3,1));
		grRatio->GetYaxis()->SetLabelSize(0.025); 
		grRatio->GetYaxis()->SetTitle(Form("#sigma#timesBr(%dS)/#sigma#timesBr(%dS)",3,1)); 
		grRatio->GetXaxis()->SetRangeUser(0,100);
		
	}//ups 1S
	
	
	if(ups==2){
		
		double ptB[]={0,1,2.5,4,5.5,7,8.5,10,11.5,13,14.5,16,18,19.5,22,26,42};
		double pT[]={0.656,1.785,3.214,4.708,6.217,7.711,9.206,10.691,12.213,13.698,15.222,16.875,18.757,20.647,23.688,31.3};
		
		
		const int N=sizeof(ptB)/sizeof(double)-1;
		double ptP[N];
		double ptM[N];
		
		
		double xs[]={0.0829,0.331,0.409,0.362,0.286,0.212,0.146,0.1123,0.0765,0.0519,0.0376,0.0373,0.0159,0.0204,0.0158,0.0126};
		double xsEp[]={0.0092962358,0.0386005181,0.0382753184,0.0304630924,0.0223606798,0.0187882942,0.0125299641,0.0097529483,0.0063513778,0.0043600459,0.0034205263,0.0033600595,0.0017691806,0.0021260292,0.0017029386,0.0015620499};
		//{13,12,12,12,13,13,13,15,14};
		double xsEm[]={0.0087920419,0.0245967478,0.0299666481,0.0259422435,0.0214709106,0.0144222051,0.01,0.0074027022,0.0056222771,0.0041231056,0.0039962482,0.0032015621,0.0017691806,0.0021260292,0.0017029386,0.0014866069};
		//{10,10,10,11,12,12,12,14,13};
	
		for (int i=0; i<N; i++) {
			ptP[i]=(ptB[i+1]-pT[i]); 
			ptM[i]=(pT[i]-ptB[i]); 
			double BW=ptB[i+1]-ptB[i]; 
			xs[i]=xs[i]/BW; 
			
			double sysTot=xsEm[i];
			xsEm[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			sysTot=xsEp[i];
			xsEp[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			
			xsEm[i]=xsEm[i]/BW; 
			xsEp[i]=xsEp[i]/BW;
		}
		
		gr = new TGraphAsymmErrors(N,pT,xs,ptM,ptP,xsEp,xsEm);
		scale_graph(gr,1000*(0.25)); // scale from nb, then correction for rapidity scaling. 

		double R[]={0.22,0.23,0.27,0.27,0.33,0.37,0.50,0.45,0.47};
		double Rp[]={0.03,0.03,0.04,0.04,0.04,0.04,0.06,0.05,0.09};
		double Rm[]={0.03,0.03,0.03,0.03,0.03,0.04,0.06,0.06,0.13};
		
		grRatio = new TGraphAsymmErrors(nR,ptR,R,ptRE,ptRE,Rp,Rm); 
		
		grRatio->SetName(Form("Ratio2010_%d_%d",2,1));
		grRatio->GetYaxis()->SetLabelSize(0.025); 
		grRatio->GetYaxis()->SetTitle(Form("#sigma#timesBr(%dS)/#sigma#timesBr(%dS)",2,1)); 
		grRatio->GetXaxis()->SetRangeUser(0,100);
	}//ups 2S
	
	if(ups==3){
		double ptB[]={0,2.5,5,7.5,10,13,16,18,22,38};
		double pT[]={1.539,3.62,6.147,8.618,11.306,14.298,16.943,19.702,26.514};
		
		const int N=sizeof(ptB)/sizeof(double)-1;
		double ptP[N];
		double ptM[N];
		
		double xs[]={0.203,0.287,0.227,0.157,0.1127,0.0617,0.0227,0.0229,0.0185};
		double xsEp[]={0.0210950231,0.0281780056,0.0232594067,0.0196977156,0.0099849887,0.0051400389,0.0026248809,0.0024083189,0.0021260292};
		double xsEm[]={0.0194164878,0.0317804972,0.0205912603,0.0136014705,0.007981228,0.0052201533,0.0026248809,0.0026400758,0.0030413813};
	
		for (int i=0; i<N; i++) {
			ptP[i]=(ptB[i+1]-pT[i]); 
			ptM[i]=(pT[i]-ptB[i]); 
			double BW=ptB[i+1]-ptB[i]; 
			xs[i]=xs[i]/BW; 
			
			double sysTot=xsEm[i];
			xsEm[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			sysTot=xsEp[i];
			xsEp[i]=TMath::Sqrt(sysTot*sysTot-TMath::Power(0.04*sysTot,2)); 
			
			xsEm[i]=xsEm[i]/BW; 
			xsEp[i]=xsEp[i]/BW;
		}
		
		gr = new TGraphAsymmErrors(N,pT,xs,ptM,ptP,xsEp,xsEm);
		scale_graph(gr,1000*(0.25)); // scale from nb, then correction for rapidity scaling. 
		
		double R[]={0.49,0.44,0.46,0.57,0.60,0.69,0.61,0.63,0.69};
		double Rp[]={0.08,0.09,0.07,0.09,0.08,0.09,0.1,0.1,0.2};
		double Rm[]={0.08,0.06,0.06,0.08,0.06,0.07,0.09,0.08,0.17};
		
		grRatio = new TGraphAsymmErrors(nR,ptR,R,ptRE,ptRE,Rp,Rm); 
		grRatio->SetName(Form("Ratio2010_%d_%d",3,2));
		grRatio->GetYaxis()->SetTitle(Form("#sigma#timesBr(%dS)/#sigma#timesBr(%dS)",3,2)); 
		grRatio->GetYaxis()->SetLabelSize(0.025); 
		grRatio->GetXaxis()->SetRangeUser(0,100);
	}//ups 3S
	gr->GetXaxis()->SetTitle(x_label.c_str()); 
	gr->GetYaxis()->SetTitle(xs_y.c_str()); 
	gr->SetName(Form("cms2010_%dS",ups)); 
	
	grRatio->GetXaxis()->SetTitle("p_{T}(#mu#mu)");
	grRatio->SetMarkerStyle(8);
	grRatio->SetMarkerSize(0.5);
	gr->SetName(Form("cms2010_%dS",ups)); 
	
	output_files->cd("various_otherxs");
	gr->Write();
	output_files->cd("Ratio"); 
	grRatio->Write(); 
	delete grRatio; 
	delete gr; 
	
}

void plot_atlas(int ups){
	
	TCanvas *C = new TCanvas(Form("Atlas_%dS",ups),"",800,800);
	double SF=1000; // 1000 for fb->pb, 0.6 for rapidity region

	TGraphAsymmErrors *xs=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",ups));
	
	// Note that Atlas plots d^2sigma/dPtdy... so we have to account for this
	
	if(ups==1){
		// Plot: p8322_d6x1y2
		double p8322_d6x1y1_xval[] = { 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 
			4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 
			9.75, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 
			19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 
			29.5, 30.5, 31.5, 33.0, 35.0, 37.0, 39.0, 42.5, 47.5, 55.0, 
			65.0 };
		/*
		double p8322_d6x1y1_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
			5.0 };
		
		double p8322_d6x1y1_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
			5.0 };
		 
		 */
		double p8322_d6x1y1_yval[] = { 35392.84, 98878.16, 148372.15, 182173.1, 206775.92, 231574.2, 242509.16, 240950.58, 229797.15, 
			208978.78, 184938.25, 179761.51, 157285.63, 136729.51, 122675.42, 108822.11, 94001.57, 82014.84, 71983.66, 
			62278.51, 49124.29, 35832.87, 27727.3, 20263.43, 16302.79, 12331.68, 9516.72, 7338.76, 5701.8, 
			4476.17, 3595.25, 2801.82, 2288.65, 1867.8, 1472.76, 1240.35, 988.08, 780.36, 662.63, 
			599.29, 488.3, 374.2, 309.67, 236.99, 164.84, 126.49, 74.05, 41.24, 18.07, 
			7.28 };
		double p8322_d6x1y1_yerrminus[] = { 10040.187744644021, 29165.59695075347, 44279.714260927656, 53244.39693433479, 58505.326055196034, 66398.84614706342, 67324.27643658638, 67963.08338160136, 65933.42573497679, 
			61006.88931563139, 54157.40316778953, 52706.095821221665, 46209.00046205717, 38514.3298287196, 35650.19575819044, 32012.60791104811, 26752.9963248792, 23204.191100292635, 20099.002776060806, 
			17090.635482351146, 13269.239589588395, 9396.921362808142, 7052.913343363862, 5006.5960821500275, 3853.3447551445483, 2811.3084500104214, 2090.0416370493676, 1552.4720581382455, 1173.3702550772282, 
			886.0189593343925, 684.7933779615571, 519.8145781911085, 413.08033613329985, 327.06105760239933, 251.20286105058597, 204.92675276790973, 158.94243014374732, 122.73967655163509, 103.61023260276951, 
			89.65957282967615, 69.64123275761278, 54.57913612361413, 42.612845481145705, 31.766877403987948, 21.940068368170596, 16.76881927864929, 8.930761445699913, 5.295431993709295, 2.098952119511067, 
			1.0630616162763098 };
		double p8322_d6x1y1_yerrplus[] = { 8699.38653125035, 32576.66833363719, 76805.90114189795, 236700.99427909675, 453550.8479124082, 624837.9455057994, 712761.5664403944, 718525.2938768687, 666355.8395618561, 
			575792.9900605836, 468661.24430696777, 419350.19024839095, 327654.08364405914, 253070.89968113083, 197740.74979823153, 147822.2229745988, 104222.93167987072, 75840.32422017855, 56664.81719483528, 
			43287.58068620376, 29173.177920019272, 18089.10823410596, 11864.652619925288, 7629.401854418208, 5552.797011849073, 3822.807598022689, 2722.015197992105, 1958.100781088655, 1400.19382069055, 
			1038.064565718337, 789.4447244107722, 586.6845467540456, 455.53080170280475, 355.2579786577636, 269.27636565432175, 218.1831874824456, 168.5934043787004, 128.73072088666325, 103.87868212487103, 
			93.38044441958927, 75.79567072069486, 56.37493592014096, 43.5152789259129, 32.44108660325668, 22.11096108268476, 16.809961332495682, 9.070622911355096, 5.281098370604358, 2.1842389979120873, 
			1.0938921336219583 };
		double p8322_d6x1y1_ystatminus[] = { 936.86, 2155.91, 2587.05, 2985.8, 3191.88, 2809.87, 3533.43, 3492.8, 3058.9, 
			2937.09, 2939.08, 2576.18, 2420.64, 1824.43, 1731.36, 1984.32, 1559.85, 1322.67, 1157.24, 
			746.27, 663.51, 426.62, 319.35, 332.39, 182.0, 142.19, 126.55, 107.58, 88.57, 
			76.01, 63.94, 56.97, 51.14, 46.03, 40.61, 37.7, 32.44, 30.2, 26.21, 
			24.61, 22.76, 19.91, 12.43, 11.26, 9.21, 7.96, 4.0, 2.96, 1.38, 
			0.89 };
		double p8322_d6x1y1_ystatplus[] = { 936.86, 2155.91, 2587.05, 2985.8, 3191.88, 2809.87, 3533.43, 3492.8, 3058.9, 
			2937.09, 2939.08, 2576.18, 2420.64, 1824.43, 1731.36, 1984.32, 1559.85, 1322.67, 1157.24, 
			746.27, 663.51, 426.62, 319.35, 332.39, 182.0, 142.19, 126.55, 107.58, 88.57, 
			76.01, 63.94, 56.97, 51.14, 46.03, 40.61, 37.7, 32.44, 30.2, 26.21, 
			24.61, 22.76, 19.91, 12.43, 11.26, 9.21, 7.96, 4.0, 2.96, 1.38, 
			0.89 };
		int p8322_d6x1y1_numpoints = 50;
		
		double p8322_d6x1y1_xerrminus[]={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
			5.0 };
		double p8322_d6x1y1_xerrplus[]={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
			0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 5.0, 
			5.0 };
		
		double cmspt;
		double cmsxs=0;
		double xstmp;
		
		
		xs->GetPoint(15,cmspt,xstmp);
		cout << "x range: " << cmspt; 
		cmsxs+=xstmp;
		xs->GetPoint(16,cmspt,xstmp);
		cmsxs+=xstmp;
		xs->GetPoint(17,cmspt,xstmp);
		cmsxs+=xstmp;
		cmsxs=cmsxs/3;
		cout << "-" << cmspt << endl;
	
		int ref=0; 
		
		for (int i=0; i<p8322_d6x1y1_numpoints; i++) {
			p8322_d6x1y1_yval[i]=p8322_d6x1y1_yval[i]/SF;
			p8322_d6x1y1_yerrminus[i]=p8322_d6x1y1_yerrminus[i]/SF;
			p8322_d6x1y1_yerrplus[i]=p8322_d6x1y1_yerrplus[i]/SF;
			if(i<p8322_d6x1y1_numpoints-3){ 
				//p8322_d6x1y1_xerrminus[i]=0;
				//p8322_d6x1y1_xerrplus[i]=0;
			}
			if(p8322_d6x1y1_xval[i]>39.0 && ref==0)ref=i;
			
		}
		
		double atlasxs=((p8322_d6x1y1_yval[ref]+ p8322_d6x1y1_yval[ref+1])/2);
		
		TGraphAsymmErrors *p8322_d6x1y1 =new TGraphAsymmErrors(p8322_d6x1y1_numpoints, p8322_d6x1y1_xval, p8322_d6x1y1_yval, p8322_d6x1y1_xerrminus, p8322_d6x1y1_xerrplus, p8322_d6x1y1_yerrminus, p8322_d6x1y1_yerrplus);
		//scale_graph(p8322_d6x1y1,(cmsxs/(atlasxs*0.9))); 
		scale_graph(p8322_d6x1y1,(1.2));
		
		//p8322_d6x1y1->Draw("AP");
		output_files->cd("various_otherxs");
		p8322_d6x1y1->SetName(Form("Graph_Atlas_%dS",ups));
		p8322_d6x1y1->Write();
		delete p8322_d6x1y1;
		
		
	}//1S
	if(ups==2) {
		// Plot: p8322_d7x1y2
		double p8322_d7x1y1_xval[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 
			9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 
			29.0, 31.0, 34.0, 38.0, 45.0, 60.0 };
		/*
		double p8322_d7x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		double p8322_d7x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		 */
		double p8322_d7x1y1_yval[] = { 14099.19, 37169.61, 43895.41, 53632.18, 51016.7, 41440.45, 35410.02, 29636.38, 22632.57, 
			19053.93, 12611.54, 7939.33, 5343.06, 3309.25, 2072.13, 1337.29, 897.91, 608.07, 418.21, 
			296.1, 209.3, 128.9, 74.37, 31.83, 7.51 };
		double p8322_d7x1y1_yerrminus[] = { 3915.9007974794254, 10660.241072564917, 12853.494024929563, 14293.457798353062, 13238.90561274987, 11292.63679040905, 9961.513892185263, 8653.921212323346, 6430.445462843145, 
			5246.714875557847, 3297.1416707960852, 1967.388458693402, 1217.0259453273788, 705.8831526534685, 410.33545277004765, 246.36233255106188, 157.04912002300426, 101.24744243683392, 67.16017048816954, 
			44.87941398904402, 31.301439264033853, 17.62770830255595, 10.05802167426577, 3.9284348028190563, 1.1434159348198711 };
		double p8322_d7x1y1_yerrplus[] = { 3584.398011884283, 14539.50576821647, 48442.20717914285, 101196.61661934454, 106239.46591437947, 80524.4152455173, 58953.63710351805, 39473.663513499734, 22059.308130478617, 
			13904.259904140888, 7055.801357018209, 3287.4933701225923, 1749.2841637652814, 920.9084363279554, 497.0995806676968, 286.2217689135472, 173.65078692594514, 108.14296371008149, 70.38218595639098, 
			47.09316192399911, 32.442011343318406, 17.764475787368454, 10.074691062260918, 3.9676315353117153, 1.162497311824849 };
		double p8322_d7x1y1_ystatminus[] = { 462.5, 901.0, 1013.67, 1219.81, 1122.06, 989.2, 818.33, 1227.45, 664.19, 
			465.43, 283.93, 153.6, 82.14, 59.1, 43.89, 33.56, 27.06, 21.34, 18.39, 
			14.47, 12.56, 6.82, 5.13, 2.21, 0.67 };
		double p8322_d7x1y1_ystatplus[] = { 462.5, 901.0, 1013.67, 1219.81, 1122.06, 989.2, 818.33, 1227.45, 664.19, 
			465.43, 283.93, 153.6, 82.14, 59.1, 43.89, 33.56, 27.06, 21.34, 18.39, 
			14.47, 12.56, 6.82, 5.13, 2.21, 0.67 };
		int p8322_d7x1y1_numpoints = 25;
		
		double p8322_d7x1y1_xerrminus[]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		double p8322_d7x1y1_xerrplus[]= { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		
		int ref=0; 
		
		for (int i=0; i<p8322_d7x1y1_numpoints; i++) {
			//cout <<  	"xs before scaling: " << p8322_d7x1y1_yval[i] << endl; 
			p8322_d7x1y1_yval[i]=p8322_d7x1y1_yval[i]/SF;
			p8322_d7x1y1_yerrminus[i]=p8322_d7x1y1_yerrminus[i]/SF;
			p8322_d7x1y1_yerrplus[i]=p8322_d7x1y1_yerrplus[i]/SF;
		    if(p8322_d7x1y1_xval[i]>39.0 && ref==0)ref=i;	
			//cout <<  	"xs: " << p8322_d7x1y1_yval[i] << endl; 
		}
		
		double atlasxs=p8322_d7x1y1_yval[ref];
		
		
		
		double cmspt;
		double cmsxs=0;
		double xstmp;
		
		xs->GetPoint(15,cmspt,xstmp);
		cout << "x range: " << cmspt; 
		cmsxs+=xstmp;
		xs->GetPoint(16,cmspt,xstmp);
		cmsxs+=xstmp;
		xs->GetPoint(17,cmspt,xstmp);
		cmsxs+=xstmp;
		cmsxs=cmsxs/3;

	
		TGraphAsymmErrors *p8322_d7x1y1 =new TGraphAsymmErrors(p8322_d7x1y1_numpoints, p8322_d7x1y1_xval, p8322_d7x1y1_yval, p8322_d7x1y1_xerrminus, p8322_d7x1y1_xerrplus, p8322_d7x1y1_yerrminus, p8322_d7x1y1_yerrplus);
		scale_graph(p8322_d7x1y1,(cmsxs/(atlasxs*0.9))); 
		
		output_files->cd("various_otherxs");
		p8322_d7x1y1->SetName(Form("Graph_Atlas_%dS",ups));
		p8322_d7x1y1->Write();
		delete p8322_d7x1y1;

		}//2S
	if(ups==3){
		// Plot: p8322_d8x1y2
		double p8322_d8x1y1_xval[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 
			9.5, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 
			29.0, 31.0, 34.0, 38.0, 45.0, 60.0 };
		
		/*
		double p8322_d8x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		double p8322_d8x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		*/
		 
		 double p8322_d8x1y1_yval[] = { 6022.32, 14021.13, 20711.77, 20317.61, 19951.64, 18330.26, 17320.64, 14677.17, 12151.44, 
			9357.4, 6510.87, 4765.1, 3157.02, 1992.37, 1334.81, 946.0, 610.91, 437.57, 331.78, 
			208.81, 161.25, 85.68, 53.61, 19.27, 5.49 };
		double p8322_d8x1y1_yerrminus[] = { 1572.0742719413736, 3795.0331310016254, 5680.281874528763, 5146.731171403846, 4966.7803433008785, 4885.1679069403535, 5283.628766387737, 4121.8282818550315, 3441.0830164644385, 
			2577.152806102114, 1724.769824179447, 1227.0599322771484, 718.3371423642244, 425.4963911715351, 266.156395564713, 178.58454944367386, 107.63820371968309, 73.39877110687888, 55.09720319580659, 
			34.224321760993305, 24.980154122823183, 12.657551105960426, 8.39753535270915, 2.8718460961548757, 1.1015443704181869 };
		double p8322_d8x1y1_yerrplus[] = { 1487.2148013316703, 4576.14641170931, 14303.538135066443, 28691.751653766278, 33960.546032234226, 31241.88712183373, 25458.392024750898, 18944.708991259275, 11234.75530593791, 
			6797.286830199238, 3500.615076811502, 1952.7710345557668, 1037.0850676776713, 552.0350199942029, 322.970171068475, 204.0989519816307, 119.8237192712695, 79.175352225298, 57.12247368593205, 
			33.601421100899884, 25.62148707628033, 12.793592146070626, 7.988372800514507, 2.885186302476843, 1.1015443704181869 };
		double p8322_d8x1y1_ystatminus[] = { 357.98, 537.42, 639.56, 820.97, 750.21, 762.74, 698.45, 583.39, 587.97, 
			401.56, 233.58, 135.86, 73.7, 53.92, 40.56, 31.02, 24.73, 20.24, 17.8, 
			13.49, 12.29, 6.38, 4.99, 2.05, 0.71 };
		double p8322_d8x1y1_ystatplus[] = { 357.98, 537.42, 639.56, 820.97, 750.21, 762.74, 698.45, 583.39, 587.97, 
			401.56, 233.58, 135.86, 73.7, 53.92, 40.56, 31.02, 24.73, 20.24, 17.8, 
			13.49, 12.29, 6.38, 4.99, 2.05, 0.71 };
		const int p8322_d8x1y1_numpoints = 25;
		
		double p8322_d8x1y1_xerrminus[]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		double p8322_d8x1y1_xerrplus[]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
			0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
			1.0, 1.0, 2.0, 2.0, 5.0, 10.0 };
		
		
		double cmspt;
		double cmsxs=0;
		double xstmp;
		
		xs->GetPoint(15,cmspt,xstmp);
		cout << "x range: " << cmspt; 
		cmsxs+=xstmp;
		xs->GetPoint(16,cmspt,xstmp);
		cmsxs+=xstmp;
		xs->GetPoint(17,cmspt,xstmp);
		cmsxs+=xstmp;
		cmsxs=cmsxs/3;
		
		int ref=0; 
		for (int i=0; i<p8322_d8x1y1_numpoints; i++) {
			p8322_d8x1y1_yval[i]=p8322_d8x1y1_yval[i]/SF;
			p8322_d8x1y1_yerrminus[i]=p8322_d8x1y1_yerrminus[i]/SF;
			p8322_d8x1y1_yerrplus[i]=p8322_d8x1y1_yerrplus[i]/SF;
			if(i<p8322_d8x1y1_numpoints-3){
				//p8322_d8x1y1_xerrminus[i]=0; 
				//p8322_d8x1y1_xerrplus[i]=0; 
			}
			if(p8322_d8x1y1_xval[i]>39.0 && ref==0)ref=i;	

		}
		double atlasxs=p8322_d8x1y1_yval[ref];

		TGraphAsymmErrors *p8322_d8x1y1 =new TGraphAsymmErrors(p8322_d8x1y1_numpoints, p8322_d8x1y1_xval, p8322_d8x1y1_yval, p8322_d8x1y1_xerrminus, p8322_d8x1y1_xerrplus, p8322_d8x1y1_yerrminus, p8322_d8x1y1_yerrplus);

		scale_graph(p8322_d8x1y1,(cmsxs/(atlasxs*0.9))); 
		
		
		cout << "cmsxs: " << cmsxs << endl;
		cout << "atlasxs: " << atlasxs << endl; 
		cout << "R: " << cmsxs/(atlasxs*0.9) << endl; 
		output_files->cd("various_otherxs");
		p8322_d8x1y1->SetName(Form("Graph_Atlas_%dS",ups));
		p8322_d8x1y1->Write();
		
		delete p8322_d8x1y1;
		
		}//3S

 }

void compare_atlas(int ups){
	int iy=0;
	TGraphAsymmErrors *xs=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",ups)); 
	TGraphAsymmErrors *xs_stat=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_stat_graph_%dS",ups)); 
	
	TGraphAsymmErrors *gr_atlas =(TGraphAsymmErrors*)output_files->FindObjectAny(Form("Graph_Atlas_%dS",ups));
	
	gr_atlas->SetMarkerStyle(7);
	gr_atlas->SetTitle("");
	xs->GetXaxis()->SetTitle(x_label.c_str()); 
	xs->GetYaxis()->SetTitle(xs_y.c_str()); 
	
	
	set_CMS2011(xs,"",ups);
	set_CMS2011(xs_stat,"stat",ups); 
	xs->SetTitle(""); 
	
	set_ATLAS(gr_atlas); 
	
	TLegend *L = new TLegend(0.15,0.2,0.45,0.35,"","NDC");
	L->SetTextSize(label_text_size*0.9);
	L->SetFillColor(10); 
	L->SetLineColor(10);
	L->AddEntry(xs,"CMS 2011 |y(#mu#mu)|<0.6","LP"); 	
	//L->AddEntry(xs_stat,"Statistical Uncertainty","L"); 
	//L->AddEntry(xs,"Total Uncertainty","L"); 
	L->AddEntry(gr_atlas, "ATLAS 2011 |y(#mu#mu)|<1.2, 1.8 fb^{-1}","fp"); 
	
	gr_atlas->SetName("Atlas"); 
	TCanvas *comp = new TCanvas(Form("cms_atlas_comp_%dS",ups),"",800,800);
	
	double atlas_max=gr_atlas->Eval(4)+gr_atlas->GetErrorYhigh(gr_atlas->GetXaxis()->FindBin(4));
	//cout << "ATLAS MAX: " << atlas_max << endl; 
	double max = TMath::Max(xs->Eval(10.),atlas_max);
	
	
	
	xs->SetMaximum(max*1.3); 
	gPad->SetLogy();
	gr_atlas->GetXaxis()->SetRangeUser(0,100);
	xs->GetXaxis()->SetRangeUser(0,100); 
	//xs->SetMarkerStyle(7);
	
	xs->Draw("AP same");
	xs_stat->Draw("|| same"); 
	gr_atlas->Draw("p2 same");
	
	/*
	TMultiGraph *gRM=new TMultiGraph();
	
	gRM->Add(xs,"ap"); 
	gRM->Add(xs_stat,"||");
	gRM->Add(gr_atlas,"ap");
	
	gRM->Draw("ap");
	*/
	//gr_atlas->Draw("P same");
	L->Draw("same"); 
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	TLatex UpsPeak(0.5,label_y-3*label_text_size,Form("#Upsilon(%dS)",ups));
	UpsPeak.SetNDC(kTRUE);
	UpsPeak.SetTextSize(label_text_size); 
	UpsPeak.Draw();
	//UpsPeak.DrawLatex(0.47,label_y-4.1*label_spacing,rapidity_text(iy).c_str()); 
	UpsPeak.DrawLatex(0.35,label_y-4.5*label_text_size,"ATLAS points normalized to ");
	UpsPeak.DrawLatex(0.35,label_y-5.5*label_text_size,"CMS 40<p_{T}(#mu#mu)<50 GeV");
	UpsPeak.DrawLatex(0.35,label_y-6.7*label_text_size,lumi_exc);
	
	CName["atlas_canvas"]->cd(ups);
	gPad->SetLogy();
	gr_atlas->GetXaxis()->SetRangeUser(0,100);
	xs->GetXaxis()->SetRangeUser(0,100); 
	xs->DrawClone("AP");
	xs_stat->DrawClone("|| same"); 
	gr_atlas->DrawClone("P same");
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	
	
	output_files->cd("comparisons"); 
	comp->Write();
	comp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_atlas_comp%dS.pdf",ups)); 
	comp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_atlas_comp%dS.jpg",ups)); 

	gr_atlas->Write();
	delete xs; 
	delete xs_stat; 
	delete gr_atlas;
}

void compare_cms2010(int ups){
	int iy=0;
	TGraphAsymmErrors *xs=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_graph_%dS",ups)); 
	TGraphAsymmErrors *xs_stat=(TGraphAsymmErrors*)output_files->FindObjectAny(Form("xs_stat_graph_%dS",ups));
	xs->GetXaxis()->SetTitle(x_label.c_str()); 
	xs->GetYaxis()->SetTitle(xs_y.c_str()); 
	xs->SetTitle("");
	TGraphAsymmErrors *cms2010=(TGraphAsymmErrors*)output_files->FindObjectAny(Form(Form("cms2010_%dS",ups)));
	
	double max=TMath::Max(cms2010->Eval(4.),xs->Eval(10.)); 
	xs->SetMaximum(max*1.4);
	
	//cms2010->GetXaxis()->SetRangeUser(0,100);
	//xs->GetXaxis()->SetRangeUser(0,100);

	set_CMS2011(xs,"",1);
	set_CMS2011(xs_stat,"stat",1);
	set_CMS2010(cms2010,"");
	
	TLegend *L = new TLegend(0.15,0.2,0.45,0.35,"","NDC");
	L->SetTextSize(label_text_size*0.9);
	L->SetFillColor(10); 
	L->SetLineColor(10);
	L->AddEntry(xs,"CMS 2011 |y(#mu#mu)|<0.6","LP"); 	
	//L->AddEntry(xs_stat,"Statistical Uncertainty","L"); 
	//L->AddEntry(xs,"Total Uncertainty","L"); 
	L->AddEntry(cms2010, "CMS 2010 |y(#mu#mu)|<2.4 #times0.25, 36 pb^{-1}","LP"); 
	
	//xs_stat->GetXaxis()->SetRangeUser(10,100); 
	
	TCanvas *C = new TCanvas(Form("Compare_cms2010_%dS",ups),"",800,800); 
	gPad->SetLogy();
	
	TMultiGraph *grM=new TMultiGraph(); 
	grM->Add(xs);
	grM->Add(xs_stat,"||");
	grM->Add(cms2010); 
	grM->SetMaximum(max*1.4); 
	grM->Draw("ap"); 
	grM->GetXaxis()->SetRangeUser(0,100); 
	set_titles(grM); 
	grM->Draw("ap");

	L->Draw("same"); 
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	TLatex UpsPeak(0.5,label_y-3*label_text_size,Form("#Upsilon(%dS)",ups));
	UpsPeak.SetNDC(kTRUE);
	UpsPeak.SetTextSize(label_text_size); 
	UpsPeak.Draw();
	UpsPeak.DrawLatex(0.4, label_y-4*label_text_size,lumi_exc);
	//UpsPeak.DrawLatex(0.5,0.65-label_text_size,rapidity_text(iy).c_str());	

	output_files->cd("comparisons");
	C->Write(); 
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_cms2010_comp%dS.pdf",ups));
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_cms2010_comp%dS.jpg",ups));

	delete xs; 
	delete xs_stat;
	delete cms2010; 
	delete grM; 

	delete C;
	
}

void compare_xs(int ups){
	
	double a[]={4.7289,.8699,.4249};
	double aE[]={0.4079,.0981,.0652};
	
	double b[]={0.3123,0.2544,0.2373};
	double bE[]={0.0126,0.0146,0.0215};
	
	double c[]={0.0016,0.0011,0.0010};
	double cE[]={0.0004,0.0004,0.0006};

	
	
	TF1 *xs_shape = new TF1("xs_shape",expo_pol2,10,50,3);
	xs_shape->SetParameters(a[ups-1],b[ups-1],c[ups-1]); 
	TF1 *xs_shapeP = new TF1("xs_shapeP",expo_pol2,10,50,3);
	xs_shapeP->SetParameters(a[ups-1]+aE[ups-1],b[ups-1]+bE[ups-1],c[ups-1]+cE[ups-1]); 
	TF1 *xs_shapeM = new TF1("xs_shapeM",expo_pol2,10,50,3);
	xs_shapeM->SetParameters(a[ups-1]-aE[ups-1],b[ups-1]-bE[ups-1],c[ups-1]-cE[ups-1]); 
	xs_shape->SetLineStyle(kDashed);
	xs_shapeP->SetLineColor(kRed);
	xs_shapeM->SetLineColor(kBlue);
	
	TH1D *xs=(TH1D*)output_files->FindObjectAny(Form("xs_%dS",ups)); 
		
	double scale=.25*1000;
	
	xs_shape->SetParameter(0,xs_shape->GetParameter(0)*scale); 
	xs_shapeP->SetParameter(0,xs_shapeP->GetParameter(0)*scale); 
	xs_shapeM->SetParameter(0,xs_shapeM->GetParameter(0)*scale); 
	
	
	TCanvas *C = new TCanvas("summary_comparison"); 
	gPad->SetLogy();
	xs->GetYaxis()->SetTitle(xs_y.c_str()); 
	double A=xs->Integral("width"); 
	//xs->Scale(1/A);
	xs->SetStats(kFALSE); 
	xs->SetTitle("");
	xs->Draw();
	xs_shape->Draw("same"); 
	xs_shapeP->Draw("same");
	xs_shapeM->Draw("same"); 
	
	C->Write();
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/xs/xs_2010bora_comp%dS.pdf",ups)); 
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/xs/xs_2010bora_comp%dS.jpg",ups)); 

	
	delete xs;
	delete C;
}

void print_tables(int iy, string mode){
		//print the yields
	
	cout << "Latex table " << endl;
	
	TH1D *y1 = (TH1D*)output_files->FindObjectAny(yield_histogram(iy,1,"actual").c_str());
	TH1D *y2 = (TH1D*)output_files->FindObjectAny(yield_histogram(iy,2,"actual").c_str());
	TH1D *y3 = (TH1D*)output_files->FindObjectAny(yield_histogram(iy,3,"actual").c_str());

	TH1D *yE1=(TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],1));
	TH1D *yE2=(TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],2));
	TH1D *yE3=(TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],3));
	
	yE1->Scale(1/100.);
	yE2->Scale(1/100.);
	yE3->Scale(1/100.);
	
	TH1D *sysP1 = (TH1D*)output_files->FindObjectAny(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],1));
	TH1D *sysM1 = (TH1D*)output_files->FindObjectAny(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],1));
	
	TH1D *sysP2 = (TH1D*)output_files->FindObjectAny(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],2));
	TH1D *sysM2 = (TH1D*)output_files->FindObjectAny(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],2));
	
	TH1D *sysP3 = (TH1D*)output_files->FindObjectAny(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],3));
	TH1D *sysM3 = (TH1D*)output_files->FindObjectAny(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],3));
	
	TH1D * A1=(TH1D*)output_files->FindObjectAny("acceptance_1S"); 
	TH1D * A2=(TH1D*)output_files->FindObjectAny("acceptance_2S"); 
	TH1D * A3=(TH1D*)output_files->FindObjectAny("acceptance_3S"); 
	
	TH1D * A1Ep=(TH1D*)output_files->FindObjectAny("acceptance_Ep_1S"); 
	TH1D * A2Ep=(TH1D*)output_files->FindObjectAny("acceptance_Ep_2S"); 
	TH1D * A3Ep=(TH1D*)output_files->FindObjectAny("acceptance_Ep_3S"); 
	
	TH1D * A1Em=(TH1D*)output_files->FindObjectAny("acceptance_Em_1S"); 
	TH1D * A2Em=(TH1D*)output_files->FindObjectAny("acceptance_Em_2S"); 
	TH1D * A3Em=(TH1D*)output_files->FindObjectAny("acceptance_Em_3S"); 
	
	TH1D *chi2_min = (TH1D*)output_files->FindObjectAny(Form("chi2_min_y%d",iy)); 
	TH1D *NDOF_min = (TH1D*)output_files->FindObjectAny(Form("NDOF_min_y%d",iy));
	
	//$<p_T>$ & Yield & $\mathcal{A}$ & Yield & $\mathcal{A}$ & Yield &
	//$\mathcal{A}$   \\ \hline
	
	string S=" $\\pm$ "; 
	
	ofstream output_table; 
	output_table.open("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/latex_table.txt");
	
	
	for (int ipt=1; ipt<=y1->GetNbinsX(); ipt++) {
		
		double y[]={y1->GetBinContent(ipt),y2->GetBinContent(ipt),y3->GetBinContent(ipt)};
		double yE[]={y[0]*yE1->GetBinContent(ipt),y[1]*yE2->GetBinContent(ipt),y[2]*yE3->GetBinContent(ipt)};

		double yEp[]={y[0]*sysP1->GetBinContent(ipt)/100,y[1]*sysP2->GetBinContent(ipt)/100,y[2]*sysP3->GetBinContent(ipt)/100};
		double yEm[]={y[0]*sysM1->GetBinContent(ipt)/100,y[1]*sysM2->GetBinContent(ipt)/100,y[2]*sysM3->GetBinContent(ipt)/100};
		double chi2=chi2_min->GetBinContent(ipt); 
		int NDOF=NDOF_min->GetBinContent(ipt); 
		
		int Nd=rint(y[0]/10); 
		int Nstat=rint(yE[0]/10);
		
		
		output_table << Form("%.1f",get_bin_center(ipt-1,1) ) << " & " << Form("%.1f",y1->GetBinContent(ipt)) << S<< Form("%.1f",yE[0]) << Form(" $^{+%.1f}_{-%.1f}$",yEp[0],yEm[0]);
		output_table << " & " << Form("%.1f",y2->GetBinContent(ipt) ) << S << Form("%.1f",yE[1] ) <<Form(" $^{+%.1f}_{-%.1f}$",yEp[1],yEm[1]);
		
		output_table << " & " << Form("%.1f",y3->GetBinContent(ipt) ) << S << Form("%.1f",yE[2] ) << Form(" $^{+%.1f}_{-%.1f}$",yEp[2],yEm[2]);
		output_table << " & " << Form("%.0f & %d",chi2, NDOF) << " \\\\ " << endl; 

		output_table << " & & & & &  \\\\" << endl; 
	}
	
	output_table << endl << endl << endl; 
	
	output_table << "xs table:  Units fb/GeV " << endl << endl << endl; 

	
	for (int ipt=1; ipt<=y1->GetNbinsX(); ipt++) {

		double sigma1=1000*((TH1D*)output_files->FindObjectAny(Form("xs_%dS%s",1,mode.c_str())))->GetBinContent(ipt); // xs 
		double sigma2=1000*((TH1D*)output_files->FindObjectAny(Form("xs_%dS%s",2,mode.c_str())))->GetBinContent(ipt);
		double sigma3=1000*((TH1D*)output_files->FindObjectAny(Form("xs_%dS%s",3,mode.c_str())))->GetBinContent(ipt);

		double sigma1E=sigma1*yE1->GetBinContent(ipt);
		double sigma2E=sigma2*yE2->GetBinContent(ipt);
		double sigma3E=sigma3*yE3->GetBinContent(ipt);
		
		double sigma1S[]={(sysP1->GetBinContent(ipt)/100)*sigma1,(sysM1->GetBinContent(ipt)/100)*sigma1}; // stat error +/- excluding lumi 
		double sigma2S[]={(sysP2->GetBinContent(ipt)/100)*sigma2,(sysM2->GetBinContent(ipt)/100)*sigma2};
		double sigma3S[]={(sysP3->GetBinContent(ipt)/100)*sigma3,(sysM3->GetBinContent(ipt)/100)*sigma3};

		
		string pt_bin=Form("%.0f-%.0f & ",fPTbin2[ipt-1],fPTbin2[ipt]);
		string ptm=Form("%.1f",get_bin_center(ipt-1,1) );
		string xs1=Form("%.1f",sigma1);
		string xs_stat1=Form("%.2f",sigma1E);
		string xs_sys1=Form(" $^{+%.2f}_{-%.2f}$",sigma1S[0],sigma1S[1]);
		
		string xs2=Form("%.1f",sigma2);
		string xs_stat2=Form("%.2f",sigma2E);
		string xs_sys2=Form(" $^{+%.2f}_{-%.2f}$",sigma2S[0],sigma2S[1]);

		
		string xs3=Form("%.1f",sigma3);
		string xs_stat3=Form("%.2f",sigma3E);
		string xs_sys3=Form(" $^{+%.2f}_{-%.2f}$",sigma3S[0],sigma3S[1]);

		
		output_table <<pt_bin << ptm << " & " << xs1 << S<< xs_stat1 << xs_sys1;
		output_table << " & " << xs2 << S << xs_stat2 << xs_sys2;
		output_table << " & " << xs3 << S << xs_stat3 << xs_sys3 << " \\\\" << endl;
		output_table << " & & & & \\\\ " << endl; 
	
	}
	
	output_table << "Acceptance Plot: " << endl << endl; 
	
	for (int ipt=1; ipt<=A1->GetNbinsX(); ipt++) {
		output_table << Form("%.1f",get_bin_center(ipt-1,1) ) << " & " << Form("%.2f",A1->GetBinContent(ipt) ) << S << Form("%.3f",A1->GetBinContent(ipt)-A1Ep->GetBinContent(ipt));
		double deltaA=100*(((TH1D*)output_files->FindObjectAny(Form("acceptance_reweighted_%dS",1)))->GetBinContent(ipt)-A1->GetBinContent(ipt))/A1->GetBinContent(ipt); 
		deltaA=TMath::Abs(deltaA);
		output_table << " & " << Form("%.2f",deltaA); 
		double A_p0=((TH1D*)output_files->FindObjectAny(Form("acceptance_nopol_%dS",1)))->GetBinContent(ipt);
		output_table << " & " << Form("%.2f",A_p0) << " & "; 
		
		output_table << Form("%.2f",A2->GetBinContent(ipt) ) << S << Form("%.3f & ",A2->GetBinContent(ipt)-A2Ep->GetBinContent(ipt));
		deltaA=100*(((TH1D*)output_files->FindObjectAny(Form("acceptance_reweighted_%dS",1)))->GetBinContent(ipt)-A2->GetBinContent(ipt))/A2->GetBinContent(ipt); 
		deltaA=TMath::Abs(deltaA);
	//	output_table << " & " << Form("%.2f",deltaA) << " "; 
		
		
		output_table << Form("%.2f",A3->GetBinContent(ipt) ) << S << Form("%.3f",A3->GetBinContent(ipt)-A3Ep->GetBinContent(ipt));
		deltaA=100*(((TH1D*)output_files->FindObjectAny(Form("acceptance_reweighted_%dS",1)))->GetBinContent(ipt)-A3->GetBinContent(ipt))/A3->GetBinContent(ipt); 
		deltaA=TMath::Abs(deltaA);
	//	output_table << " & " << Form("%.2f",deltaA) << " "; 
		
		output_table << " \\\\" << endl;
		//output_table << " & & & & \\\\ " << endl;
		
		
	}
	
	output_table << " cross section fits: " << endl << endl; 
	
	output_table << " Upsilon << A << C << alpha << chi2 << NDF " << endl; 
	
	for (int ups=1; ups<=3; ups++) {
		TF1 *xs_fit = (TF1*)output_files->FindObjectAny(Form("xs_fit_fcn_%dS",ups));
		output_table << "$\\" << Form("Upsilon(%dS)$ ",ups) << " & " << Form("%.2f %s %.3f",xs_fit->GetParameter(0),S.c_str(), xs_fit->GetParError(0)) << " & ";
		output_table << Form("%.2f %s %.3f",xs_fit->GetParameter(2),S.c_str(), xs_fit->GetParError(2)) << " & "; 
		output_table << Form("%.2f %s %.3f",xs_fit->GetParameter(1),S.c_str(), xs_fit->GetParError(1)) << " & "; 
		output_table << Form("%.0f",xs_fit->GetChisquare()) << " & ";
		output_table << Form("%d",xs_fit->GetNDF())  << "\\\\" << endl; 
		delete xs_fit;

	}
	output_table << endl << endl; 
	
	output_table.close();
	
	delete y1; 
	delete y2;
	delete y3; 
	
	delete yE1; 
	delete yE2;
	delete yE3;

	delete sysP1;
	delete sysM1; 
	delete sysP2;
	delete sysM2; 
	delete sysP3;
	delete sysM3; 
	
	delete A1; 
	delete A2;
	delete A3; 
	
	delete A1Em; 
	delete A1Ep; 
	delete A2Ep;
	delete A2Em; 
	delete A3Ep; 
	delete A3Em;
	

}

double expo_pol2(double *x, double *par){
	return par[0]*TMath::Exp(-par[1]*x[0]+par[2]*x[0]*x[0]);
}

double power_law(double *x, double*par){
	return par[0]/(par[2]+TMath::Power(x[0]/20,par[1])); 
}

double xs_full_func(double *x, double *par){
	double Exp=TMath::Exp(par[0]+par[1]*x[0]);
	double power=par[2]/(par[3]+TMath::Power(x[0]/par[4],par[5])); 
	
	if(x[0]<par[4]) return Exp; 
	else return power; 

}

void compute_yields(){
	CName["fits"]->Divide(4,4);
	CName["fits2"]->Divide(3,3); 
	int iy=0; 
	TH1D *chi2_min = new TH1D(Form("chi2_min_y%d",iy), "Minimum #Chi^{2}", fNpt2, fPTbin2); 
	TH1D *NDOF_min = new TH1D(Form("NDOF_min_y%d",iy), "NDOF for Minimum #Chi^{2}", fNpt2, fPTbin2); 
		
	for(int ipt=0; ipt<fNpt2; ipt++){
		//if(ipt!=0) continue;
		cout << "Pt bin: " << ipt << " PT: " << fPTbin2[ipt] << endl; 
		cout << "get dm_m distribution." << endl;
		get_dm_m(0,ipt,""); // For efficiency weighting
		cout << "get mass distribution." << endl; 
		make_mdata(0,ipt);
		cout << "Get lineshape." << endl; 
		get_LS(0, ipt, "W");
		fit(0,ipt);
		//make_fit_plot(0,ipt,""); 
		double chi2=0; 
		int NDOF=0; 
		compute_mass_sys(0, ipt,chi2,NDOF);
		chi2_min->SetBinContent(ipt+1, chi2); 
		NDOF_min->SetBinContent(ipt+1, NDOF); 
		//fit_sys_error(0, ipt);
		cout << endl << endl << endl; 
	}
	
	chi2_min->Write();
	NDOF_min->Write();
	
	cout << "Fill Toy Yields: " << endl; 
	//fill_yields(0,""); 
	fill_toy_yield(1,0); 
	fill_toy_yield(2,0); 
	fill_toy_yield(3,0); 
	fill_ratio_yield(0,3);
	fill_ratio_yield(0,2);
	
	cout << "Make fit plots: " << endl; 
	
	for(int ipt=0; ipt<fNpt2; ipt++){
		make_fit_plot(0,ipt,"");

	}
	
	yield_GeV(0,"");
	yield_GeV(0,"fit_sys");

	cout << "Write Histos: " << endl; 
	write_histos();
	//cout << "Draw Bkg errors: " << endl; 
	//draw_background_errors();
	
	
}

double compute_resid(int iy, int ipt, TH1D *h, TF1 *fit_shape, int &NDOF){
	//Compute residuals and return chi2 for entire range 
	
	//TH1D *resid=(TH1D*)mdata->Clone("resid");
	//resid->GetYaxis()->SetTitle("Residuals"); 
	
	//cout << "hbins: " << h->GetNbinsX() << endl; 
	//cout << "eval: " << fit_shape->Eval(9.4603) << endl; 
	
	
	TH1D *resid = new TH1D("resid","pull;M_{#mu#mu} GeV; pull=(Fit-mdata)/#sigma",h->GetNbinsX(),8.7,11.2);
	
	TH1D *chi_hist =new TH1D("chi_hist","Histogram of pull;pull,Number of events",100,-15,15); 
	
	double chi=0; // chi for eachbin 
	double chi2_tot=0;
	double error=0; 
	
	bool sideband=false; 
	double M1; 
	double M2;
	
	if(fit_shape->GetNpar()==7){ 
		M1=fit_shape->GetParameter(5);
		M2=fit_shape->GetParameter(6); 
	}
	else {
		M1=fit_shape->GetParameter(10);
		M2=fit_shape->GetParameter(11);
	}

	if(M1>0 || M2>0) sideband=true;
	//cout << "M1: " << M1 << " M2: " << M2 << endl; 
	
	NDOF=fit_shape->GetNDF();
	
	//cout <<"NDOF: " << NDOF << endl; 	
	
	for(int i=h->FindBin(8.7); i<h->FindBin(11.2); i++)	{
		error=h->GetBinError(i); //Get bin errors from data 
		if(h->GetBinContent(i)<1){
			error=TMath::Sqrt(h->GetBinContent(i)+1); 
		}
		
		chi=(fit_shape->Eval(h->GetBinCenter(i))-h->GetBinContent(i))/error;
		if(h->GetBinContent(i)<5 && sideband==false){
				chi=0; 
				NDOF-=1; 
		}//SB false
		if(sideband==true){
			if (h->GetBinCenter(i)>M1 && h->GetBinCenter(i)<M2) {
				chi=0; 
			}
			else if(h->GetBinContent(i)<5){
				//cout <<"Empty bin:" << h->GetBinCenter(i)<< " Bin Content: " << h->GetBinContent(i) <<  endl; 
				chi=0; 
				if(NDOF>1)NDOF-=1; 
			}
		}//SB true
				
		chi_hist->Fill(chi); 
		chi2_tot+=chi*chi; // chi2
		resid->SetBinContent(i,h->GetBinCenter(i),chi);
		resid->SetBinError(i,h->GetBinCenter(i),error);
		
	}
	
	//cout << "Findal NDOF: " << NDOF << endl; 
	
	//Return total Chi2. 
	
	TCanvas *CResid=new TCanvas(Form("Residuals_pt%.1f-%.1f",fPTbin2[ipt],fPTbin2[ipt+1]),"",800,800);
	CResid->Divide(2,2);
	
	TPaveText *txt = new TPaveText(0.1,0.1,0.9,0.9); 
	
	
	txt->AddText(Form("#chi^{2}=%.1f",chi2_tot)); 
	txt->AddText(Form("NDOF=%d",NDOF)); 
	txt->AddText(Form("#chi^{2}/NDOF=%.1f",chi2_tot/static_cast<double>(NDOF)));
	txt->AddText(Form("Mean #chi^{2}=%.1f +/- %.2f",chi_hist->GetMean(),chi_hist->GetMeanError()));
	
	double up=TMath::Abs(resid->GetBinContent(resid->GetMaximumBin())+resid->GetBinError(resid->GetMaximumBin()));
	double down=TMath::Abs(resid->GetBinContent(resid->GetMinimumBin())) +TMath::Abs(resid->GetBinError(resid->GetMinimumBin())); 
	
	double range=TMath::Max(up,down); 
	
	range=range+range*0.1;
	
	TF1	 *mean = new TF1("mean","pol0",8.7,11.2); 
	mean->SetParameter(0,chi_hist->GetMean()); 
	mean->SetLineColor(kRed); 
	
	resid->SetAxisRange(-range,range,"Y");
	resid->SetStats(kFALSE);
	CResid->cd(1);
	resid->Draw("E1");
	mean->Draw("same"); 
	txt->DrawClone("same"); 
	CResid->cd(2);
	chi_hist->Draw(); 
	CResid->cd(3);
	txt->DrawClone("same"); 
	
	output_files->cd("mass_shape_pull");
	if(!sideband) CResid->Write();
	
	delete resid; 
	delete CResid;
	delete chi_hist;
	delete mean; 
	delete txt;
	
	return chi2_tot; 
}

void set_name_all(TH1D *h){
	h->SetTitle("N_{gen}(|y(#mu#mu)|<0.6)"); 
	h->GetYaxis()->SetTitle("Events"); 
	h->GetXaxis()->SetTitle("p_{T}(#mu#mu) [GeV]");
	
}

void set_name_reco(TH1D *h){
	h->SetTitle("N_{gen}(|#eta^{#mu}|<0.9 && P^{#mu}_{T}>4.5 && Muons reconstructed)"); 
	h->GetYaxis()->SetTitle("Events"); 
	h->GetXaxis()->SetTitle("p_{T}(#mu#mu) [GeV]");
}

void set_acceptance_info(TH1D *h,int ups){

	h->GetYaxis()->SetTitle(Form("A(%dS)",ups));
	h->SetStats(kFALSE);
	h->SetTitle(""); //Acceptance |y(#mu#mu)|<0.5
}

void acceptance_2D(int ups){
	TH2D *den=(TH2D*)gDirectory->FindObjectAny(Form("UG_AllGenRes_y_pt_%dS",ups));
	TH2D *A2d=(TH2D*)gDirectory->FindObjectAny(Form("UG_RecoGenRes_y_pt_%dS",ups));
	A2d->Divide(den); 
	TCanvas *A2DCan=new TCanvas("","",800,800); 
	A2d->SetStats(kFALSE);
	A2d->GetXaxis()->SetTitle("p_{T}(#mu#mu) [GeV]"); 
	A2d->GetYaxis()->SetTitle(Form("A(%dS)",ups)); 
	A2d->Draw("COLZ"); 
	A2DCan->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/acceptance/A2d_%dS.pdf",ups));

	delete A2DCan; 
	delete A2d;
	delete den; 
}

void acceptance(int ups){
	gROOT->SetStyle("Plain");
	cout << "Compute acceptance: " << ups << endl; 
	if(ups==1) acceptance_file1S->cd("acceptance");
	if(ups==2) acceptance_file2S->cd("acceptance");
	if(ups==3) acceptance_file3S->cd("acceptance"); 
	
	acceptance_2D(ups);
	
	TH1D *allnoW=(TH1D*)gDirectory->Get(Form("UG_AllGenResnoW_%dS",ups));//no weights

	TH1D *all=(TH1D*)gDirectory->Get(Form("UG_AllGenRes_%dS",ups));
	set_name_all(all);
	TH1D *all_Ep=(TH1D*)gDirectory->Get(Form("UG_AllGenRes_Ep_%dS",ups));
	set_name_all(all_Ep); 
	TH1D *all_Em=(TH1D*)gDirectory->Get(Form("UG_AllGenRes_Em_%dS",ups));
	set_name_all(all_Em); 
	
	TH1D *reconoW=(TH1D*)gDirectory->Get(Form("UG_RecoGenResnoW_%dS",ups));
	TH1D *reco=(TH1D*)gDirectory->Get(Form("UG_RecoGenRes_%dS",ups));
	set_name_reco(reco);
	TH1D *reco_Ep=(TH1D*)gDirectory->Get(Form("UG_RecoGenRes_Ep_%dS",ups));
	set_name_reco(reco_Ep); 
	TH1D *reco_Em=(TH1D*)gDirectory->Get(Form("UG_RecoGenRes_Em_%dS",ups));
	set_name_reco(reco_Em); 
	
	//cout << "contents: " << all->GetBinContent(1) << " " << reco->GetBinContent(1) << endl; 
	//cout << "error: " << all->GetBinError(1) << " " << reco->GetBinError(1) << endl; 
	
	double alpha=reco->GetBinError(1)/all->GetBinError(1);
	double err=(1/(TMath::Sqrt(alpha+1)*alpha))*(1/TMath::Sqrt(all->GetBinContent(1)));
	double A=reco->GetBinContent(1)/all->GetBinContent(1); 
	double AEcalc=A*TMath::Sqrt(TMath::Power((reco->GetBinError(1)/reco->GetBinContent(1)),2)+TMath::Power((all->GetBinError(1)/all->GetBinContent(1)),2));
	
	
	TH1D *acceptance=(TH1D*)reco->Clone();
	TH1D *acceptancenoW=(TH1D*)reconoW->Clone();
	TH1D *acceptance_Ep=(TH1D*)reco_Ep->Clone();
	TH1D *acceptance_Em=(TH1D*)reco_Em->Clone(); 
	
	set_acceptance_info(acceptance,ups);
	set_acceptance_info(acceptance_Ep,ups); 
	set_acceptance_info(acceptance_Em,ups); 
	
	acceptance->Divide(all);
	acceptance_Ep->Divide(all_Ep); 
	acceptance_Em->Divide(all_Em); 
	acceptancenoW->Divide(allnoW); 
	
	acceptance_Ep->SetMarkerColor(kRed); 
	acceptance_Ep->SetLineColor(kRed); 
	acceptance_Em->SetMarkerColor(kGreen);
	acceptance_Em->SetLineColor(kGreen); 
	acceptancenoW->SetLineColor(kBlue);
	
	TH1D *error=(TH1D*)reco->Clone();
	TH1D *error_Ep=(TH1D*)reco->Clone();
	TH1D *error_Em=(TH1D*)reco->Clone();
	
	error->SetStats(kFALSE);
	error->SetTitle("Acceptance Percent Error"); 
	error->GetYaxis()->SetTitle("Error [%]"); 
	
	error_Ep->SetStats(kFALSE);
	error_Ep->SetTitle("Acceptance Percent Error"); 
	error_Ep->GetYaxis()->SetTitle("Error [%]"); 
	
	error_Em->SetStats(kFALSE);
	error_Em->SetTitle("Acceptance Percent Error"); 
	error_Em->GetYaxis()->SetTitle("Error [%]"); 
	
	
	error->SetMarkerStyle(7); 
	error_Ep->SetMarkerStyle(7); 
	error_Em->SetMarkerStyle(7); 
	
	error_Ep->SetLineColor(kRed); 
	error_Ep->SetMarkerColor(kRed); 
	error_Em->SetLineColor(kGreen); 
	error_Em->SetMarkerColor(kGreen);

	cout << "Accpetance(" << ups << "S)" << endl; 
	
	for (int i=1; i<=acceptance->GetNbinsX();i++ ) {
		double percent_error=acceptance->GetBinError(i)*100/acceptance->GetBinContent(i);
		double percent_error_Ep=TMath::Abs((acceptance_Ep->GetBinContent(i)-acceptance->GetBinContent(i))*100/acceptance->GetBinContent(i));
		double percent_error_Em=TMath::Abs((acceptance->GetBinContent(i)-acceptance_Em->GetBinContent(i))*100/acceptance->GetBinContent(i));
		
		//cout << "ipT: " << i << " Ep: " << percent_error_Ep << " Em: " << percent_error_Em << endl; 
		
		error->SetBinContent(i,percent_error);
		error_Ep->SetBinContent(i,percent_error_Ep); 
		error_Em->SetBinContent(i,percent_error_Em); 
		
		error_Ep->SetBinError(i,0);
		error_Em->SetBinError(i,0);
		acceptance->SetBinError(i,0); 
		acceptance_Em->SetBinError(i,0);
		acceptance_Ep->SetBinError(i,0);
		acceptancenoW->SetBinError(i,0); 
		
	}
	
	TCanvas *C=new TCanvas(Form("acceptance_summary_canvas_%dS",ups),"",800,800);
	C->Divide(2,2);
	C->cd(1);
	gPad->SetLogy();
	all->Draw();
	C->cd(2);
	gPad->SetLogy();
	reco->Draw();
	C->cd(3); 
	acceptance->SetMinimum(0);
	acceptance->SetMaximum(1);
	acceptance->DrawCopy();
	acceptance_Ep->DrawCopy("same"); 
	acceptance_Em->DrawCopy("same"); 
	acceptancenoW->DrawCopy("same");
	
	C->cd(4);
	error->SetMinimum(0);
	error->SetMaximum(5); 
	error->Draw("histo P");
	error_Ep->Draw("same histo P"); 
	error_Em->Draw("same histo P"); 
	
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/acceptance/acceptance_%dS.pdf",ups));
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/acceptance/acceptance_%dS.jpg",ups));

	
	TCanvas *CC = new TCanvas(Form("acceptance_canvas_%dS",ups),"",800,800);
	

	output_files->cd("acceptance");
	
	cout << "Writing Acceptance Histograms." << endl; 
	
	C->Write();
	CC->Write();
	
	acceptance->SetName(Form("acceptance_1GeV%dS",ups)); 
	acceptance_Ep->SetName(Form("acceptance_1GeV_Ep_%dS",ups)); 
	acceptance_Em->SetName(Form("acceptance_1GeV_Em_%dS",ups)); 
	
	error_Ep->SetName(Form("acceptance_Error_1GeV_Ep_%dS",ups));
	error_Em->SetName(Form("acceptance_Error_1GeV_Em_%dS",ups)); 
	
	acceptance->Write(); 
	acceptance_Em->Write();
	acceptance_Ep->Write();
	
	TH1D *aR=(TH1D*)acceptance->Rebin(fNpt2, Form("acceptance_%dS",ups), fPTbin2); 
	TH1D *aR0=(TH1D*)acceptancenoW->Rebin(fNpt2, Form("acceptance_nopol_%dS",ups), fPTbin2);
	TH1D *aR_Ep=(TH1D*)acceptance_Ep->Rebin(fNpt2, Form("acceptance_Ep_%dS",ups), fPTbin2); 
	TH1D *aR_Em=(TH1D*)acceptance_Em->Rebin(fNpt2, Form("acceptance_Em_%dS",ups), fPTbin2); 

	aR->SetName(Form("acceptance_%dS",ups)); 
	aR_Ep->SetName(Form("acceptance_Ep_%dS",ups)); 
	aR_Em->SetName(Form("acceptance_Em_%dS",ups)); 
	
	aR0->Divide(bin_width);
	aR->Divide(bin_width);
	aR_Ep->Divide(bin_width);
	aR_Em->Divide(bin_width);
	
	TH1D *errorR_Ep=(TH1D*)acceptance_Ep->Rebin(fNpt2, Form("acceptance_error_Ep_%dS",ups), fPTbin2); 
	TH1D *errorR_Em=(TH1D*)acceptance_Em->Rebin(fNpt2, Form("acceptance_error_Em_%dS",ups), fPTbin2); 
	errorR_Ep->SetName( Form("acceptance_error_Ep_%dS",ups)); 
	errorR_Em->SetName(Form("acceptance_error_Em_%dS",ups)); 
	
	cout << "rebinned: " << endl;
	
	for (int i=1; i<=aR->GetNbinsX();i++ ) {
		double percent_error_Ep=TMath::Abs((aR_Ep->GetBinContent(i)-aR->GetBinContent(i))*100/aR->GetBinContent(i));
		double percent_error_Em=TMath::Abs((aR_Em->GetBinContent(i)-aR->GetBinContent(i))*100/aR->GetBinContent(i)); 
		errorR_Ep->SetBinContent(i,percent_error_Ep); 
		errorR_Em->SetBinContent(i,percent_error_Em); 
		
		//cout << "pT: " << aR->GetBinCenter(i) << "A: " << aR->GetBinContent(i) << " Ep: " << percent_error_Ep << " Em: " << percent_error_Em << endl; 
		aR->SetBinError(i,0);
		aR0->SetBinError(i,0);
		errorR_Ep->SetBinError(i,0);
		errorR_Em->SetBinError(i,0);
		
	}

	aR->Write(); 
	aR0->Write();
	aR_Em->Write();
	aR_Ep->Write();
	errorR_Ep->Write();
	errorR_Em->Write();
	
	CC->cd();
	aR->Draw();
	aR_Ep->Draw("same");
	aR_Em->Draw("same");
	prelimTextSim.SetNDC(kTRUE);
	prelimTextSim.Draw();
	
	CC->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/acceptance/acceptance_only_%dS.pdf",ups)); 
	CC->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/acceptance/acceptance_only_%dS.jpg",ups)); 

	
	error_Ep->Write();
	error_Em->Write(); 
	
	delete acceptance;
	delete acceptance_Ep;
	delete acceptance_Em;
	
	delete aR;
	delete aR_Em;
	delete aR_Ep;
	
	delete error;
	delete error_Ep;
	delete error_Em;
	
	delete C; 
	delete CC;
}

void acceptance_summary_plot(){
	TH1D *A1=(TH1D*)output_files->FindObjectAny(Form("acceptance_%dS",1)); 
	TH1D *A1Ep=(TH1D*)output_files->FindObjectAny(Form("acceptance_Ep_%dS",1));
	TH1D *A1Em=(TH1D*)output_files->FindObjectAny(Form("acceptance_Em_%dS",1));
	
	TH1D *A2=(TH1D*)output_files->FindObjectAny(Form("acceptance_%dS",2)); 
	TH1D *A2Ep=(TH1D*)output_files->FindObjectAny(Form("acceptance_Ep_%dS",2));
	TH1D *A2Em=(TH1D*)output_files->FindObjectAny(Form("acceptance_Em_%dS",2));
	
	TH1D *A3=(TH1D*)output_files->FindObjectAny(Form("acceptance_%dS",3)); 
	TH1D *A3Ep=(TH1D*)output_files->FindObjectAny(Form("acceptance_Ep_%dS",3));
	TH1D *A3Em=(TH1D*)output_files->FindObjectAny(Form("acceptance_Em_%dS",3));
	
	TCanvas *summary = new TCanvas("summary3_acceptance");
	
	summary->Divide(2,2);
	
	
	TLegend * L = new TLegend(0.25,0.5,0.9,0.9); 
	
	L->AddEntry(A1); 
	L->AddEntry(A1Ep); 
	L->AddEntry(A1Em); 
	
	summary->cd(1);	
	A1Em->Draw("histo P");
	A1Ep->Draw("histo P same");
	A1->Draw("histo P same"); 
	
	summary->cd(2);
	A2Em->Draw("histo P");
	A2Ep->Draw("histo P same");
	A2->Draw("histo P same");

	summary->cd(3);
	A3Em->Draw("histo P");
	A3Ep->Draw("histo P same"); 
	A3->Draw("histo P same");
	
	summary->cd(4);
	L->Draw();
	
	summary->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/acceptance/all3_accept.pdf");
	summary->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/acceptance/all3_accept.jpg");

	output_files->cd("acceptance");
	cout << "Writing summary acceptance plots. " << endl;
	summary->Write();
	delete summary; 

	
}

void setBW(){
	for (int i=0; i<fNpt2; i++) {
		bin_width->SetBinContent(i+1,fPTbin2[i+1]-fPTbin2[i]);
		bin_width->SetBinError(i+1,0);
	}	
}

void compare_bkg_systematics(int iy, int ups){
	
	cout << "Comparison of Systematics for BKG1 and BKG2. " << endl; 
	cout << "In Method 1 we fix deltaM and dm_scale and in method 2 we leave them floating. " << endl;  
	
	for (int ipt=1; ipt<=fNpt2; ipt++) {
		double y=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinContent(ipt);
		
		double pt=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinCenter(ipt);
		
		double BKG1=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"BKGP_SF").c_str()))->GetBinContent(ipt); 
		double deltaM1=((TF1*)output_files->FindObjectAny(fit_function(iy,ipt-1,"BKGP_SF").c_str()))->GetParameter("deltaM"); 	
		double dm_scale1=((TF1*)output_files->FindObjectAny(fit_function(iy,ipt-1,"BKGP_SF").c_str()))->GetParameter("dm_scale"); 	
		
		double BKG2=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"BKGP").c_str()))->GetBinContent(ipt); 
		double deltaM2=((TF1*)output_files->FindObjectAny(fit_function(iy,ipt-1,"BKGP").c_str()))->GetParameter("deltaM"); 	
		double dm_scale2=((TF1*)output_files->FindObjectAny(fit_function(iy,ipt-1,"BKGP").c_str()))->GetParameter("dm_scale"); 
		
		cout << setprecision(4) << " pt: " << pt << " sys1: " << BKG1-y << " deltaM: " << deltaM1*1000 << "MeV dm_scale: " << dm_scale1 << " sys2: "<< BKG2-y << " deltaM: " << deltaM2*1000 << "MeV dm_scale: " << dm_scale2 << endl;
		
	}
}

void plot_total_error(int iy,int ups){

		
	TH1D *stat=(TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups));
	TH1D *sysP=(TH1D*)output_files->FindObjectAny(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups));
	TH1D *sysM=(TH1D*)output_files->FindObjectAny(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups));
	
	TH1D *AEp=(TH1D*)output_files->FindObjectAny(Form("acceptance_error_Ep_%dS",ups));
	TH1D *AEm=(TH1D*)output_files->FindObjectAny(Form("acceptance_error_Em_%dS",ups)); 

	TH1D *EffP=(TH1D*)output_files->FindObjectAny(Form("EffP_per_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups)); 
	TH1D *EffM=(TH1D*)output_files->FindObjectAny(Form("EffM_per_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups));
	
	TH1D *Fit_sys = (TH1D*)output_files->FindObjectAny(Form("Fit_per_y%d_%dS",iy,ups));
	
	TH1D *rho = (TH1D*)output_files->FindObjectAny(Form("rho_per_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups)); 
	
	AEp->SetLineColor(kBlue); 
	AEm->SetLineColor(kBlue);
	AEm->SetLineStyle(kDashed); 
	
	double SM[]={stat->GetBinContent(stat->GetMaximumBin()),sysP->GetBinContent(sysP->GetMaximumBin()),sysM->GetBinContent(sysM->GetMaximumBin())};
	double max=TMath::MaxElement(3,SM); 
	
	stat->SetMaximum(max*1.1); 
	stat->SetTitle(Form("#Upsilon(%dS)",ups));
	stat->SetStats(kFALSE);
	
	CName["error_summary"]->cd(ups);
	
	stat->GetXaxis()->SetTitle("p_{T}(#mu#mu) [GeV]"); 
	stat->GetYaxis()->SetTitle("Error [%]"); 
	
	stat->SetLineColor(kRed); 
	stat->DrawCopy("");
	
	sysP->SetLineColor(kBlack); 
	sysP->DrawCopy("same"); 
	
	sysM->SetLineColor(kBlack); 
	sysM->SetLineStyle(kDashed); 
//	sysM->DrawCopy("same"); 
	
	rho->SetLineColor(kOrange); 
	rho->DrawCopy("same"); 
	
	//EffM->Draw("same"); 

	EffP->SetLineStyle(kDashed); 
	EffP->SetLineColor(kGreen);
	EffP->DrawCopy("same");
	
	Fit_sys->SetLineColor(kRed);
	Fit_sys->SetLineStyle(kDotted);
	Fit_sys->DrawCopy("same"); 
	
	AEp->SetLineColor(kBlue);
	AEp->DrawCopy("same");
	//AEm->DrawCopy("same"); 
	
	
	TH1D *Lumi = new TH1D("Lumi", "Lumi [%]", fNpt2,fPTbin2); 
	for (int ipt=1; ipt<=Lumi->GetNbinsX(); ipt++) {
		Lumi->SetBinContent(ipt,4); 
		Lumi->SetBinError(ipt,0); 
	}
	Lumi->SetLineColor(kBlue); 
	
	if(ups==1){
		error_legend->SetFillColor(10);
		error_legend->SetLineColor(10);
		error_legend->SetShadowColor(10);
		
		error_legend->AddEntry((TH1D*)stat->Clone("S"),"Statistical Uncertainty","L");
		error_legend->AddEntry((TH1D*)sysP->Clone("SYSP"),"Total Sys","L"); 
		//error_legend->AddEntry((TH1D*)sysM->Clone("SYSM"),"Total Sys-"); 
		
		error_legend->AddEntry((TH1D*)EffP->Clone("EFFP"),"#epsilon_{#mu} Sys","L"); 
		//error_legend->AddEntry((TH1D*)EffM->Clone("EFFM"),"Eff Sys-"); 
		
		error_legend->AddEntry((TH1D*)Fit_sys->Clone("Fit"),"Fit Sys","L"); 
		
//		error_legend->AddEntry((TH1D*)Lumi->Clone(),"Lumi Sys Error"); 
		AEp->SetLineColor(kBlue);
		error_legend->AddEntry((TH1D*)Lumi->Clone(),"Polarization Sys","L"); 
		//error_legend->AddEntry((TH1D*)AEm->Clone(),"Acceptance sys-"); 
		error_legend->AddEntry((TH1D*)rho->Clone(),"#rho sys","L"); 
	}
	
	
	//Lumi->DrawCopy("same");

	delete Lumi;
	delete stat;
	delete sysP; 
	delete sysM; 
	delete Fit_sys;
	delete AEm;
	delete AEp; 
	delete rho; 

}

void stat_error(int iy, int ups){

	TH1D *stat=new TH1D(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups), "Total Statistical Error [%]", fNpt2,fPTbin2); 
	for(int ipt=1; ipt<=stat->GetNbinsX(); ipt++){
		double y=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinContent(ipt);
		double Err=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinError(ipt);
		stat->SetBinContent(ipt,100*Err/y); 
		}
	output_files->cd("yield_histograms"); 
	stat->Write(); 
	delete stat; 

}

void total_systematic(int iy, int ups){
	
	ofstream output_table; 
	output_table.open(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/systematic_table_%d.txt",ups));

	TFile *output_filesEp=new TFile("output_feb11Ep.root","READ");
	TFile *output_filesEm=new TFile("output_feb11Em.root","READ");

	TH1D *total_sysP_per=new TH1D(Form("total_sysP_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups),"Total Systematic Positive [%] total yield", fNpt2, fPTbin2); 
	TH1D *total_sysM_per=new TH1D(Form("total_sysM_y%.1f_%.1f_%dS_percent",fYbin2[iy],fYbin2[iy+1],ups),"Total Systematic Negative [%] total yield", fNpt2, fPTbin2);

	TH1D *Fitsys= new TH1D(Form("Fit_per_y%d_%dS",iy,ups), Form("Fit Sys #Upsilon(%dS)",ups), fNpt2, fPTbin2); 
	
	TH1D *EffP_per = new TH1D(Form("EffP_per_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups),"Shape syst+", fNpt2, fPTbin2); 
	TH1D *EffM_per = new TH1D(Form("EffM_per_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups),"Shape syst-", fNpt2, fPTbin2); 

	TH1D *rhoP_per = new TH1D(Form("rho_per_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups),"#rho sys+", fNpt2, fPTbin2); 
	
	TCanvas * C = new TCanvas(Form("Fit_systematics_summary_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups),"",800,800);
	C->Divide(1,2); 
	
	TFile *SG_file = new TFile("/Users/carlsonbt1/results/data_files/SG_efficiency.root","READ"); 
	TF1 *SG_eff=(TF1*)SG_file->Get("SG_eff"); 
	output_files->cd();
	
	TLegend *L = new TLegend(0.1,0.1,0.9,0.9); 
	
	cout << "Filling Systematic Histogram: " << endl; 
	
	TH1D *rho_uncertaintyP = new TH1D("rho_uncertaintyP", "#rho Uncertainty percent", fNpt2, fPTbin2); 
	TH1D *rho_uncertaintyM = new TH1D("rho_uncertaintyM", "#rho Uncertainty percent", fNpt2, fPTbin2); 
	output_table << "Y(" << ups << "S)" << endl; 
	output_table << "pT: " << " fit: " << " Eff: " << " rho: " << " A: " << " vtx: " << " sg: " << endl; 
	
	for (int ipt=1; ipt<=rho_uncertaintyP->GetNbinsX(); ipt++) {
		double pT=rho_uncertaintyP->GetBinCenter(ipt);
		
		if(pT<50) {
			rho_uncertaintyP->SetBinContent(ipt,2); 
			rho_uncertaintyM->SetBinContent(ipt,2);
		}
		if(pT>50 && pT<=59)rho_uncertaintyM->SetBinContent(ipt,2.5); 
		if(pT>59 && pT<=70) rho_uncertaintyM->SetBinContent(ipt,6); 
		if(pT>70 && pT<=100) rho_uncertaintyM->SetBinContent(ipt,15); 
		
		if(pT>50 && pT<=59)rho_uncertaintyP->SetBinContent(ipt,2.5); 
		if(pT>59 && pT<=70) rho_uncertaintyP->SetBinContent(ipt,6); 
		if(pT>70 && pT<=100) rho_uncertaintyP->SetBinContent(ipt,15); 

	} 
	
	for (int ipt=1; ipt<=total_sysP_per->GetNbinsX();ipt++) {
		
		double y=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinContent(ipt); 
		double yE=((TH1D*)output_files->FindObjectAny(Form("stat_error_y%.1f_%.1f_%dS",fYbin2[iy],fYbin2[iy+1],ups)))->GetBinContent(ipt);
		double EEp=0; // systematic error plus 
		double EEm=0; // error minus
		
		//fit systematic (symmetric) due to method
		double FitP=100*((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"fit_sys").c_str()))->GetBinContent(ipt)/y; 
		double FitM=FitP; 
		
		output_table << total_sysP_per->GetBinCenter(ipt) << " & "; 
		
		output_table << Form("%.2f & ",FitP);
		
		Fitsys->SetBinContent(ipt,FitP); 
		Fitsys->SetBinError(ipt,0); 
		
		//cout << "FitP: " << FitP << " " << Fitsys->GetBinContent(ipt) << endl; 
		
		EEp+=TMath::Power(FitP,2); 
		EEm+=TMath::Power(FitM,2); 
		
        //efficiency error 		
		double EffP=100*(((TH1D*)output_filesEp->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinContent(ipt)-y)/y;
		double EffM=100*(((TH1D*)output_filesEm->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinContent(ipt)-y)/y;

		double rho_SysP=rho_uncertaintyP->GetBinContent(ipt);
		double rho_SysM=rho_uncertaintyM->GetBinContent(ipt);
		
		output_table << Form("%.1f(%.1f)",rho_SysP, rho_SysM) << " & ";
		
		EEp+=TMath::Power(rho_SysP,2);
		EEm+=TMath::Power(rho_SysM,2);
		
		double EffP_rho=TMath::Power(EffP*EffP+rho_SysP*rho_SysP,2); 
		double EffM_rho=TMath::Power(EffM*EffM+rho_SysM*rho_SysM,2); 

		rhoP_per->SetBinContent(ipt,rho_SysP); 
		
		EffP_per->SetBinContent(ipt,TMath::Abs(EffP) ); 
		EffM_per->SetBinContent(ipt,TMath::Abs(EffM) ); 
		
		EEp+=TMath::Power(EffP,2);
		EEm+=TMath::Power(EffM,2); 
		
		output_table << Form("%.1f(%.1f)",TMath::Abs(EffP), TMath::Abs(EffM)) << " & ";

		//Acceptance Error [%] already
		double Ap=((TH1D*)output_files->FindObjectAny(Form("acceptance_error_Ep_%dS",ups)))->GetBinContent(ipt);
		double Am=((TH1D*)output_files->FindObjectAny(Form("acceptance_error_Em_%dS",ups)))->GetBinContent(ipt);
		
		EEp+=TMath::Power(Ap,2); 
		EEm+=TMath::Power(Am,2); 
		
		output_table << Form("%.1f(%.1f)",TMath::Abs(Ap), TMath::Abs(Am)) << " & ";

		//Vertex eff systematic of 1%
		
		EEp+=TMath::Power(1,2); 
		EEm+=TMath::Power(1,2); 
		
		output_table << "1 "<< " & ";

		
		//SG eff systematic 

		double pt=((TH1D*)output_files->FindObjectAny(yield_histogram(iy,ups,"").c_str()))->GetBinCenter(ipt);
		double SGeffP=SG_eff->GetParameter(0)+SG_eff->GetParError(0) + pt*(SG_eff->GetParameter(1)+SG_eff->GetParError(1));
		double SGeffM=SG_eff->GetParameter(0)-SG_eff->GetParError(0) + pt*(SG_eff->GetParameter(1)-SG_eff->GetParError(1));

		//cout << "SG Eff: " << SG_eff->Eval(pt) << " SGeffp: " << SGeffP << endl; 
		
		//double SGp=100*TMath::Abs(SG_eff->Eval(pt)-SGeffP);
		//double SGm=100*TMath::Abs(SG_eff->Eval(pt)-SGeffM);

		double SGp=2;
		double SGm=2;
		
		output_table << Form("%.1f(%.1f)",SGp, SGm) << " & ";

		//cout << "Pt: " << pt << " BP: " << BP << " SP: " << SP << " EffP: " << EffP << " Ap: " << Ap  << " SGp: " << SGp << endl;
		
		EEp+=TMath::Power(SGp,2); 
		EEm+=TMath::Power(SGm,2); 
	//	cout << "Sg: " << EEp << endl; 

		// sqrt of sum 
		
		double Ep=TMath::Sqrt(EEp); 
		double Em=TMath::Sqrt(EEm); 
		
		output_table << Form("%.1f(%.1f) & ",Ep, Em);
		output_table << Form("%.1f",yE)<< " \\\\" << endl;

		//cout << "pt: " << pt << " Total Systematic: " << Ep << endl; 
		
		total_sysM_per->SetBinContent(ipt,Em); 
		total_sysP_per->SetBinContent(ipt,Ep); 
		
		total_sysP_per->SetBinError(ipt,0); 
		total_sysM_per->SetBinError(ipt,0); 
		//cout << "Shape: " << SP << " BKG: "<< BP << " Eff: " << EffP << " Acc: " << Ap << " Tot: " << Ep << endl;
	}
	output_table.close();
	output_files->cd();

	C->cd(1); 
	EffP_per->SetLineColor(kRed); 
	EffM_per->SetLineColor(kRed);
	EffM_per->SetLineStyle(kDashed); 
	
	Fitsys->SetLineColor(kBlue); 
		
	L->AddEntry(EffP_per,"Efficiency Sys+");
	L->AddEntry(EffM_per,"Efficiency Sys-"); 

	L->AddEntry(Fitsys,"Fit systematic error"); 
	
	double SM[]={EffP_per->GetBinContent(EffP_per->GetMaximumBin()),EffM_per->GetBinContent(EffM_per->GetMaximumBin())};
	double max = TMath::MaxElement(4,SM);//max element of N elements in an array 
	
	EffP_per->SetMaximum(max*1.1);
	EffP_per->SetMinimum(0); 
	EffP_per->GetYaxis()->SetTitle("Error [%]");
	EffP_per->GetXaxis()->SetTitle("p_{T}(#mu#mu)"); 
	EffP_per->SetTitle(Form("#Upsilon(%dS)",ups)); 
	EffP_per->SetStats(kFALSE);
	EffP_per->Draw(); 
	EffM_per->Draw("same"); 
	Fitsys->Draw("same"); 
	
	C->cd(2); 
	L->Draw(); 
	
	total_sysP_per->Write(); 
	total_sysM_per->Write();
	EffM_per->Write();
	EffP_per->Write();
	Fitsys->Write(); 
	rhoP_per->Write(); 
	C->Write(); 
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/Fit_systematics_summary_y%.1f_%.1f_%dS.pdf",fYbin2[iy],fYbin2[iy+1],ups));
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/Fit_systematics_summary_y%.1f_%.1f_%dS.jpg",fYbin2[iy],fYbin2[iy+1],ups));

	
	delete C;
	delete EffP_per; 
	delete EffM_per; 
	delete Fitsys; 
	delete L; 
	
	delete total_sysP_per; 
	delete total_sysM_per; 
	
	delete output_filesEp; 
	delete output_filesEm; 
	
	delete SG_file;
	delete SG_eff; 
}

void yield_GeV(int iy, string mode){
	TH1D *yield1S=(TH1D*)output_files->FindObjectAny(yield_histogram(iy,1,mode.c_str()).c_str());
	TH1D *yield2S=(TH1D*)output_files->FindObjectAny(yield_histogram(iy,2,mode.c_str()).c_str());
	TH1D *yield3S=(TH1D*)output_files->FindObjectAny(yield_histogram(iy,3,mode.c_str()).c_str());

	
	TH1D *yield1S_GeV=(TH1D*)yield1S->Clone(yield_histogram_GeV(iy,1,mode.c_str()).c_str());
	TH1D *yield2S_GeV=(TH1D*)yield2S->Clone(yield_histogram_GeV(iy,2,mode.c_str()).c_str());
	TH1D *yield3S_GeV=(TH1D*)yield3S->Clone(yield_histogram_GeV(iy,3,mode.c_str()).c_str());
	
	yield1S_GeV->Sumw2();
	yield2S_GeV->Sumw2();
	yield3S_GeV->Sumw2();

	yield1S_GeV->Divide(bin_width); 
	yield2S_GeV->Divide(bin_width); 
	yield3S_GeV->Divide(bin_width); 


	output_files->cd();
	yield1S_GeV->Write();
	yield2S_GeV->Write();
	yield3S_GeV->Write(); 
	
	delete yield1S;
	delete yield2S;
	delete yield3S;
	
	delete yield1S_GeV;
	delete yield2S_GeV;
	delete yield3S_GeV; 
	
}

void write_histos(){
	output_files->cd("fit_parameters"); 
	TCanvas *C = new TCanvas("deltaM_pt_canvas","",600,600);
	deltaM_pt->SetStats(kFALSE); 
	deltaM_pt->Draw();
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	
	C->Write();
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/deltaM_fit.pdf"));
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/deltaM_fit.jpg"));

	
	dm_scale_pt->SetStats(kFALSE);
	dm_scale_pt->Draw(); 
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.Draw();
	lumiText.Draw();
	C->SetName("dm_scale_pt_canvas"); 
	C->Write();
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/cw_fit.pdf"));
	C->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/cw_fit.jpg"));

	
	deltaM_pt->Write();
	dm_scale_pt->Write();
	output_files->cd();
	bin_width->Write();
	
	output_files->cd("fit_parameters"); 
	res_pt_P1S->Write();
	res_pt_M1S->Write();
	
	S1_pt->Write();
	S2_pt->Write();
	S3_pt->Write();

	
	chi2_pt->Write();
	
	Y1_bkg_4->Write();
	Y1_bkg_3->Write();
	Y1_bkg_2->Write();
	Y1_bkg_1->Write();
	eff->Write();
	
	CName["fits"]->SetFillColor(10);
	CName["fits2"]->SetFillColor(10); 
	CName["fits"]->Write(); 
	CName["fits2"]->Write();
	CName["fits"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits_summary.pdf"); 
	CName["fits"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits_summary.jpg"); 

	CName["fits2"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits_summary2.pdf"); 
	CName["fits2"]->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits_summary2.jpg"); 

	/*
	r12->Draw("E1");
	C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/r12.pdf");
	r13->Draw("E1"); 
	C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/r13.pdf");

	r12->Write();
	r13->Write();
	 */
	delete C;
	output_files->cd();

}

void draw_background_errors(){

	TCanvas *CArea=new TCanvas("Background_Function");
	CArea->Divide(2,2); 
	CArea->cd(1);
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_4"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_4"))->SetMarkerColor(kBlue); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_4"))->SetMinimum(0.3);
	((TH1D*)output_files->FindObjectAny("Y1_bkg_4"))->Draw("E1"); 

	((TH1D*)output_files->FindObjectAny("Y1_bkg_3"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_3"))->SetMarkerColor(kRed); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_3"))->Draw("E1 same"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_2"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_2"))->SetMarkerColor(kGreen); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_2"))->Draw("E1 same"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_1"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_1"))->SetMarkerColor(kBlack); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_1"))->Draw("E1 same"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_0"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_0"))->SetMarkerColor(kOrange); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_0"))->Draw("E1 same"); 
	
	TLegend *L = new TLegend(0.25,0.5,0.9,0.9); 
	L->AddEntry(((TH1D*)output_files->FindObjectAny("Y1_bkg_4")));
	L->AddEntry(((TH1D*)output_files->FindObjectAny("Y1_bkg_3")));
	L->AddEntry(((TH1D*)output_files->FindObjectAny("Y1_bkg_2")));
	L->AddEntry(((TH1D*)output_files->FindObjectAny("Y1_bkg_1")));
	L->AddEntry(((TH1D*)output_files->FindObjectAny("Y1_bkg_0")));

	CArea->cd(2); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_4"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_4"))->SetMarkerColor(kBlue); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_4"))->SetMaximum(10);
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_4"))->Draw("histo P"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_3"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_3"))->SetMarkerColor(kRed); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_3"))->Draw("histo P same"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_2"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_2"))->SetMarkerColor(kGreen); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_2"))->Draw("histo P same"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_1"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_1"))->SetMarkerColor(kBlack); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_1"))->Draw("histo P same"); 
	
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_0"))->SetMarkerStyle(7); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_0"))->SetMarkerColor(kOrange); 
	((TH1D*)output_files->FindObjectAny("Y1_bkg_chi2_0"))->Draw("histo P same"); 
	
	
	
	CArea->cd(3);
	L->Draw(); 
	
	CArea->Write();
}

int get_mem(long* vmrss_kb, long* vmsize_kb){
	 //Get the the current process' status file from the proc filesystem
	 //Enter the job PID to get the current process
	 //This code is adapted from http://www.umbc.edu/hpcf/resources-tara/checking-memory-usage.html
	 //The output is VmSize and VmRSS. VmSize is the total amount of memory required by the program, while VmRSS is the "Resident Set Size",
	 //Or the amount of memory use right now. (https://wiki.duke.edu/display/SCSC/Monitoring+Memory+Usage )
	 
	 
	int PID=gSystem->GetPid();
	FILE* procfile = fopen(Form("/proc/%d/status",PID), "r");
	
	const long to_read = 8192;
	char buffer[to_read];
	int read = fread(buffer, sizeof(char), to_read, procfile);
	
	fclose(procfile);
	
	short found_vmrss = 0;
	short found_vmsize = 0;
	char* search_result;
	
	// Look through proc status contents line by line 
	char delims[] = "\n";
	char* line = strtok(buffer, delims);
	
	while (line != NULL && (found_vmrss == 0 || found_vmsize == 0) )
    {
		search_result = strstr(line, "VmRSS:");
		if (search_result != NULL)
        {
			sscanf(line, "%*s %ld", vmrss_kb);
			found_vmrss = 1;
		}
		
		search_result = strstr(line, "VmSize:");
		if (search_result != NULL)
        {
			sscanf(line, "%*s %ld", vmsize_kb);
			found_vmsize = 1;
        }
		
		line = strtok(NULL, delims);
    }
	
	return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}

double BKG_cheb(double *x, double *par){
	//There are always 6 parameters
	// par[0]=A0
	// par[1]=A1
	// par[2]=A2
	// par[3]=A3
	
	//Returns function A0*T0+A1*T1+A2*T3+A4*T4 over the mass interval 8.7-11.2
	
	//par[4] and par[5] used for exclusions. If you don't want to exclude any of the mass region, set to <0. 
	// if you want to use the exclusion, it will exclude the mass region from M1-M2, which is assumed to be at
	// the center 
	double m=x[0];
	
	double M1=par[5];	
	double M2=par[6]; 
	
	double m0=(M2-M1)/2+M1;//center of mass range
	double m_=2*(m-m0)/(M2-M1);//redefine interval
	
	//Chebyshevs up to fourth order 
	double T0=1; 
	double T1=m_;
	double T2=2*m_*m_-1;
	double T3=4*m_*m_*m_-3*m_;
	double T4=8*TMath::Power(m_,4)-8*m_*m_+1;
	
	double BKG;
	BKG=par[0]*T0 + par[1]*T1 + par[2]*T2 + par[3]*T3 + par[4]*T4; // pol 2
	
	return BKG;
}

double sig_bkg2(double *x, double *par){

	//Fit function for 3 smeared signals (histograms) + background function
	//par[4]=0 is the signal/bkg fraction
	//parameter 5 is the background yield 
	// parameter 6 etc are for the background function  description 
	bool print=0; 
	double m=x[0]; 
	
	double S1=0; 
	double S2=0; 
	double S3=0; 
	
	double dm_scale=par[3];
	double deltaM=par[4];
	
	if(print) cout << "Before eval. " << endl;
	S1=LS1_shape->Eval(m-deltaM, dm_scale); 
	S2=LS2_shape->Eval(m-par[12], par[13]); 
	S3=LS3_shape->Eval(m-par[12], par[13]); 
	
	S1=par[0]*S1;
	S2=par[0]*par[1]*S2; 
	S3=par[0]*par[2]*S3; 
	//S2=par[1]*S2; 
	//S3=par[2]*S3;  
	
//	double F=par[5]; 
//	double ys=par[5]; // normalization factor 
	if(print) cout << "par[0]: " << par[0] << " par[1]: " << par[1] << " par[2]: " << par[2] << endl; 
	return (S1+S2+S3)+BKG_cheb(x,&par[5]);
	
}

double sig_bkg(double *x, double *par){
	
	//Fit function for 3 smeared signals (histograms) + background function
	//par[4]=0 is the signal/bkg fraction
	//parameter 5 is the background yield 
	// parameter 6 etc are for the background function  description 
	bool print=0; 
	double m=x[0]; 
	
	double S1=0; 
	double S2=0; 
	double S3=0; 
	
	double dm_scale=par[3];
	double deltaM=par[4];
	
	if(print) cout << "Before eval. " << endl;
	S1=LS1_shape->Eval(m-deltaM, dm_scale); 
	S2=LS2_shape->Eval(m-deltaM, dm_scale); 
	S3=LS3_shape->Eval(m-deltaM, dm_scale); 
	
	S1=par[0]*S1;
	S2=par[0]*par[1]*S2; 
	S3=par[0]*par[2]*S3; 
	//S2=par[1]*S2; 
	//S3=par[2]*S3;  
	
	//	double F=par[5]; 
	//	double ys=par[5]; // normalization factor 
	if(print) cout << "par[0]: " << par[0] << " par[1]: " << par[1] << " par[2]: " << par[2] << endl; 
	return (S1+S2+S3)+BKG_cheb(x,&par[5]);
	
}

double Cheb_int(double a,double b, int order){
	double Int=0; 
	TF1 *T=new TF1("T",BKG_cheb,8.7,11.2,7); 
	T->FixParameter(5,8.7); 
	T->FixParameter(6,11.2); 
	T->FixParameter(4,0);
	T->FixParameter(3,0); 
	T->FixParameter(2,0); 
	T->FixParameter(1,0); 
	T->FixParameter(0,0);
	
	T->FixParameter(order,1); 
	Int=T->Integral(a,b); 
	
	delete T; 
	
	return Int; 
	
}

double integral_error(TF1 *FF, double a, double b){
	
	double Err=0; 
	
	for (int i=0; i<4; i++) {
		double Int=Cheb_int(a,b,i);
		Err+=TMath::Power(Int*FF->GetParError(i)*FF->GetParameter(i),2);
	}
	
	return TMath::Sqrt(Err);
	
}

double signal_shape(double *x, double *par){
	double m=x[0]; 
	
	int ups = static_cast<int>(par[0]);
	
	double dm_scale=par[2];
	double deltaM=par[3];
	
	double S=0; 
	
	if(ups==1) S=LS1_shape->Eval(m-deltaM, dm_scale); 
	if(ups==2) S=LS2_shape->Eval(m-deltaM, dm_scale); 
	if(ups==3) S=LS3_shape->Eval(m-deltaM, dm_scale); 
	
	S=par[1]*S; 
	
	return S; 
}

double fit_bkg_N(TH1D *mdata,TF1 *BKG_fixed, int ipt){
	TH1D *dataN=(TH1D*)mdata->Clone();
	
	double M1=9.0;
	double M2=10.6;
	
	double N1=dataN->Integral(dataN->FindBin(8.7),dataN->FindBin(M1));
	double N2=dataN->Integral(dataN->FindBin(M2),dataN->FindBin(11.2));
	
	
	BKG_fixed->SetParName(0,"A0"); 
	BKG_fixed->SetParName(1,"A1"); 
	BKG_fixed->SetParName(2,"A2"); 
	BKG_fixed->SetParName(3,"A3"); 
	BKG_fixed->SetParName(4,"A4"); 
	
	BKG_fixed->SetParName(5,Form("M1: m<%.1f",M1)); 
	BKG_fixed->SetParName(6,Form("M2: m>%.1f",M2)); 
	
	BKG_fixed->FixParameter(5,M1);//Select mass region
	BKG_fixed->FixParameter(6,M2);
		
	string fit_method="QLIS0";
	if(fPTbin2[ipt]>LL_thres) fit_method="QSLS0"; 
	else fit_method="QIS0"; 
	
	TFitResultPtr RR=dataN->Fit(BKG_fixed,fit_method.c_str()); 
	
	double chi2= RR->Chi2();
	double NDOF= static_cast<double> (RR->Ndf());
	
	BKG_fixed->FixParameter(5,8.7); 
	BKG_fixed->FixParameter(6,11.2); 
	delete dataN;
	double ys=(BKG_fixed->Integral(8.7,M1)+BKG_fixed->Integral(M2,11.2))/BKG_fixed->Integral(8.7,11.2);
	return ys; 
}

void estimate_yield(double *y){
	
	y[0]=mdata->Integral(y1m[0],y1m[1]);
	y[1]=mdata->Integral(y2m[0],y2m[1]);
	y[2]=mdata->Integral(y3m[0],y3m[1]);

}

double PDF_shape(double *x, double *par){
	
	double f=0; 
	double neventd=static_cast<double> (N_event);
	
	for(int i=0; i<N_event; i++){
		int im=2*i;//index for mean
		int is=2*i+1;//index for sigma
		double meanprime=par[im]; // +x[1]
		double sigmaprime=par[is]*x[1];
		double norm=1./( TMath::Sqrt(2*TMath::Pi())*sigmaprime); // Normalization
		f+=norm*TMath::Gaus(x[0],meanprime,sigmaprime);// Gaussian at point x, with position mean, with width sigma, normalized by 1/2pi
		
	}
	
	return f/neventd; 
}

double get_bin_center(int ipt,int ups){
	//cout << "Get Bin Center: " << endl; 
	
	TH1D *pt_dist=(TH1D*)data->FindObjectAny(Form("Y%dPt",ups)); 
	
	TH1D *pt_tmp=(TH1D*)pt_dist->Clone("tmp");
	
	pt_tmp->SetAxisRange(fPTbin2[ipt],fPTbin2[ipt+1]-0.01); 
	double BC=pt_tmp->GetMean(); 
	
	//cout << "Bin Center: " << BC << endl; 
	
	delete pt_dist;
	delete pt_tmp;
	
	if(ipt==21) BC=82; 
	
	return BC; 
}

void check_background_model(int iy, int ipt){
	TF1 *BKG_fixed=new TF1("BKG_fixed",BKG_cheb,8.7,11.2,7); 
	
	double bin_center1S=get_bin_center(ipt,1); //(fPTbin2[ipt+1]+fPTbin2[ipt])/2;
	BKG_fixed->FixParameter(5,8.7);
	BKG_fixed->FixParameter(6,11.2);
	for (int i=5; i>0; i--) {
		if(i<5)BKG_fixed->FixParameter(i,0);
		fit_bkg_N(mdata,BKG_fixed,ipt); 
		double Int=BKG_fixed->Integral(y1m[0],y1m[1])/mdata->GetBinWidth(1); 
		double IntE=integral_error(BKG_fixed,y1m[0],y1m[1]);
		
		BKG_fixed->FixParameter(5,9.0);
		BKG_fixed->FixParameter(6,10.6);
		
		int NDF=BKG_fixed->GetNDF();
		
		double chi2=compute_resid(iy,ipt,mdata,BKG_fixed, NDF); 
		double chi2_NDOF=chi2/static_cast<double>(NDF);//BKG_fixed->GetChisquare()/static_cast<double>(BKG_fixed->GetNDF());
	
		BKG_fixed->FixParameter(5,8.7);
		BKG_fixed->FixParameter(6,11.2);
		
	
		//cout << "order: " << i-1 << " integral: " << Int<< " +/- " << IntE <<endl; 
		//cout << "order: " << i-1 << " chi2/NDOF: "<< chi2_NDOF << endl << endl; 
		if(i==5){
			Y1_bkg_4->Fill(bin_center1S,Int); 
			Y1_bkg_4->SetBinError(Y1_bkg_4->FindBin(bin_center1S),IntE); 
			Y1_bkg_chi2_4->Fill(bin_center1S,chi2_NDOF);
		}
		if(i==4){
			Y1_bkg_3->Fill(bin_center1S,Int); 
			Y1_bkg_3->SetBinError(Y1_bkg_4->FindBin(bin_center1S),IntE); 
			Y1_bkg_chi2_3->Fill(bin_center1S,chi2_NDOF);
		}
		if(i==3){
			Y1_bkg_2->Fill(bin_center1S,Int); 
			Y1_bkg_2->SetBinError(Y1_bkg_4->FindBin(bin_center1S),IntE); 
			Y1_bkg_chi2_2->Fill(bin_center1S,chi2_NDOF);
		}
		if(i==2){
			Y1_bkg_1->Fill(bin_center1S,Int); 
			Y1_bkg_1->SetBinError(Y1_bkg_4->FindBin(bin_center1S),IntE); 
			Y1_bkg_chi2_1->Fill(bin_center1S,chi2_NDOF);
		}
		
		if(i==1){
			Y1_bkg_0->Fill(bin_center1S,Int); 
			Y1_bkg_0->SetBinError(Y1_bkg_4->FindBin(bin_center1S),IntE); 
			Y1_bkg_chi2_0->Fill(bin_center1S,chi2_NDOF);
			//cout << "par: " << BKG_fixed->GetParameter(0) << " +/- " << BKG_fixed->GetParError(0) << endl;
		}
		
		
	}
	
	delete BKG_fixed;
}

void set_bkg_order(int n,TF1 *LS_bkg){
	
	
	LS_bkg->ReleaseParameter(7);//n=2
	LS_bkg->ReleaseParameter(8);//n=3
	LS_bkg->ReleaseParameter(9);//n=4
	
	
	if(n==1) {
		LS_bkg->FixParameter(7,0);
		LS_bkg->FixParameter(8,0);
		LS_bkg->FixParameter(9,0); 
	}
	
	if(n==2){
		LS_bkg->FixParameter(8,0);
		LS_bkg->FixParameter(9,0); 
	}
	if(n==3){
		LS_bkg->FixParameter(9,0); 
	}
	
	if(n==4) return; 
	
}

void fit_mass_range(TF1 *FF, int ipt, double *Y, double *YE,double *R, double *RE, double &chi2, int &NDOF, string fit_method){
	//cout << "Fit Mass Range. " << endl; 
	TF1 *LS_bkg;
	int Npar=FF->GetNpar();
	if(Npar==13) LS_bkg=new TF1("LS_bkg", sig_bkg2, 8.5,11.5,Npar);//13 
	else LS_bkg=new TF1("LS_bkg", sig_bkg, 8.5,11.5,Npar); 
	
	LS_bkg->SetNpx(10000); 
	
	double BW=mdata->GetBinWidth(1);
	
	int order=-1; 
	
	for (int i=0; i<Npar; i++) {
		double par=FF->GetParameter(i);
		double parE=FF->GetParError(i);
		LS_bkg->ReleaseParameter(i);
		if(parE==0)LS_bkg->FixParameter(i,par); 
		else LS_bkg->SetParameter(i,par); 
		LS_bkg->SetParName(i,FF->GetParName(i)); 
		TString parName=FF->GetParName(i); 
		TString keyword1("A"); 
		if(parName.Contains(keyword1) && parE==0 && order<0){
			LS_bkg->ReleaseParameter(i);
			order=i; 
		}
	}
	
	//n_opt[ipt]=4;
	set_bkg_order(n_opt[ipt],LS_bkg); 
	
	bool good_fit=false; 
	
	TRandom3 ML(0); 
	TRandom3 MH(0); 
	
	int counter=0; 
	
	while (good_fit==false){
		//cout << "Fitting. " << endl; 
		double mrL=ML.Uniform(fit_rangesL[0],fit_rangesL[1]); 
		double mrH=MH.Uniform(fit_rangesH[0],fit_rangesH[1]); 
		counter++; 
		if(counter>10) break;
		LS_bkg->SetRange(mrL, mrH); 
		LS_bkg->FixParameter(10,mrL);
		LS_bkg->FixParameter(11,mrH); 
		LS_bkg->Update();
		double temp[]={LS_bkg->GetParameter(0),LS_bkg->GetParameter(1)*LS_bkg->GetParameter(0),LS_bkg->GetParameter(2)*LS_bkg->GetParameter(0)};
		TFitResultPtr RR=mdata->Fit(LS_bkg, fit_method.c_str());

		TMatrixDSym cor = RR->GetCorrelationMatrix();
		//cout << "cov[0][0] " << cor[0][0] << endl;
		//	RR->Print("V");
		
	//	cout << "ups error: " << ups_errors(LS_bkg,cor,2) << endl; 
	//	cout << "ups error: " << ups_errors(LS_bkg,cor,3) << endl;
		
		if(temp[0]==LS_bkg->GetParameter(0) || temp[1]==LS_bkg->GetParameter(1) || temp[2]==LS_bkg->GetParameter(2)){
			cout << "Same Parameters. " << endl; 
			cout << "Y1: " << LS_bkg->GetParameter(0)/BW << endl; 
			continue; 
		}
		
		if(LS_bkg->GetParError(1)/LS_bkg->GetParameter(1)>1 || LS_bkg->GetParError(2)/LS_bkg->GetParameter(2)>1){
			continue;
		}
		
		TString Report=gMinuit->fCstatu; 
		if(Report.Contains("RESET")==1) cout <<"Reset" << endl;  
		if(Report.Contains("CONVERGED")==1 && LS_bkg->GetParameter(0)>0 && LS_bkg->GetParameter(1)>0 && LS_bkg->GetParameter(2)>0){
			if(gMinuit->fStatus!=0)cout << "Fit Status: " << gMinuit->fStatus << endl; 
			good_fit=true; 
			//cout <<"Fit Converged. " << endl; 
			//if(ipt>=15 && ipt < 18) cout << LS_bkg->GetParameter(0)/BW << endl; 
			double Y1=LS_bkg->GetParameter(0)/BW; 
			double Y2=LS_bkg->GetParameter(0)*LS_bkg->GetParameter(1)/BW; 
			double Y3=LS_bkg->GetParameter(0)*LS_bkg->GetParameter(2)/BW;
			
			R[0]=LS_bkg->GetParameter(1);
			R[1]=LS_bkg->GetParameter(2);
			RE[0]=LS_bkg->GetParError(1);
			RE[1]=LS_bkg->GetParError(2); 
			
			Y[0]=Y1;
			Y[1]=Y2;
			Y[2]=Y3;
			
			double S1E=(LS_bkg->GetParError(0)/BW);
			double S2E=(ups_errors(LS_bkg,cor,2)/BW);
			double S3E=(ups_errors(LS_bkg,cor,3)/BW);
			
			YE[0]=S1E; 
			YE[1]=S2E;
			YE[2]=S3E;
			double C=0; 
			int N=0; 
			compute_chi2(mdata,LS_bkg, 9.1,10.6,C,N); 
			chi2=C;
			NDOF=N; 
			//chi2=LS_bkg->GetChisquare()/static_cast<double>(LS_bkg->GetNDF() ) ; 
			
		}//Converged 
	}//while false
	
	delete LS_bkg; 
	
}

void compute_mass_sys(int iy, int ipt, double &chi2_min, int &NDOF_min){
	const int N=1000; 
	cout << "Computing Mass Systematic for pT: " << fPTbin2[ipt] << "-" << fPTbin2[ipt+1] << endl; 
	
	TF1 *FF=(TF1*)output_files->FindObjectAny(fit_function(iy,ipt,"").c_str());
	//TF1 *FFnom=(TF1*)nom->FindObjectAny(fit_function(iy,ipt,"").c_str());
	
	string fit_method; 
	if(fPTbin2[ipt]>LL_thres) fit_method="QLRS0"; 
	else if(fPTbin2[ipt]<70) fit_method="QRS0";
	else fit_method="QRS0I";


	double BW=mdata->GetBinWidth(1);
	double YL=(FF->GetParameter(0)-FF->GetParameter(0)*0.5)/BW; 
	double YH=(FF->GetParameter(0)+FF->GetParameter(0)*0.5)/BW;

   // cout << "Yield Range (1S): " << YL << "-" << YH << endl; 
	
	
	TH1D *Ychi2=new TH1D(Form("chi2_y%d_pt%d",iy,ipt),"#chi^{2}",30,0,3); 
	TH1D *failing_mass=new TH1D(Form("failing_mass_y%d_pt%d",iy,ipt),"Failing Mass Window", 1250,8.5,11.5); 
	TH1D *failing_deltaM=new TH1D(Form("failing_deltaM_y%d_pt%d",iy,ipt),"Failing Mass Window",500,1.5,3.); 
	TH1D *passing_deltaM=new TH1D(Form("passing_deltaM_y%d_pt%d",iy,ipt),"Passing Mass Window",500,1.5,3.0); 
	

	TStopwatch t; 
	t.Start(); 
	int i=0; 
	
	double Y[]={-9,-9,-9};
	double YE[3]; 
	double R[]={-9,-9};
	double RE[]={-9,-9};
	double chi2=0;
	int ND=0; 
	
	double YR1[N];
	double YR2[N];
	double YR3[N];
	
	double YR1E[N];
	double YR2E[N];
	double YR3E[N];
	
	double R21[N];
	double R31[N];
	double R21E[N];
	double R31E[N];
	
	double chi2_array[N];

	while (i<N) {
		if(i%100 ==0) cout << "Event: " << i << endl; 
		//if(i>0) fit_method=fit_method+"Q";
		fit_mass_range(FF, ipt,Y,YE,R,RE,chi2,ND,fit_method);
		if(Y[0]>0 && Y[1]>0 && Y[2]>0){
			//cout << Y[0] << " " << Y[1] << " "<< Y[2] << endl; 
			YR1[i]=Y[0];
			YR2[i]=Y[1];
			YR3[i]=Y[2];
			
			YR1E[i]=YE[0];
			YR2E[i]=YE[1];
			YR3E[i]=YE[2];
			
			R21[i]=R[0];
			R31[i]=R[1];

			R21E[i]=RE[0];
			R31E[i]=RE[1];

			chi2_array[i]=chi2;
			NDOF_min=ND; 
			
			i++;
		}
	}//while loop 
	
	YL=TMath::MinElement(N,YR1); 
	YH=TMath::MaxElement(N,YR1); 

	chi2_min=TMath::MinElement(N,chi2_array); 
	
	double YLE=TMath::MinElement(N,YR1E); 
	double YHE=TMath::MaxElement(N,YR1E);
	
	YL-=YL*0.001;
	YH+=YH*0.001;
	YLE-=YLE*0.01;
	YHE+=YHE*0.01;

	TH2D *yield_E1= new TH2D(Form("yield1S_E_y%d_pt%d",iy,ipt),"YieldE vs Yield; Yield(1S); #sigma_{Yield}(1S)", 100, YL, YH, 100, YLE,YHE); 
	
	double RL=TMath::MinElement(N,R21); 
	double RH=TMath::MaxElement(N,R21);
	double RLE=TMath::MinElement(N,R21E);
	double RHE=TMath::MaxElement(N,R21E); 
	
	TH2D *R21_E= new TH2D(Form("r21_E_y%d_pt%d",iy,ipt),"Ratio vs. Ratio E; r21; #sigma_{r21}", 100, RL, RH, 100, RLE,RHE); 
	
	
	YL=TMath::MinElement(N,YR2); 
	YH=TMath::MaxElement(N,YR2); 
	YLE=TMath::MinElement(N,YR2E); 
	YHE=TMath::MaxElement(N,YR2E);
	YL-=YL*0.001;
	YH+=YH*0.001;
	YLE-=YLE*0.01;
	YHE+=YHE*0.01;
	
	TH2D *yield_E2= new TH2D(Form("yield2S_E_y%d_pt%d",iy,ipt),"YieldE vs Yield; Yield(2S); #sigma_{Yield}(2S)", 100, YL, YH, 100, YLE,YHE); 
	
	YL=TMath::MinElement(N,YR3); 
	YH=TMath::MaxElement(N,YR3); 
	YLE=TMath::MinElement(N,YR3E); 
	YHE=TMath::MaxElement(N,YR3E);
	YL-=YL*0.001;
	YH+=YH*0.001;
	YLE-=YLE*0.01;
	YHE+=YHE*0.01;
	
	TH2D *yield_E3= new TH2D(Form("yield3S_E_y%d_pt%d",iy,ipt),"YieldE vs Yield; Yield(3S); #sigma_{Yield}(3S)", 100, YL, YH, 100, YLE,YHE); 
		
	 RL=TMath::MinElement(N,R31); 
	 RH=TMath::MaxElement(N,R31);
	 RLE=TMath::MinElement(N,R31E);
	 RHE=TMath::MaxElement(N,R31E); 
	
	TH2D *R31_E= new TH2D(Form("r31_E_y%d_pt%d",iy,ipt),"Ratio vs. Ratio E; r31; #sigma_{r31}", 100, RL, RH, 100, RLE,RHE); 
	
	for(int i=0; i<N; i++){
		yield_E1->Fill(YR1[i],YR1E[i]);
		yield_E2->Fill(YR2[i],YR2E[i]);
		yield_E3->Fill(YR3[i],YR3E[i]);
		R21_E->Fill(R21[i],R21E[i]);
		R31_E->Fill(R31[i],R31E[i]);
		}
	
	
	t.Stop(); 
	cout << "Real Time: " << t.RealTime() << endl; 
	
	output_files->cd("range_sys"); 

	yield_E1->Write();
	yield_E2->Write();
	yield_E3->Write();

	R31_E->Write();
	R21_E->Write(); 
	
	Ychi2->Write();
	failing_mass->Write();
	failing_deltaM->Write(); 
	passing_deltaM->Write(); 
	
	delete Ychi2; 
	
	delete yield_E1; 
	delete yield_E2;
	delete yield_E3;
	
	delete R21_E; 
	delete R31_E;
	
	delete failing_deltaM;
	delete failing_mass;
	delete passing_deltaM; 
	delete FF;
	
}

void fill_ratio_yield(int iy, int num){
	if(num<2 || num >3 ) return; 
	//only ompute the ratio for 2/1 and 3/1
	TH1D *Rn1pT = new TH1D(ratio_histogram(iy,num,"").c_str(),Form("R%d1",num),fNpt2, fPTbin2); 
	TH1D *Rn1EpT = new TH1D(ratio_histogram(iy,num,"fit_sys").c_str(),Form("R%d1 systematic",num),fNpt2, fPTbin2); 
	
	for (int ipt=1; ipt<=Rn1pT->GetNbinsX(); ipt++) {
		double pt=Rn1pT->GetBinCenter(ipt); 
		TH2D *Rn1_E=(TH2D*)output_files->FindObjectAny(Form("r%d1_E_y%d_pt%d",num,iy,ipt-1)); 
		if(ipt==1 && iy==0){
			TCanvas * tmp = new TCanvas("","",800,800); 
			Rn1_E->SetTitle(BT(iy,ipt-1).c_str()); 
			Rn1_E->SetStats(kFALSE);
			Rn1_E->Draw("COLZ"); 
			tmp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/R%d1_E_y.pdf",num));
			tmp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/R%d1_E_y.jpg",num)); 

			delete tmp; 
		}
		TH1D *sigma_Rn1=(TH1D*)Rn1_E->ProjectionY(); 
		if(ipt==1 && iy==0){
			TCanvas * tmp = new TCanvas("","",800,800); 
			sigma_Rn1->SetTitle(BT(iy,ipt-1).c_str()); 
			sigma_Rn1->SetStats(kFALSE); 
			sigma_Rn1->Draw(); 
			tmp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/R%d1E_dist.pdf",num)); 
			tmp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/R%d1E_dist.jpg",num)); 

			delete tmp; 
		}
		
		double mean_sigma=sigma_Rn1->GetMean();//mean fit error sigma_Rn1->GetBinCenter(sigma_Rn1->GetMaximumBin());
		double RMS_sigma=sigma_Rn1->GetRMS(); //spread in fit errors 
		
		TH1D *R=(TH1D*)Rn1_E->ProjectionX("R",sigma_Rn1->FindBin(mean_sigma-RMS_sigma),sigma_Rn1->FindBin(mean_sigma+RMS_sigma)); 
		
		if(ipt==1 && iy==0){
			TCanvas * tmp = new TCanvas("","",800,800); 
			R->SetTitle(BT(iy,ipt-1).c_str()); 
			R->SetStats(kFALSE); 
			TLatex LT(0.2,0.8,Form("%.3f<#sigma_{R%d1}<%.3f",mean_sigma-RMS_sigma,num,mean_sigma+RMS_sigma));
			LT.SetNDC(kTRUE);
			R->Draw(); 
			LT.Draw();
			tmp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/R%d1_dist.pdf",num)); 
			tmp->Print(Form("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/R%d1_dist.jpg",num)); 

			delete tmp; 
		}
		
		double Ratio=R->GetMean();//R->GetBinCenter(R->GetMaximumBin())
		double R_RMS=R->GetRMS();
		double R_stat=mean_sigma;
		double Rsys=TMath::Sqrt(R_RMS*R_RMS+RMS_sigma*RMS_sigma);
		
		Rn1pT->SetBinContent(ipt,Ratio); //nominal ratio 
		Rn1pT->SetBinError(ipt,R_stat);//statistical error  
		Rn1EpT->SetBinContent(ipt,Rsys*100/Ratio); //systematic error 
		Rn1EpT->SetBinError(ipt,0); // no bin error, this keeps track of fit systematics only 
		
		delete Rn1_E; 
		delete sigma_Rn1;
		delete R;
	}//end pT loop 
	output_files->cd("yield_histograms");
	Rn1pT->Write(); 
	Rn1EpT->Write();
	
	delete Rn1pT;
	delete Rn1EpT; 
	
	
}

void fill_toy_yield(int ups,int iy){
//	yield_histogram(iy,1,mode.c_str()).c_str()
//	Form("toy_yield_%dS_y%d",ups,iy)
	
	TH1D *toy_yield = new TH1D(yield_histogram(iy,ups,"").c_str(), Form("Yield(%dS)",ups), fNpt2,fPTbin2);
	TH1D *yield_actual = new TH1D(yield_histogram(iy,ups,"actual").c_str(), Form("Actual Yield(%dS)",ups), fNpt2,fPTbin2);

	TH1D *toy_RMS = new TH1D(Form("toy_RMSfrac_%dS_y%d",ups,iy),Form("RMS/Yield(%dS);p_{T}(#mu#mu); RMS/y ['%']",ups), fNpt2,fPTbin2);
	TH1D *toy_RMS_stat = new TH1D(Form("toy_RMS_stat_%dS_y%d",ups,iy),Form("RMS/stat %dS; p_{T}(#mu#mu); RMS/stat",ups),fNpt2, fPTbin2); 
	TH1D *toy_yield_sys = new TH1D(yield_histogram(iy,ups,"fit_sys").c_str(), Form("Yield sys(%dS)",ups), fNpt2,fPTbin2);

	
	TFile *SG_file = new TFile("/Users/carlsonbt1/results/data_files/SG_efficiency.root","READ"); 
	TF1 *SG_eff=(TF1*)SG_file->Get("SG_eff"); 	
	output_files->cd();
	
	for (int ipt=1; ipt<=toy_yield->GetNbinsX(); ipt++) {
		double pt=toy_yield->GetBinCenter(ipt); 
		double eff_SG_vtx=1/(SG_eff->Eval(pt)*0.99);//seagull*vertex efficiency
		
		TH2D *y_E=(TH2D*)output_files->FindObjectAny(Form("yield%dS_E_y%d_pt%d",ups,iy,ipt-1)); 
		if(ipt==1 && ups==1 && iy==0){
			TCanvas * tmp = new TCanvas("","",800,800); 
			y_E->SetTitle(BT(iy,ipt-1).c_str()); 
			y_E->SetStats(kFALSE);
			y_E->Draw("COLZ"); 
			tmp->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/yield_E_y.pdf"); 
			tmp->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/yield_E_y.jpg"); 

			delete tmp; 
		}
		TH1D *sigma_y=(TH1D*)y_E->ProjectionY(); 
		if(ipt==1 && ups==1 && iy==0){
			TCanvas * tmp = new TCanvas("","",800,800); 
			TH1D *sy=(TH1D*)sigma_y->Clone("sy_tmp"); 
			sy->SetTitle(BT(iy,ipt-1).c_str()); 
			sy->SetStats(kFALSE); 
			sy->GetYaxis()->SetTitle("frequency");
			//sigma_y->GetYaxis()->SetTitleSize(0.03); 
			sy->SetAxisRange(sy->GetMean()-3*sy->GetRMS(),sy->GetMean()+3*sy->GetRMS());
			sy->Draw(); 
			tmp->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/yieldE_dist.pdf"); 
			tmp->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/yieldE_dist.jpg"); 

			delete tmp; 
			delete sy;
		}
		
		double mean_sigma=sigma_y->GetMean();//mean fit error
		double RMS_sigma=sigma_y->GetRMS(); //spread in fit errors 
		
		TH1D *y=(TH1D*)y_E->ProjectionX("y",sigma_y->FindBin(mean_sigma-RMS_sigma),sigma_y->FindBin(mean_sigma+RMS_sigma)); 
		
		if(ipt==1 && ups==1 && iy==0){
			TCanvas * tmp = new TCanvas("","",800,800); 
			y->SetTitle(BT(iy,ipt-1).c_str()); 
			y->SetStats(kFALSE); 
			y->GetYaxis()->SetTitle("frequency"); 
		//	y->GetYaxis()->SetTitleSize(0.03); 
			TLatex LT(0.55,0.75,Form("%.1f<#sigma_{Yield}<%.1f",mean_sigma-RMS_sigma,mean_sigma+RMS_sigma));
			LT.SetNDC(kTRUE);
			LT.SetTextSize(0.035); 
			y->Draw(); 
			LT.Draw();
			LT.DrawLatex(0.6,0.71,Form("N_{opt}=%d",n_opt[ipt-1]));
			tmp->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/yield_dist.pdf"); 
			tmp->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/yield_dist.jpg"); 

			delete tmp; 
		}
		
		mean_sigma=mean_sigma*eff_SG_vtx;
		RMS_sigma=RMS_sigma*eff_SG_vtx;
		
		double yield=y->GetMean()*eff_SG_vtx;
		double yield_RMS=y->GetRMS()*eff_SG_vtx;
		double yield_stat=sigma_y->GetMean()*eff_SG_vtx;
		double sys=TMath::Sqrt(yield_RMS*yield_RMS+RMS_sigma*RMS_sigma);//total systematic error including vertex and seagull efficiency

		
		toy_yield->SetBinContent(ipt,yield); //nominal yield 
		toy_yield->SetBinError(ipt,yield_stat);//statistical error  
		toy_RMS->SetBinContent(ipt,sys*100/yield); 
		toy_RMS_stat->SetBinContent(ipt,yield_RMS/mean_sigma); 
		
		double YA=y->GetMean()/eff->GetBinContent(ipt); //actual yield before correction
		
		yield_actual->SetBinContent(ipt,YA); 
		yield_actual->SetBinError(ipt,0);
									
		toy_yield_sys->SetBinContent(ipt,sys); 
		
		delete y; 
		delete sigma_y;
		delete y_E;
		
	}
	output_files->cd("yield_histograms");
	toy_yield->Write(); 
	toy_RMS->Write(); 
	toy_RMS_stat->Write();
	toy_yield_sys->Write();
	
	delete toy_RMS_stat;
	delete toy_yield; 
	delete toy_RMS;
	delete toy_yield_sys; 
	
	SG_file->Close();
	delete SG_file;
	delete SG_eff;

	
}

void fit_sys_error(int iy, int ipt){
	//fit functions
	//TFile *nom =new TFile("output_feb11.root","READ"); 
	
	output_files->cd("yield_fit"); 
	TF1 *FF=(TF1*)output_files->FindObjectAny(fit_function(iy,ipt,"").c_str());
	//TF1 *FFnom=(TF1*)nom->FindObjectAny(fit_function(iy,ipt,"").c_str());

	int Npar=FF->GetNpar();
	
	TF1 *LS_bkg;
	if(Npar==13) LS_bkg=new TF1("LS_bkg", sig_bkg2, 8.65,11.35,Npar);//13 
	else LS_bkg=new TF1("LS_bkg", sig_bkg, 8.65,11.35,Npar); 

	LS_bkg->SetNpx(10000); 

	
	int order=-1; 
	
	for (int i=0; i<Npar; i++) {
		double par=FF->GetParameter(i);
		double parE=FF->GetParError(i);
		LS_bkg->ReleaseParameter(i);
		if(parE==0)LS_bkg->FixParameter(i,par); 
		else LS_bkg->SetParameter(i,par); 
		LS_bkg->SetParName(i,FF->GetParName(i)); 
		TString parName=FF->GetParName(i); 
		TString keyword1("A"); 
		if(parName.Contains(keyword1) && parE==0 && order<0){
			LS_bkg->ReleaseParameter(i);
			order=i; 
		}
	}
	
	LS_bkg->FixParameter(10,8.5);
	LS_bkg->FixParameter(11,11.5);
	
	string fit_method; 
	if(fPTbin2[ipt]>LL_thres) fit_method="LQRIS0"; 
	else fit_method="RQIS0";
	
	TFitResultPtr RR_FULLm=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"FULLM").c_str()); 
	LS_bkg->Write(); 
	
	LS_bkg->SetRange(9.0,10.8);
	TFitResultPtr RR_NARROWm=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"NARROWM").c_str()); 
	LS_bkg->Write(); 
	
	LS_bkg->SetRange(8.7,11.2); 

	
	cout << "Order: " << order << " " << FF->GetParName(order) << endl; 
	cout << "Background order +1: " << endl; 
	
	LS_bkg->SetParameter(order,-1.0); 	

	TFitResultPtr RR_bkgorder_m=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"BKGP_SF").c_str()); 
	LS_bkg->Write(); 
	
	
	LS_bkg->SetRange(8.65,11.35);
	TFitResultPtr RR_bkgsf_FULLm=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"BKGP_SF_FULLM").c_str()); 
	LS_bkg->Write(); 
	
	LS_bkg->SetRange(9.0,10.8);
	TFitResultPtr RR_bkgsf_NARROWm=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"BKGP_SF_NARROWM").c_str()); 
	LS_bkg->Write(); 
	
	
	cout << "Background Order -1: " << endl; 
	LS_bkg->SetRange(8.7,11.2); 
	LS_bkg->FixParameter(order,0); 
	LS_bkg->FixParameter(order-1,0); 
	mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"BKGM_SF").c_str()); 
	LS_bkg->Write(); 
		
	LS_bkg->SetRange(8.65,11.35);
	TFitResultPtr RR_bkgMsf_FULLm=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"BKGM_SF_FULLM").c_str()); 
	LS_bkg->Write(); 
	
	LS_bkg->SetRange(9.0,10.8);
	TFitResultPtr RR_bkgMsf_NARROWm=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"BKGM_SF_NARROWM").c_str()); 
	LS_bkg->Write(); 
	
	LS_bkg->ReleaseParameter(7);
	LS_bkg->FixParameter(8,0);
	LS_bkg->FixParameter(9,0); 
	
	LS_bkg->SetRange(8.8,11.1);
	TFitResultPtr n2junk=mdata->Fit(LS_bkg, fit_method.c_str());
	LS_bkg->SetName(fit_function(iy,ipt,"n2junk").c_str()); 
	LS_bkg->Write(); 
	
	delete LS_bkg; 
	delete FF;
	
}

void compute_chi2(TH1D *h, TF1 *F, double a, double b, double &chi2, int &NDOF){
	//cout << "Compute Chi2: " << endl; 
	chi2=0; 
	NDOF=0; 
	NDOF=(h->FindBin(b)-h->FindBin(a))-F->GetNumberFreeParameters(); 
	//cout << "chi2 Init: " << chi2 << " NDOF: " << NDOF << endl; 
	for (int i=1; i<=h->GetNbinsX(); i++) {
		double m=h->GetBinCenter(i); 
		if(m<a || m>b) continue;  
		double error=1; 
		if(h->GetBinContent(i)>0) error=h->GetBinError(i); 
		else {
			error=TMath::Sqrt(h->GetBinContent(i)+1); 
			NDOF-=1; 
		}

		chi2+=TMath::Power((h->GetBinContent(i)-F->Eval(m))/error,2); 
	}

}

double ups_errors(TF1 *f,TMatrixDSym cor, int ups){
	double Err=0; 
	double S1=f->GetParameter(0); 
	double R12=f->GetParameter(1); 
	double R13=f->GetParameter(2); 
	
	double S1E=f->GetParError(0); 
	double R12E=f->GetParError(1); 
	double R13E=f->GetParError(2);
	

//	cout << "cor[0,0]" << cor[0][0] << endl; 
	
	if(ups==2) {
		//cout << "Error: " << TMath::Power(S1*R12E*cor[1][1],2) + TMath::Power((R12*S1E*cor[0][0]),2) << endl; 
		Err=TMath::Power(S1*R12E*cor[1][1],2) + TMath::Power((R12*S1E*cor[0][0]),2) + 2*S1*R12*S1E*R12E*cor[0][1];
		//cout << "Including Cov: " << Err << endl; 
	}
	if(ups==3) Err=TMath::Power(S1*R13E*cor[2][2],2) + TMath::Power((R13*S1E*cor[0][0]),2) + 2*S1*R13*S1E*R13E*cor[0][2];
	Err=TMath::Sqrt(Err);
	
	
	
	return Err; 
															  
}

void fit(int iy, int ipt){
	double mL=8.7; 
	double mH=11.2;
	
	TF1 *BKG_fixed=new TF1("BKG_fixed",BKG_cheb,mL,mH,7);
	BKG_fixed->FixParameter(4,0);
	if(fPTbin2[ipt]>=pol_thres){
		BKG_fixed->FixParameter(2,0);
		BKG_fixed->FixParameter(3,0);
	}
	double ys = fit_bkg_N(mdata,BKG_fixed,ipt); // fit edges of normalized background
	double BKG_est=mdata->GetBinContent(mdata->FindBin(8.8)); 
	
	BKG_fixed->FixParameter(5,mL);//remove limited region of BKG function to now define over the entire region 
	BKG_fixed->FixParameter(6,mH); 

	int Npar;
	if(fPTbin2[ipt]>=add_params)Npar=14; 
	else Npar=12;

	double dm_scale_max;
	double dm_scale_min;
	if(fPTbin2[ipt]>=add_params){
		dm_scale_max=1.0+dm_scale_width_last[0];
		dm_scale_min=1.0-dm_scale_width_last[1];
	}
	else {
		dm_scale_min=1.0-dm_scale_width; 
		dm_scale_max=1.0+dm_scale_width; 
	}
	
	
    TF1 *LS_bkg;
	if(fPTbin2[ipt]>=add_params) LS_bkg=new TF1("LS_bkg", sig_bkg2, mL,mH,Npar);//13 
	else LS_bkg=new TF1("LS_bkg", sig_bkg, mL,mH,Npar); 
	LS_bkg->SetParameter(4,0);
	if(fPTbin2[ipt]>=add_params){
		LS_bkg->SetParLimits(12,-.1,.1); 
		LS_bkg->SetParameter(12,-0.029); 
		
		LS_bkg->SetParLimits(13,dm_scale_min,dm_scale_max);
		LS_bkg->SetParameter(13,.8);
		LS_bkg->SetParName(12,"deltaM2");
		LS_bkg->SetParName(13,"dm_scale2"); 
	}
	

	LS_bkg->SetParNames("Y1","r12","r13","dm_scale","deltaM");

	LS_bkg->SetNpx(1000); 
	//yields 
	
	double yield_est[]={0,0,0};
	estimate_yield(yield_est);

	LS_bkg->SetParameter(0,yield_est[0]);
	LS_bkg->SetParLimits(0,0,1000);

	LS_bkg->SetParameter(1,yield_est[1]/yield_est[0]);
	LS_bkg->SetParLimits(1,0,1000);

	LS_bkg->SetParameter(2,yield_est[2]/yield_est[0]);
	LS_bkg->SetParLimits(2,0,1000);

	//dmscale / shift
	LS_bkg->SetParameter(3,1); 
	LS_bkg->SetParLimits(3,dm_scale_min,dm_scale_max);
	LS_bkg->SetParLimits(4,-0.1,0.1);
	LS_bkg->SetParameter(4,-0.002);

	LS_bkg->SetParameter(5,BKG_est*ys);
	for (int i=0; i<7; i++) {
		if(i>=4 || (i<4 && i>1 && fPTbin2[ipt]>=pol_thres)) LS_bkg->FixParameter(i+5,BKG_fixed->GetParameter(i)); 
		else LS_bkg->SetParameter(i+5,BKG_fixed->GetParameter(i));
		LS_bkg->SetParName(i+5,BKG_fixed->GetParName(i));
	}

	set_bkg_order(n_opt[ipt],LS_bkg);
	
	string fit_method; 
	if(fPTbin2[ipt]>LL_thres) fit_method="LRS0"; 
	else fit_method="RS0";
	
	cout << "Fit Method: " << fit_method << endl; 
	LS_bkg->SetRange(8.7,11.2);
	TFitResultPtr RR=mdata->Fit(LS_bkg, fit_method.c_str());
	 TMatrixDSym cor = RR->GetCorrelationMatrix();
	//cout << "cov[0][0] " << cor[0][0] << endl;
	RR->Print("V");
	
cout << "ups error: " << ups_errors(LS_bkg,cor,2) << endl; 
cout << "ups error: " << ups_errors(LS_bkg,cor,3) << endl;
	
	
	
	//cout << "Fit Status: " << gMinuit->fStatus << endl; 
	cout << "fCStatu: " <<  gMinuit->fCstatu << endl; 

	
	double bin_center1S = get_bin_center(ipt,1); 
	double bin_center2S = get_bin_center(ipt,2);
	double bin_center3S = get_bin_center(ipt,3); 
	
	deltaM_pt->Fill(bin_center1S,1000*LS_bkg->GetParameter(4));
	dm_scale_pt->Fill(bin_center1S,LS_bkg->GetParameter(3));
	
	deltaM_pt->SetBinError(deltaM_pt->FindBin(bin_center1S),1000*LS_bkg->GetParError(4)); 
	dm_scale_pt->SetBinError(dm_scale_pt->FindBin(bin_center1S),LS_bkg->GetParError(3));
	
	S1_pt->Fill(bin_center1S,LS_bkg->GetParameter(0)); 
	S1_pt->SetBinError(S1_pt->FindBin(bin_center1S), LS_bkg->GetParError(0)); 
	
	S2_pt->Fill(bin_center2S,LS_bkg->GetParameter(1)); 
	S2_pt->SetBinError(S2_pt->FindBin(bin_center2S), LS_bkg->GetParError(1)); 

	S3_pt->Fill(bin_center3S,LS_bkg->GetParameter(2)); 
	S3_pt->SetBinError(S3_pt->FindBin(bin_center3S), LS_bkg->GetParError(2)); 
	
	cout << "Compute Residuals: " << endl; 
	
	int NDF; 
	
	double chi2=compute_resid(0,ipt,mdata,LS_bkg,NDF);
	double NDOF=static_cast<double>(NDF);
	double chi2_NDF=chi2/NDOF;
	
	cout << "Over 8.7-11.2: " << endl; 
	cout << "Chi2: " << LS_bkg->GetChisquare() << endl; 
	cout << "NDOF: " << LS_bkg->GetNDF() << endl; 
	/*
	LS_bkg->SetRange(9.1,10.6); 
	cout << "Over Restricted Mass Range: " << endl; 
	cout << "Chi2: " << LS_bkg->GetChisquare() << endl; 
	cout << "NDF: " << LS_bkg->GetNDF() << endl; 
	 */
	/*
	r12->Fill(bin_center1S,LS_bkg->GetParameter(1)); 
	r13->Fill(bin_center1S,LS_bkg->GetParameter(2));
	r12->SetBinError(r12->FindBin(bin_center1S),LS_bkg->GetParError(1));
	r13->SetBinError(r13->FindBin(bin_center1S),LS_bkg->GetParError(2));
*/
	chi2_pt->Fill(bin_center1S,chi2_NDF); 
	
	TH1D *mdata_clone = (TH1D*)mdata->Clone(fit_histogram(iy,ipt).c_str());
	/*
	cout << "My method: " << endl; 
	compute_chi2(mdata_clone, LS_bkg, 9.1,10.6,chi2,NDF);
	cout << "Chi2: " << chi2<< endl; 
	cout << "NDOF: " << NDF << endl; 
	*/
	LS_bkg->SetName(fit_function(iy,ipt,"").c_str());
	output_files->cd("yield_fit");
	mdata_clone->Write();
	LS_bkg->Write();
	
	//delete pointers
	delete mdata_clone;
	delete LS_bkg; 
	delete BKG_fixed;
	

}

double resolutionP(TF1 *signal, int ups ){

	double pm=PDG_mass[ups-1];
	double step=0.001; // 1MeV step size 
	double upper=pm+step; // upper edge
	double tot_area = signal->Integral(8.7,11.); 
	for (int i=0; i<250; i++) {
		double area=signal->Integral(pm,upper)/tot_area; 
		//cout << "area: " << area << endl; 
		if(area>=0.34) continue; 
		upper+=step; 
	}
	return TMath::Abs(upper-pm); 
}

double resolutionM(TF1 *signal, int ups ){
	
	double pm=PDG_mass[ups-1];
	double step=0.001; // 1MeV step size 
	double lower=pm-step; // upper edge
	double tot_area = signal->Integral(8.7,11.); 
	for (int i=0; i<250; i++) {
		double area=signal->Integral(lower,pm)/tot_area; 
		//cout << "area: " << area << endl; 
		if(area>=0.34) continue; 
		lower-=step; 
	}
	return TMath::Abs(lower-pm); 
}

void make_sample_LS(){

	TF1 *sampleLS1=new TF1("sampleLS",signal_shape,8.7,11.2,4);
	TF1 *sampleLS2=new TF1("sampleLS",signal_shape,8.7,11.2,4);
	TF1 *sampleLS3=new TF1("sampleLS",signal_shape,8.7,11.2,4);

	sampleLS1->SetNpx(10000);
	sampleLS2->SetNpx(10000);
	sampleLS3->SetNpx(10000);
	sampleLS1->SetParameters(1,1,1.0,0.0); 
	
	TCanvas *LScanvas = new TCanvas("LScanvas","",800,800);
	LScanvas->cd();
	sampleLS1->SetLineWidth(1);
	sampleLS2->SetLineWidth(1);
	sampleLS3->SetLineWidth(1);

	sampleLS1->SetLineStyle(1);
	sampleLS2->SetLineStyle(2);
	sampleLS3->SetLineStyle(4);
	
	sampleLS1->Draw(); 	
	sampleLS2->Draw();
	sampleLS3->Draw();
	
	sampleLS1->Draw();
	sampleLS1->GetHistogram()->GetXaxis()->SetTitle("M_{#mu#mu} [GeV]"); 
	sampleLS1->GetHistogram()->GetYaxis()->SetTitle("arb. units"); 

	LScanvas->Modified();	
	sampleLS2->SetParameters(1,1,1.2,0); 
	sampleLS2->SetLineColor(kRed); 
	sampleLS2->Draw("same");
	sampleLS3->SetParameters(1,1,1.,-.1); 
	sampleLS3->SetLineColor(kBlue);
	sampleLS3->Draw("same");
	
	TLegend *L = new TLegend(0.55,0.4,0.85,0.6,"","NDC"); 
	L->SetFillColor(10);
	L->SetLineColor(10); 
	
	L->AddEntry(sampleLS1,"#delta m=0, c_{w}=1.0","L"); 
	L->AddEntry(sampleLS2,"#delta m=0, c_{w}=1.2","L"); 
	L->AddEntry(sampleLS3,"#delta m=-100 MeV, c_{w}=1.0","L"); 
	
	prelimText.SetNDC(kTRUE); 
	lumiText.SetNDC(kTRUE);
	prelimText.DrawLatex(label_x+0.15,label_y,cms_pre);
	lumiText.DrawLatex(label_x+0.15,label_y-label_text_size,lumi_string);
	L->Draw("same"); 
	LScanvas->Modified();
	
	LScanvas->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/sample_LS.pdf");
	LScanvas->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/sample_LS.jpg");


	
	delete sampleLS1;
	delete sampleLS2;
	delete sampleLS3;
	delete L;
	delete LScanvas; 
	
	
}

void fill_yields(int iy, string mode){
	cout << "Fill Yields: " << mode << endl; 
	
	string name; 
	
	TH1D *yield1S = new TH1D(yield_histogram(iy,1,mode.c_str()).c_str(), "yield1S", fNpt2,fPTbin2); 
	TH1D *yield2S = new TH1D(yield_histogram(iy,2,mode.c_str()).c_str(), "yield2S", fNpt2,fPTbin2); 
	TH1D *yield3S = new TH1D(yield_histogram(iy,3,mode.c_str()).c_str(), "yield3S", fNpt2,fPTbin2); 
	
	TString MM(mode.c_str()); 
	
	for (int ipt=1; ipt<=yield1S->GetNbinsX(); ipt++) {
		if(fPTbin2[ipt-1]>=pol_thres && MM.Contains("BKGM")){
			//cout << "pT: " << fPTbin2[ipt-1] << endl; 
			yield1S->SetBinContent(ipt,0); 
			yield2S->SetBinContent(ipt,0); 
			yield3S->SetBinContent(ipt,0);
			yield1S->SetBinError(ipt,0); 
			yield2S->SetBinError(ipt,0); 
			yield3S->SetBinError(ipt,0);

			continue;
		}
		
		double BW=((TH1D*)output_files->FindObjectAny(fit_histogram(iy,ipt-1).c_str()))->GetBinWidth(1);
		string name = fit_function(iy,ipt-1,mode); 
		
		TF1 *LS_bkg=(TF1*)output_files->FindObjectAny(name.c_str());
		
		yield1S->SetBinContent(ipt,LS_bkg->GetParameter(0)/BW); 
		yield2S->SetBinContent(ipt,LS_bkg->GetParameter(1)*LS_bkg->GetParameter(1)/BW); 
		yield3S->SetBinContent(ipt,LS_bkg->GetParameter(2)*LS_bkg->GetParameter(2)/BW); 
		yield1S->SetBinError(ipt,LS_bkg->GetParError(0)/BW); 
		double R1E=LS_bkg->GetParError(0);
		double R2E=LS_bkg->GetParError(1);
		double R3E=LS_bkg->GetParError(2); 
		
		
		yield2S->SetBinError(ipt,TMath::Sqrt(R1E*R1E+R2E*R2E)/BW); // used to just use LS_bkg->GetParError
		yield3S->SetBinError(ipt,TMath::Sqrt(R2E*R2E+R3E*R3E)/BW); 

		delete LS_bkg;

	
	}//for loop over ipt
	//Scale yields by Seagull Efficiency
	
	TFile *SG_file = new TFile("/Users/carlsonbt1/results/data_files/SG_efficiency.root","READ"); 
	TF1 *SG_eff=(TF1*)SG_file->Get("SG_eff"); 
	cout << "Divide by Seagull Efficiency. " << endl; 
	for(int ipt=1; ipt<=yield1S->GetNbinsX(); ipt++){
		double pt=yield1S->GetBinCenter(ipt); 
		double eff_SG_vtx=1/(SG_eff->Eval(pt)*0.99);
	
		double y1=eff_SG_vtx*yield1S->GetBinContent(ipt); 
		double y1E=eff_SG_vtx*yield1S->GetBinError(ipt); 
		yield1S->SetBinContent(ipt,y1); 
		yield1S->SetBinError(ipt,y1E); 
		
		double y2=eff_SG_vtx*yield2S->GetBinContent(ipt); 
		double y2E=eff_SG_vtx*yield2S->GetBinError(ipt); 
		yield2S->SetBinContent(ipt,y2); 
		yield2S->SetBinError(ipt,y2E);
		
		double y3=eff_SG_vtx*yield3S->GetBinContent(ipt); 
		double y3E=eff_SG_vtx*yield3S->GetBinError(ipt); 
		yield3S->SetBinContent(ipt,y3); 
		yield3S->SetBinError(ipt,y3E); 
		
	}
	
	cout << "Done Filling Yields: " << endl; 
	output_files->cd("yield_histograms");
	yield1S->Write();
	yield2S->Write();
	yield3S->Write();
	
	delete yield1S; 
	delete yield2S;
	delete yield3S;
	
	delete SG_eff; 
	delete SG_file; 
	
}

void make_fit_plot(int iy, int ipt, string mode){
	
	double mL=8.7;
	double mH=11.2;

	TLatex prelimTextM(label_x+0.1, label_y, cms_pre);
	TLatex lumiTextM(label_x+0.1, label_y-label_text_size, lumi_string);
	
	prelimTextM.SetTextSize(label_text_size*1.25); 
	lumiTextM.SetTextSize(label_text_size*1.25);
	
	string name=fit_function(iy,ipt,mode);
	//cout << "Get fit function: " << name << endl; 
	
	TF1 *LS_bkg=(TF1*)output_files->FindObjectAny(name.c_str());
	
	if(ipt==0 && mode=="") make_sample_LS();
	
	
	TF1 *signal1S = new TF1("signal1S",signal_shape,8.7,11.2,4); 
	signal1S->SetParNames("Ups","Yield 1S","dm_scale","deltaM");
	signal1S->SetNpx(10000);
	signal1S->SetParameters(1,1,LS_bkg->GetParameter(3),LS_bkg->GetParameter(4)); 
	//cout << "Integral of normalized shape: " << signal1S->Integral(8.7,11.2) << endl; 
	signal1S->SetParameters(1,LS_bkg->GetParameter(0),LS_bkg->GetParameter(3),LS_bkg->GetParameter(4)); 
	signal1S->SetLineColor(kRed);
	
	//cout << "Yield Integral: " << signal1S->Integral(8.7,11.2)/mdata->GetBinWidth(1) << endl; 
		
	TF1 *signal2S = new TF1("signal2S",signal_shape,mL,mH,4); 
	signal2S->SetParNames("Ups","Yield 2S","dm_scale","deltaM");
	signal2S->SetNpx(1000);
	int par[]={3,4};
	if(fPTbin2[ipt]>=add_params){
		par[0]=13;
		par[1]=12;
	}
	
	signal2S->SetParameters(2,LS_bkg->GetParameter(1)*LS_bkg->GetParameter(0),LS_bkg->GetParameter(par[0]),LS_bkg->GetParameter(par[1])); 
	signal2S->SetLineColor(kBlue);
	
	TF1 *signal3S = new TF1("signal3S",signal_shape,mL,mH,4); 
	signal3S->SetParNames("Ups","Yield 3S","dm_scale","deltaM");
	signal3S->SetNpx(1000);
	signal3S->SetParameters(3,LS_bkg->GetParameter(2)*LS_bkg->GetParameter(0),LS_bkg->GetParameter(par[0]),LS_bkg->GetParameter(par[1])); 
	signal3S->SetLineColor(kGreen);
	
	TF1 *bkg = new TF1("bkg",BKG_cheb,mL,mH,7);
	bkg->SetParameter(5,mL); 
	bkg->SetParameter(6,mH); 
	bkg->SetParNames("A0","A1","A2","A3","A4","no exclusion","no exclusion"); 
	for(int i=0; i<5;i++)  
		bkg->SetParameter(i,LS_bkg->GetParameter(i+5)); 				
	bkg->SetLineColor(kBlue); 
	bkg->SetLineStyle(kDotted); 
	
	TLegend *L = new TLegend(0.5,0.5,0.9,0.9); 
	L->SetHeader(summary_plot(iy,ipt).c_str()); 
	
	if(mode==""){
		res_pt_P1S->Fill(get_bin_center(ipt,1),resolutionP(signal1S,1)*1000); 
		res_pt_M1S->Fill(get_bin_center(ipt,1),resolutionM(signal1S,1)*1000); 
	}
	
	
	TCanvas *C = new TCanvas("Mass Shape");
	string plot_sumN=summary_plot(iy,ipt).c_str()+mode;
	C->SetName(plot_sumN.c_str());
	C->cd();
	TH1D *mdata_clone =(TH1D*)output_files->FindObjectAny(fit_histogram(iy,ipt).c_str());
	mdata_clone->SetStats(kFALSE); 
	mdata_clone->SetMarkerStyle(20);
	mdata_clone->SetMarkerSize(0.4);
	mdata_clone->SetAxisRange(8.65,11.35);
	mdata_clone->SetTitle(BT(iy,ipt).c_str()); 
	mdata_clone->GetXaxis()->SetTitle("M_{#mu#mu} [GeV]"); 
	mdata_clone->GetYaxis()->SetTitle(Form("Events/%.0f MeV",mdata_clone->GetBinWidth(1)*1000)); 	
	
	LS_bkg->SetLineWidth(1);
	signal1S->SetLineWidth(1);
	signal2S->SetLineWidth(1);
	signal3S->SetLineWidth(1);
	bkg->SetLineWidth(1);
	
	mdata_clone->Draw("E1P");
	LS_bkg->Draw("same");
	signal1S->Draw("same"); 
	signal2S->Draw("same"); 
	signal3S->Draw("same");
	bkg->Draw("same"); 
	prelimTextM.SetNDC(kTRUE); 
	lumiTextM.SetNDC(kTRUE);
	prelimTextM.Draw();
	lumiTextM.Draw();
	lumiTextM.DrawLatex(label_x+0.1,label_y-2.6*label_text_size,BT(iy,ipt).c_str()); 
	
	
	output_files->cd("summary_plots");
	string file_name="/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/"+summary_plot(iy,ipt)+mode+".pdf";
	C->Print(file_name.c_str());
	file_name="/Users/carlsonbt1/results/AN_Fig/UpsAN/FigJPG/fits/"+summary_plot(iy,ipt)+mode+".jpg";
	C->Print(file_name.c_str());
	C->Write();
	
	if(ipt==0 && mode==""){
		L->AddEntry(mdata_clone,"Data","lep");
		L->AddEntry(bkg, "Background"); 
		L->AddEntry(signal1S, "Signal #Upsilon(1S)"); 
		L->AddEntry(signal2S, "Signal #Upsilon(2S)"); 
		L->AddEntry(signal3S, "Signal #Upsilon(3S)"); 
		L->Draw();
		prelimText.SetNDC(kTRUE); 
		lumiText.SetNDC(kTRUE);
		prelimText.Draw();
		lumiText.Draw();
		
		C->Print("/Users/carlsonbt1/results/AN_Fig/UpsAN/Fig/fits/Legend.pdf"); 
	}
	
	if(ipt+1<=16) CName["fits"]->cd(ipt+1); 
	else CName["fits2"]->cd(ipt+1-16); 
	mdata_clone->SetTitle(BT(iy,ipt).c_str());
	mdata_clone->DrawCopy("E1"); 
	LS_bkg->DrawCopy("same"); 
	signal1S->DrawClone("same"); 
	signal2S->DrawClone("same"); 
	signal3S->DrawClone("same");
	bkg->DrawClone("same"); 
	
	if(ipt==fNpt2-1) {
		L->SetHeader(""); 
		CName["fits2"]->cd(9);
		L->DrawClone("same");
	}
	
	delete signal1S;
	delete signal2S;
	delete signal3S; 
	delete mdata_clone;
	delete C; 
	delete L;
	
}

void get_dm_m(int iy, int ipt, string mode){
	//Get dm_m histogram for a given bin, in some mode. example modes are:
	//"W" - weighted
	//"NW" - no weight
	//WEP - Weighted-error positive
	//weighted - error neg
	
	
	string hist_name; 
	hist_name=dm_m_hist_name_weighted(iy,ipt,mode); 
	
	//cout << "Getting Histogram: " << hist_name << endl; 
	
	dm_m_data=(TH2D*)data->FindObjectAny(hist_name.c_str()); 
	
	string hist_name_unweighted;
	hist_name_unweighted=dm_m_hist_name(iy,ipt); 
	dm_m_data_unweighted=(TH2D*)data->FindObjectAny(hist_name_unweighted.c_str());
	
}

void get_LS(int iy, int ipt, string mode){
	string LS_name1; 
	string LS_name2;
	string LS_name3; 
	
	if(mode=="W"){
		LS_name1=LS_namew(1, iy,ipt); 
		LS_name2=LS_namew(2, iy,ipt);
		LS_name3=LS_namew(3, iy,ipt);
	}
	
	//LS1_shape=(TF2*)lineshape_file->FindObjectAny(LS_name1.c_str());
	
	LS1_shape=(TF2*)lineshape_file->FindObjectAny(LS_name1.c_str()); 
	LS2_shape=(TF2*)lineshape_file->FindObjectAny(LS_name2.c_str());
	LS3_shape=(TF2*)lineshape_file->FindObjectAny(LS_name3.c_str());
}

double adjust_uncertainties(TH1D *h, TH1D *hun_w){

	double w=h->Integral()/hun_w->Integral(); 
	//cout <<"w: " << w << endl; 
	for (int im=1; im<=h->GetNbinsX(); im++) {
		double Sw=h->GetBinContent(im); 
		
		double N=hun_w->GetBinContent(im); 
		double sigmaN=hun_w->GetBinError(im); 
		
		double correctedError=sigmaN*w;
		
	//	if(h->FindBin(9.46)==im) cout << "Sw: " << Sw << " N: " << N << " sigma " << correctedError << endl; 
		
		h->SetBinError(im, correctedError); 
		
	}
	return w;
	
}

void make_mdata(int iy, int ipt){
	//Get mass PDF for data '
	//cout << " Make data. " << endl; 
	TH1D *dm =(TH1D*)dm_m_data->ProjectionY("dm"); 
	//cout << "Projection." << endl; 
	TH1D *mdata_tmp = (TH1D*)dm_m_data->ProjectionX("mdata_tmp",0,dm->FindBin(dm_max)); 
	TH1D *mdata_unw = (TH1D*)dm_m_data_unweighted->ProjectionX("mdata_un",0,dm->FindBin(dm_max)); 
	
	cout << "Weighted Integral: " << mdata_tmp->Integral(mdata_tmp->FindBin(9.26),mdata_tmp->FindBin(9.66)) << endl; 
	cout << "Unweighted Integral:" << mdata_unw->Integral(mdata_unw->FindBin(9.26),mdata_unw->FindBin(9.66)) << endl; 
	
	cout << "Weighted Integral: " << mdata_tmp->Integral(mdata_tmp->FindBin(y2m[0]),mdata_tmp->FindBin(y2m[1])) << endl; 
	cout << "Unweighted Integral:" << mdata_unw->Integral(mdata_unw->FindBin(y2m[0]),mdata_unw->FindBin(y2m[1])) << endl; 
	
	double w=adjust_uncertainties(mdata_tmp,mdata_unw); 
	eff->SetBinContent(ipt+1,w); 
	eff->SetBinError(ipt+1,0);
	
	if(fPTbin2[ipt]>rebin_thres && fPTbin2[ipt]<40) mdata_tmp->Rebin(rebin); 
	if(fPTbin2[ipt]>=40 && fPTbin2[ipt]<50) mdata_tmp->Rebin(3);
	if(fPTbin2[ipt]>=50 && fPTbin2[ipt]<70) mdata_tmp->Rebin(4);
	if(fPTbin2[ipt]>=70) mdata_tmp->Rebin(6);
		
	mdata=(TH1D*)mdata_tmp->Clone("mdata");
	
	delete mdata_tmp;
	delete mdata_unw; 
	delete dm; 

}

