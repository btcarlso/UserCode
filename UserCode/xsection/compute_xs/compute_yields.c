/*
 *  compute_yields.c
 *  
 *
 *  Created by Benjamin Carlson on 8/7/13.
 *  Copyright 2013 Carnegie Mellon University. All rights reserved.
 *
 */

#include "compute_yields.h"


void compute_yields(){
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); // get rid of messages from roofit below warning level significance. 
	RooMsgService::instance().setSilentMode(kTRUE); // get rid of messages from roofit below warning level significance. 

	gROOT->SetBatch();
	/*
	TFile *data_file = new TFile("/Users/carlsonbt1/results/tmp/onia2MuMu_tree_11_1_dgY.root","READ"); 
	if(data_file->IsOpen()!=1) return; 
	OniaTree=(TTree*)data_file->FindObjectAny("data"); 
	*/
	
	hist_file = new TFile("hist_file.root","READ"); 
	yields_file = new TFile("yields_file.root","RECREATE"); 
	readHisto();
	fTree= (TTree*)hist_file->Get("fTree"); 
	bookHisto();
	
//	for(int iy=0; iy<fNy; iy++){
//	
	//int iy=0; 
	for(int iy=0; iy<2; iy++){
		for(int ipt=0; ipt<fNpt; ipt++){
			iY=iy; 
			iPT=ipt; 
			MODE="";
			
			bin_weight=1./prName[eff_name(iy,1,"")]->GetBinContent(ipt+1);
			cout << Form("Efficiency weight[ipt=%d]",ipt) << bin_weight << " +/- " << prName[eff_name(iy,1,"")]->GetBinError(ipt+1) << endl;
			
			if(ipt!=0) continue; 
			get_LS(iy,ipt);
			if(ipt==0 && iy==0)  make_sample_LS(); 
			get_data(iy,ipt,""); // get histograms / unbinned data 

			fit(); 

		}
}
	//average_yields(iy);
	//print_table1();
	writeHisto();
	
}

void readHisto(){
	cout << "Read Histograms: " << endl; 
	TH1F *TH1F_names = (TH1F*)hist_file->FindObjectAny("TH1F_names"); 
	TH1F *TH2F_names = (TH1F*)hist_file->FindObjectAny("TH2F_names"); 
	TH1F *TProfile_names = (TH1F*)hist_file->FindObjectAny("TProfile_names"); 
	cout << "list of names acquired: " << endl; 
	
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		cout << name << endl; 
	}
	cout << "acquiring..." << endl; 
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names" || name.Contains("rho_dRPtE")) continue; 
		TH1F *h=(TH1F*)hist_file->FindObjectAny(name); 
		cout << "Loading histogram: " << h->GetName() << endl; 
		hName[h->GetName()]=h; 
	}
	
	for (int i=1; i<=TProfile_names->GetNbinsX(); i++) {
		TString name=TProfile_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names"|| name=="TProfile_names") continue; 
		TProfile *pr=(TProfile*)hist_file->FindObjectAny(name); 
		cout << "Loading profile: " << pr->GetName() << endl; 
		prName[pr->GetName()]=pr; 
	}
	/*
	for (int i=1; i<=TH2F_names->GetNbinsX(); i++) {
		TString name=TH2F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names") continue; 
		TH2F *h=(TH2F*)hist_file->FindObjectAny(name); 
		hName2D[h->GetName()]=h; 
	}*/
	
	
}

void bookHisto(){
	double xLow, yLow; 
	double xHigh, yHigh; 
	double BW;
	int nBinsX, nBinsY;
	TString name, title, xLabel, yLabel;
	title="";
	for(int iy=0; iy<fNy; iy++){
			for(int ups=1; ups<=3; ups++){
				for(int order=0; order<=5; order++){
					name=yield_histogram(iy,ups,Form("T%d",order),""); 
					xLabel = "p_{T} [GeV]";
					yLabel = "Yield";
					
					CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
					
					name=ratio_histogram(iy,ups,Form("T%d",order),""); 
					yLabel="Ratio"; 
					if(ups>1)CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
				
				}//loop over orders 


				name=yield_histogram(iy,ups,"best",""); 
				xLabel = "p_{T} [GeV]";
				yLabel = "Yield";
				CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
				
				name=ratio_histogram(iy,ups,"best",""); 
				yLabel="Ratio"; 
				if(ups>1)CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
				
				name=background_order_systematic(iy,Form("%dS",ups)); 
				yLabel="Yield Systematic, background order"; 
				CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
				

			}//ups loop
		
		
		name=background_order(iy,""); 
		//cout << "Create Histogram: "<< name << endl; 
		yLabel="Chebychev order";
		CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
		name=background_order(iy,"w"); 
		//cout << "Create Histogram: " << name << endl; 
		CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 		
		
		name=background_order_systematic(iy,"R21"); 
		yLabel="Yield Systematic, background order"; 
		CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
		
		name=background_order_systematic(iy,"R31"); 
		yLabel="Yield Systematic, background order"; 
		CreateHistogram(name,title,xLabel,yLabel,fNpt,fPTbin); 
		
	}//iy loop 
	cout << "Create unbinned histograms: " << endl; 
	name="order_probability_"; 
	xLabel="Probability";
	yLabel="Fits"; 
	CreateHistogram(name, title,xLabel,yLabel,20,0,1); 

}

void writeHisto(){
	yields_file->cd(); 
	yields_file->mkdir("graphs");
	yields_file->mkdir("histograms");
	yields_file->mkdir("plots");

	yields_file->cd("histograms"); 
	for (std::map<TString,TH1F*>::iterator it=h_Name.begin(); it!=h_Name.end(); it++) {
		h_Name[it->first]->Write();
		//out << "TH1F *h=(TH1F*)hist_file.FindObjectAny(" << it->first << ")"; 
		
	}
	yields_file->cd("graphs"); 
	for (std::map<TString,TGraph*>::iterator it=grName.begin(); it!=grName.end(); it++) {
		grName[it->first]->Write();
		//out << "TH1F *h=(TH1F*)hist_file.FindObjectAny(" << it->first << ")"; 
		
	}
	yields_file->cd("plots"); 
	for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
		TString plot_name =it->first+".pdf";
		if(it->first.Contains("fit_summary"))  plot_name="fits/"+it->first + ".pdf";  
		CName[it->first]->Print("/uscms_data/d3/btcarlso/UpsAN/Figures/"+plot_name);
		//out << "TH1F *h=(TH1F*)hist_file.FindObjectAny(" << it->first << ")"; 
		
	}
	
	
	
}


void make_sample_LS(){
	cout << "make sample LS " << endl; 
	TF1 *sampleLS1=new TF1("sampleLS1",signal_shape1,8.7,11.2,2);
	TF1 *sampleLS2=new TF1("sampleLS2",signal_shape1,8.7,11.2,2);
	TF1 *sampleLS3=new TF1("sampleLS3",signal_shape1,8.7,11.2,2);
	cout << "set parameters: " << endl; 
	sampleLS1->SetNpx(10000);
	sampleLS2->SetNpx(10000);
	sampleLS3->SetNpx(10000);
	sampleLS1->SetParameters(0.0,1.0);  
	sampleLS2->SetParameters(0.0,1.2); 
	sampleLS3->SetParameters(-0.1, 1.0); 

	
	cout << "Create canvas and set line styles: " << endl; 
	CreateCanvas("LSCanvas","",700,600); 
	CName["LSCanvas"]->cd();
	sampleLS1->SetLineWidth(1);
	sampleLS2->SetLineWidth(1);
	sampleLS3->SetLineWidth(1);
	
	sampleLS1->SetLineStyle(1);
	sampleLS2->SetLineStyle(2);
	sampleLS3->SetLineStyle(4);
	
	sampleLS1->SetLineColor(kBlack); 
	sampleLS2->SetLineColor(kRed); 
	sampleLS3->SetLineColor(kBlue);
	
	cout << "Fill legend: " << endl; 
	
	TLegend *L = new TLegend(0.55,0.4,0.85,0.6,"","NDC"); 
	L->SetFillColor(10);
	L->SetLineColor(10); 
	
	L->AddEntry(sampleLS1,"#delta m=0, c_{w}=1.0","L"); 
	L->AddEntry(sampleLS2,"#delta m=0, c_{w}=1.2","L"); 
	L->AddEntry(sampleLS3,"#delta m=-100 MeV, c_{w}=1.0","L"); 
	
	cout << "Draw: " << endl; 
	
	sampleLS1->Draw(); 	
	sampleLS2->Draw();
	sampleLS3->Draw();
	
	sampleLS1->Draw();
	sampleLS1->GetHistogram()->GetXaxis()->SetTitle("M_{#mu#mu} [GeV]"); 
	sampleLS1->GetHistogram()->GetYaxis()->SetTitle("arb. units"); 

	sampleLS1->DrawClone();
	sampleLS2->DrawClone("same");
	sampleLS3->DrawClone("same");
	

	
	L->DrawClone("same"); 
	draw_header();
	
	
	delete sampleLS1;
	delete sampleLS2;
	delete sampleLS3;
	
	
}


void average_yield(double *par, double *parE, double *prob, double *average, string variable_name){
	bool print=false; 
	if(print)cout << "average yields: " << endl; 
	
	ofstream output(Form("/uscms_data/d3/btcarlso/UpsAN/tables/table_weighting_%s_y%d_pt%d.tex",variable_name.c_str(),iY,iPT)); 
	output <<"\\begin{tabular}{"; 
	int N=Ncheb; 
	for(int i=0; i<=N; i++) output << "c"; 
	output << "}\\hline" << endl;
	output << "order & Y & $\\sigma_{Y}$ & prob \\\\" <<endl; 
	
	
	double Mean=0; 
	double sigma=0; 
	double sys=0; 
	for (int order=0; order<=Ncheb; order++) {
		Mean+=par[order]*prob[order]; 
		sigma+=parE[order]*prob[order]; 
		output  << order << " & " << par[order] << " & " << parE[order] << " & " << prob[order] << " \\\\" << endl; 
	}	
	for (int order=0; order<=4; order++) {
		sys+=TMath::Power(prob[order]*(par[order]-Mean),2); 
	}
	output << "hline" << endl; 
	output << "Average & $\\tilde{Y}$ & $\\sigma_{\\tilde{Y}}$ & $\\sigma_{sys}$ \\\\ hline" << endl; 
	sys=TMath::Sqrt(sys); 
	average[0]=Mean;
	average[1]=sigma; 
	average[2]=sys; 

	output  << Mean << " & " << sigma << " & " << sys << "\\\\" << endl; 
	
	output.close(); 
}

void estimate_yield(TH1F *h){
	
	yield_[0]=h->Integral(h->FindBin(y1m[0]),h->FindBin(y1m[1]));
	yield_[1]=h->Integral(h->FindBin(y2m[0]),h->FindBin(y2m[1]));
	yield_[2]=h->Integral(h->FindBin(y3m[0]),h->FindBin(y3m[1]));
	
}
void get_data(int iy, int ipt,string mode){
	bool print = false; 
	if(print) cout << "Get data iy: "<< iy << " ipt: "<< ipt << " mode: " << mode << endl; 
	RooRealVar pt("UPt","p_{T} [GeV]",fPTbin[ipt],fPTbin[ipt+1]);
	RooRealVar yp("URapidity","y",fYbin[iy],fYbin[iy+1]); 
	RooRealVar ym("URapidity","y",-fYbin[iy+1],-fYbin[iy]); 

	if(print){
		yp.Print("V");
		ym.Print("V"); 
	}
	RooRealVar eff("Eff","Eff",0,1); 
	

	RooDataSet *data(0);
	if(mode=="") {
		RooDataSet tmp("tmp",dataSet_name(iy,ipt,mode),RooArgSet(m,pt,ym),Import(*fTree)); //make unbinned data set
		data= new RooDataSet(dataSet_name(iy,ipt,mode),dataSet_name(iy,ipt,mode),RooArgSet(m,pt,yp),Import(*fTree)); //make unbinned data set
		data->append(tmp); 
	}
	if(mode=="w"){
		RooDataSet tmp("tmp",dataSet_name(iy,ipt,mode),fTree,RooArgSet(m,pt,ym,eff),0,eff.GetName());
		data= new RooDataSet(dataSet_name(iy,ipt,mode),dataSet_name(iy,ipt,mode),fTree,RooArgSet(m,pt,yp,eff),0,eff.GetName()); //make unbinned data set
		data->append(tmp); 
	}
	if(print) data->Print("V"); 
	
	dataSet[dataSet_name(iy,ipt,mode)]=data; 
	
	if(print)cout << "Get Histogram: " << h_m(iy,ipt,mode) << endl; 
	
	m.setBins((11.3-8.7)/hName[h_m(iy,ipt,mode)]->GetBinWidth(1)); 
	RooDataHist *dHist = new RooDataHist(dataHist_name(iy,ipt,mode),dataHist_name(iy,ipt,mode),m,Import(*hName[h_m(iy,ipt,mode)]));
	if(print)dHist->Print("V"); 
	dataHist[dataHist_name(iy,ipt,mode)]=dHist; 
	
}//get data

int compute_AIC_prob(double *AIC, double *P){
	bool print = false; 
	if(print) cout << "compute_AIC prob: " << endl; 
	//int N=sizeof(AIC)/sizeof(double); //size of array of AIC variables 
	int N=Ncheb; 
	int minindex=TMath::LocMin(N,AIC); 
	double min=TMath::MinElement(N, AIC); 
	if(print) cout << "N models evaluated: " << N << endl; 
	if(print) cout << "min: " <<  min << endl; 

	double totalweight=0; 
	for(int i=0; i<N; i++){
		double DeltaI=(AIC[i]-min);
		double prob=TMath::Exp(-DeltaI/2);
		totalweight+=prob; 
	}
	
	for(int i=0; i<N; i++){
		double DeltaI=(AIC[i]-min);		
		double prob=TMath::Exp(-DeltaI/2)/totalweight;
		P[i]=prob; 
		if(print) cout << "weight: " << prob <<  endl; 
		if(MODE=="")h_Name["order_probability_"]->Fill(prob);
	}
	return minindex; 
}

/*
void pt_bin(int iy, int ipt, string mode){
	iY=iy; 
	iPT=ipt; 
	MODE=mode; 

	
	//RooFormulaVar wFunc("w","event weight",Form("%.1f",weight),m); 
	
	RooDataSet data("data","data",RooArgSet(m,pt,y),Import(*fTree)); //make unbinned data set
	//RooRealVar *w = (RooRealVar*)data.addColumn(wFunc); 
	
	//RooDataSet wdata(data.GetName(), data.GetTitle(),&data,*data.get(),0,w->GetName()); 
	
//	wdata.Print(); 
	

	
	estimate_yield(hName[h_m(iy,ipt,mode)]); 
	
	RooDataHist data_binned("data_binned","data",m,Import(*hName[h_m(iy,ipt,mode)])); 
//	RooDataHist data_binned("data_binned","data",RooArgSet(m),wdata); //make unbinned data set
	//data_binned.Print();
	RooPlot *frame = m.frame(Title("Invariant Mass")); 
	
	data_binned.plotOn(frame); 
//	if(iPT!=10) return;
	fit(data_binned,frame); 
	CreateCanvas(Form("Fit_y%d_pt%d_%s",iY,iPT,MODE.c_str()),"Fit",700,600); 
	CName[Form("Fit_y%d_pt%d_%s",iY,iPT,MODE.c_str())]->cd();
	frame->Draw();
	
	//	TF1 *fitfunc = (TF1*)sum.asTF(RooArgList(m),RooArgList(A0,cw, cw2,deltaM, deltaM2,nbkg,nsig,r2,r3),m);
	//fitfunc->Draw(); 
	//	sum.fitTo(data_binned);
	//	sum.chi2FitTo(data_binned);
	
}
*/
void set_order(int order){
	
	A7.setConstant(kFALSE);
	A6.setConstant(kFALSE);
	A5.setConstant(kFALSE);
	A4.setConstant(kFALSE); 
	A3.setConstant(kFALSE); 
	A2.setConstant(kFALSE); 
	A1.setConstant(kFALSE); 
	
	p1.setConstant(kFALSE); 
	
	if(order==8){
		Aexp=1; 
		Acheb=0; 
		A7.setConstant(kTRUE);
		A6.setConstant(kTRUE);
		A5.setConstant(kTRUE);
		A4.setConstant(kTRUE); 
		A3.setConstant(kTRUE); 
		A2.setConstant(kTRUE); 
		A1.setConstant(kTRUE); 
		A7.setVal(0); 
		A6.setVal(0); 
		A5.setVal(0); 
		A4.setVal(0); 
		A3.setVal(0); 
		A2.setVal(0); 
		A1.setVal(0); 
		

	}
	
	if(order<=7){
		Acheb=1; 
		Aexp=0; 
		A7.setConstant(kTRUE); 
		A7.setVal(0); 
		p1.setConstant(kTRUE);
		p1.setVal(0); 
		
	}
	
	if(order<=5){
		A6.setConstant(kTRUE); 
		A6.setVal(0); 
	}
	
	if(order<=4){
		A5.setConstant(kTRUE); 
		A5.setVal(0); 
	}
	
	if(order<=3){
		A4.setConstant(kTRUE); 
		A4.setVal(0); 
	}
	
	if(order<=2){
		A3.setConstant(kTRUE); 
		A3.setVal(0); 
	}
	
	if(order<=1){
		A2.setConstant(kTRUE); 
		A2.setVal(0); 
	}
	if(order<=0){
		A1.setConstant(kTRUE); 
		A1.setVal(0); 
	}
	
}

void fill_yields(double *Y1,double *Y1E,double *Y2, double *Y2E, double *Y3, double *Y3E, double *r21, double *r21E,double *r31,double *r31E, double *prob){

	for (int order=0; order<=4; order++) {
		h_Name[yield_histogram(iY,1,Form("T%d",order),MODE)]->SetBinContent(iPT+1,Y1[order]); 
		h_Name[yield_histogram(iY,1,Form("T%d",order),MODE)]->SetBinError(iPT+1,Y1E[order]); 
		
		h_Name[ratio_histogram(iY,2,Form("T%d",order),MODE)]->SetBinContent(iPT+1,r21[order]); 
		h_Name[ratio_histogram(iY,2,Form("T%d",order),MODE)]->SetBinError(iPT+1,r21E[order]); 
		
		h_Name[ratio_histogram(iY,3,Form("T%d",order),MODE)]->SetBinContent(iPT+1,r31[order]); 
		h_Name[ratio_histogram(iY,3,Form("T%d",order),MODE)]->SetBinError(iPT+1,r31[order]); 
	}

	double average[3];
	string order="best"; 
	average_yield(Y1,Y1E,prob,average,"Y1S");
	h_Name[yield_histogram(iY,1,order,MODE)]->SetBinContent(iPT+1,average[0]); 
	h_Name[yield_histogram(iY,1,order,MODE)]->SetBinError(iPT+1,average[1]); 
	if(MODE=="") {
		h_Name[background_order_systematic(iY,"1S")]->SetBinContent(iPT+1,average[2]); 
	}
	
	average_yield(Y2,Y2E,prob,average,"Y2S");
	h_Name[yield_histogram(iY,2,order,MODE)]->SetBinContent(iPT+1,average[0]); 
	h_Name[yield_histogram(iY,2,order,MODE)]->SetBinError(iPT+1,average[1]); 
	if(MODE=="") {
		h_Name[background_order_systematic(iY,"2S")]->SetBinContent(iPT+1,average[2]); 
	}
	
	average_yield(Y3,Y3E,prob,average,"Y2S");
	h_Name[yield_histogram(iY,3,order,MODE)]->SetBinContent(iPT+1,average[0]); 
	h_Name[yield_histogram(iY,3,order,MODE)]->SetBinError(iPT+1,average[1]); 
	if(MODE=="") {
		h_Name[background_order_systematic(iY,"3S")]->SetBinContent(iPT+1,average[2]); 
	}
	
	average_yield(r21,r21E,prob,average,"R21");
	h_Name[ratio_histogram(iY,2,order,MODE)]->SetBinContent(iPT+1,average[0]); 
	h_Name[ratio_histogram(iY,2,order,MODE)]->SetBinError(iPT+1,average[1]); 
	if(MODE=="") {
		h_Name[background_order_systematic(iY,"R21")]->SetBinContent(iPT+1,average[2]); 
	}
	average_yield(r31,r31E,prob,average,"R31");
	h_Name[ratio_histogram(iY,3,order,MODE)]->SetBinContent(iPT+1,average[0]); 
	h_Name[ratio_histogram(iY,3,order,MODE)]->SetBinError(iPT+1,average[1]); 
	if(MODE=="") {
		h_Name[background_order_systematic(iY,"R31")]->SetBinContent(iPT+1,average[2]); 
	}

}

double CalcAIC(RooAddPdf pdf, RooFitResult *FR, RooDataHist data, int N, int k ){
	bool print = false; 
	//Take PDF and data and return AICc

	RooChi2Var chi2_var("chi2","chi2",pdf,data); 
	
	int offset=2*k+2*k*(k+1)/(N-k-1);
	//int offset=k+k*(k+1)/(N-k-1); 
	if(print) cout << "chi2: " << chi2_var.getVal() << " k " << k << " N " << N << endl; 
	
//	double AIC = chi2_var.getVal() + offset; 
	double AIC= -2*FR->minNll() + offset; 
	
	return AIC; 
}

void draw_header(){
	
	L1.SetNDC(kTRUE); 
	L2.SetNDC(kTRUE); 
	
	L1.DrawLatex(0.17,0.92, cms_pre); 
	L2.DrawLatex(0.38,0.93, lumi); 
	
}

void fit(){
	bool print = 1; 
	if(print) cout << " fit: " << endl; 
	if(print) cout << "FIT MODE: " << MODE << endl; 
	estimate_yield(hName[h_m(iY,iPT,MODE)]); 
	
	RooRealVar deltaM("deltaM","deltaM",0,-0.025,0.025);
	RooRealVar cw("cw","cw",1.0,0.75,1.75); 
	
	RooRealVar deltaM2("deltaM2","deltaM",0,-0.025,0.025);
	RooRealVar cw2("cw2","cw",1.0,0.75,1.25); 
	

	RooRealVar nsig("nsig","# sig", 1,0,5000000);
	
	RooRealVar nsig1("nsig1","# sig", 1,0,5000000);
	RooRealVar nsig2("nsig2","# sig", 1,0,5000000);
	RooRealVar nsig3("nsig3","# sig", 1,0,5000000);

	
	RooRealVar r2("r2","r2",0.5,0,1);
	RooRealVar r3("r3","r3",0.5,0,1);
	//Cheb coefficients: A0*T0+A1*T1...An*Tn where Tn is a Cheb poly and A0 is automatically included 
	
	
	TF1 *TF1_LS1=new TF1("TF1_LS1",signal_shape1, 8.5,11.5,2);
	TF1 *TF1_LS2=new TF1("TF1_LS2",signal_shape2, 8.5,11.5,2);
	TF1 *TF1_LS3=new TF1("TF1_LS3",signal_shape3, 8.5,11.5,2);
	
	TF1_LS1->SetParameters(0.0,1.0);
	TF1_LS2->SetParameters(0.0,1.0);
	TF1_LS3->SetParameters(0.0,1.0);
	
	
	RooAbsReal *LS1fcn = bindFunction(TF1_LS1,m, RooArgList(deltaM,cw));
	RooAbsReal *LS2fcn = bindFunction(TF1_LS2,m, RooArgList(deltaM,cw));
	RooAbsReal *LS3fcn = bindFunction(TF1_LS3,m, RooArgList(deltaM,cw));
	
	//define constant to make PDF
	RooRealVar tmp("tmp","tmp",0);
	RooRealVar c1("c1","c1",1);
	tmp.setConstant(kTRUE);
	RooChebychev K("K","K",m, RooArgList(tmp)); 
	
	RooAbsPdf *LS1Pdf = new RooRealSumPdf("LS1Pdf","LS1Pdf",*LS1fcn,K,c1);
	RooAbsPdf *LS2Pdf = new RooRealSumPdf("LS2Pdf","LS2Pdf",*LS2fcn,K,c1);
	RooAbsPdf *LS3Pdf = new RooRealSumPdf("LS3Pdf","LS3Pdf",*LS3fcn,K,c1);
	
	
	set_order(8);
	

	RooRealVar p2("p2","p2",-0.1,-100,100); 

	RooExponential exp("exp","bkg",m,p1); 
	
//	RooGenericPdf bkg("gexp","exp(p1*UM+p2*UM*UM)",RooArgSet(m,p1,p2));

	RooChebychev cheb("cheb","Background",m,RooArgSet(A1,A2,A3,A4,A5,A6,A7));
	
	RooAddPdf bkg("bkg","Cheb or exp", RooArgList(exp,cheb),RooArgList(RooConst(Aexp),RooConst(Acheb))); 
	
	RooRealVar nbkg("nbkg","# bkg", 1,0,5000000); 
	nsig.setVal(yield_[0]);
	nsig1.setVal(yield_[0]); 
	nsig2.setVal(yield_[1]); 
	nsig3.setVal(yield_[2]); 

	
	if(print){
		cout << "Y1 Estimate: " << yield_[0] << endl; 
		cout << "r21 Estimate: " << yield_[1]/yield_[0] << endl; 
		cout << "r31 Estimate: " << yield_[2]/yield_[0] << endl; 
	}
	
	r2.setVal(yield_[1]/yield_[0]); 
	r3.setVal(yield_[2]/yield_[0]); 
	m.setRange(8.7,11.3);
	RooAddPdf sig("signal","LS1+LS2+LS3", RooArgList(*LS1Pdf, *LS2Pdf, *LS3Pdf),RooArgList(RooConst(1),r2,r3));

	//RooArgSet obs(A1,A2,A3,A4,r1,r2,r3,nsig,nbkg);

	RooAddPdf sumR("sumR","s+bkg",RooArgList(sig, bkg),RooArgList(nsig,nbkg));
	RooAddPdf sum("sum","s+bkg",RooArgList(*LS1Pdf,*LS2Pdf,*LS3Pdf,bkg),RooArgList(nsig1,nsig2,nsig3,nbkg));
	

	double AIC[Ncheb+1]; 
	double AIC_r[Ncheb+1]; 
	double AIC_LL[Ncheb+1]; 

	double Y1[Ncheb+1];
	double Y2[Ncheb+1];
	double Y3[Ncheb+1];

	double r21[Ncheb+1];
	double r31[Ncheb+1]; 	
	double r21E[Ncheb+1];
	double r31E[Ncheb+1];
	
	double Y1E[Ncheb+1];
	double Y2E[Ncheb+1];
	double Y3E[Ncheb+1];
	
	double chi2[Ncheb+1];
	double Ndf[Ncheb+1];
	double prob[Ncheb+1]; 
	
	for(int order=0; order<=Ncheb; order++){
		if(print) cout << "order: " << order << endl; 
		RooPlot *plot = m.frame(Bins(m.getBins())); 
		plot->SetTitle(""); 
		if(print) cout << "RooPlot made: "<< endl; 
		
		if(print)cout << "HistName: " << dataHist_name(iY,iPT,MODE) << endl; 
		
		if(fPTbin[iPT]<50) dataHist[dataHist_name(iY,iPT,MODE)]->plotOn(plot); 
		else dataSet[dataSet_name(iY,iPT,MODE)]->plotOn(plot, Binning(m.getBins())); 
		if(print) cout << "plots successfully made: " << endl; 		
		
		CreateCanvas(fit_summary(iY,iPT,order,MODE),"",800,600); 
		int k=order+6; 
		int N=m.getBins(); 
		if(print) cout << "Number of bins/events: " << N << endl; 
		
		set_order(order);
		if(print) cout << "Perform fit: " << endl; 
		RooFitResult *fit_test(0);
		bool fit_qual=false;
		int PL=-1; 

		int Nfit=0; 
		do{
			if(fPTbin[iPT]<50)fit_test= sum.fitTo(*dataHist[dataHist_name(iY,iPT,MODE)],SumW2Error(kFALSE),Save(),Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(PL),Verbose(kFALSE));
			else fit_test= sum.fitTo(*dataSet[dataSet_name(iY,iPT,MODE)],SumW2Error(kFALSE),Save(),Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(PL),Verbose(kFALSE));
			if(fit_test->status()==0 && fit_test->covQual()==3) fit_qual=true;
			if(print || !fit_qual)cout << Form("Fit quality[%d]: ",order) << fit_qual << endl; 
			Nfit++; 
		}while (fit_qual==false && Nfit<10);

		if(print || !fit_qual)cout << Form("Fit quality[%d]: ",order) << fit_qual << endl;
		if(print) fit_test->Print(); 

		if(print) cout << "Y1: " << nsig.getVal() << " +/- " << nsig.getError() << endl; 
		if(print) cout << "r21: " << r2.getVal() << " +/- " << r2.getError() << endl; 
		if(print) cout << "r31: " << r3.getVal() << " +/- " << r3.getError() << endl; 
		
//		sum.SetName(fit_PDF(iY,iPT,MODE)); 
//		sum.Write(); 
		if(print) cout << "plot PDF for chi2: " << endl; 
		sum.plotOn(plot,LineColor(kBlack)); 
		
		TString name;
		if(fPTbin[iPT]<50) name=dataHist_name(iY,iPT,MODE); 
		else name=dataSet_name(iY,iPT,MODE);
		RooHist *hpull = plot->pullHist(); 

		sum.plotOn(plot,Components("bkg"), LineStyle(kDashed)); 
		sum.plotOn(plot,Components("LS1Pdf"),LineColor(kRed)); 
		sum.plotOn(plot,Components("LS2Pdf"),LineColor(kBlue)); 
		sum.plotOn(plot,Components("LS3Pdf"),LineColor(kGreen)); 		
		if(print)cout << "Create pull hist " << endl; 
		
		CName[fit_summary(iY,iPT,order,MODE)]->cd();
		if(print)cout << "Create new TPad " << endl; 
		double small = 0.1;
		CName[fit_summary(iY,iPT,order,MODE)]->Divide(1,2,small,small);
		double r  = .25;
        double sl = 1. / ( 1. - r );
		TPad* padHisto = (TPad*) CName[fit_summary(iY,iPT,order,MODE)]->cd(1);
        TPad* padResid = (TPad*) CName[fit_summary(iY,iPT,order,MODE)]->cd(2);

        padHisto->SetPad( 0., r , 1., 1. );
        padHisto->SetBottomMargin( small );
        padResid->SetPad( 0., 0., 1., r  );
        padResid->SetBottomMargin( r );
        padResid->SetTopMargin   ( small );
	
		TLatex label(0.55,0.75,BT(iY,iPT)); 
		label.SetNDC(kTRUE);
		TLatex order_label(0.55,0.65,background_shape(order)); 
		order_label.SetNDC(kTRUE); 
		
		RooPlot *pull_plot = m.frame(); 
		pull_plot->SetTitle(""); 
		pull_plot->SetYTitle("pull"); 
		
		pull_plot->addPlotable(hpull,"P"); 
		pull_plot->DrawClone(); 
		
		padHisto->cd(); 
		plot->DrawClone();
		label.DrawClone(); 
		order_label.DrawClone();
		draw_header();
		padResid->cd();
		pull_plot->DrawClone(); 
		
		AIC[order]=CalcAIC(sum,fit_test,*dataHist[dataHist_name(iY,iPT,MODE)],N,k); 
		
		Y1[order]=nsig1.getVal(); 
		Y2[order]=nsig2.getVal(); 
		Y3[order]=nsig3.getVal(); 
		
		Y1E[order]=nsig1.getError(); 
		Y2E[order]=nsig2.getError(); 
		Y3E[order]=nsig3.getError(); 

		if(print) cout <<"Now evaluate separate function N(1+r2+r3)+bkg " << endl; 
		
		do{
			if(fPTbin[iPT]<50)fit_test= sumR.fitTo(*dataHist[dataHist_name(iY,iPT,MODE)],SumW2Error(kFALSE),Save(),Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE));
			else fit_test= sumR.fitTo(*dataSet[dataSet_name(iY,iPT,MODE)],SumW2Error(kFALSE),Save(),Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE));
			if(fit_test->status()==0 && fit_test->covQual()==3) fit_qual=true;
			if(print || !fit_qual)cout << Form("Fit quality[%d]: ",order) << fit_qual << endl; 
			Nfit++; 
		}while (fit_qual==false && Nfit<10);
		if(print) fit_test->Print();
		
		AIC_r[order]=CalcAIC(sumR,fit_test,*dataHist[dataHist_name(iY,iPT,MODE)],N,k); 
		
		r21[order]=r2.getVal(); 
		r31[order]=r3.getVal(); 
		
		r21E[order]=r2.getError(); 
		r31E[order]=r3.getError(); 
		
		
		delete plot; 
		
	}//loop over orders, T0+A1T1+A2T2+A3T3+A4T4+...

	TF1_LS1->SetParameters(deltaM.getVal(),cw.getVal());
	if(print) cout << "Integral " << TF1_LS1->Integral(8.7,11.3) << endl;
	
	int best_order=compute_AIC_prob(AIC,prob); 
	h_Name[background_order(iY,"")]->SetBinContent(iPT+1,best_order); 
	if(print) cout << "best order AIC: " << best_order << " " <<  TMath::LocMin(5,AIC_LL) << endl; 
	for(int order=0; order<=Ncheb; order++){
		cout << "Yield: " << Y1[order] << endl; 
		if(order!=best_order){
			CName.erase(fit_summary(iY,iPT,order,MODE)); 
		}
		
	}
	CName[fit_summary(iY,iPT)]=CName[fit_summary(iY,iPT,best_order,MODE)]; 
	CName[fit_summary(iY,iPT)]->SetName(fit_summary(iY,iPT)); 
	CName.erase(fit_summary(iY,iPT,best_order,MODE)); 
	
	fill_yields(Y1,Y1E,Y2,Y2E,Y3,Y3E,r21,r21E, r31, r31E,prob); 
	
	if(print) cout << " Delete PDF pointers: " << endl; 

	delete LS1Pdf; 
	delete LS2Pdf; 
	delete LS3Pdf; 
	delete TF1_LS1;
	delete TF1_LS2; 
	delete TF1_LS3; 
	if(print) cout << "End of setup function: " << endl; 

}
							

void setup_fit(){
	
	/*
	cout << "Optimum Fit order: " << cheb_order<< endl; 
	cout << "Optimizer 2: " << order_opt2 << endl; 
	set_order(A1,A2,A3,A4,cheb_order);
	
	h_Name[background_order(iY,MODE)]->SetBinContent(iPT+1,cheb_order); 
	
	sum.fitTo(wdata,SumW2Error(kFALSE),Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE));
	FR=sum.fitTo(wdata,Save(),SumW2Error(kFALSE),Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE));

	FR->SetName(fit_result_name(iY,iPT,MODE));
	FR->Write();
	

	h_Name[yield_histogram(iY,1,MODE)]->SetBinContent(iPT+1,nsig.getVal()); 
	h_Name[yield_histogram(iY,1,MODE)]->SetBinError(iPT+1,nsig.getError()); 

	if(MODE=="w"){
		double Unc_uW=h_Name[yield_histogram(iY,1,"")]->GetBinError(iPT+1)/h_Name[yield_histogram(iY,1,"")]->GetBinContent(iPT+1);
		double Unc_corr=nsig.getError()/nsig.getVal()*TMath::Sqrt(bin_weight); 
		double diff=100*TMath::Abs(Unc_uW-Unc_corr); 
		h_Name[uncertainty_comparison(iY,1)]->SetBinContent(iPT+1,diff);
		h_Name[uncertainty_comparison(iY,1)]->SetBinError(iPT+1,0); 
	}
	
	const RooArgList *par = &FR->floatParsFinal(); 
	cout << "nsig: " << nsig.getVal() << " Percent Error: " << 100*nsig.getError()/nsig.getVal() << endl; 

	sum.plotOn(frame); 
	cout << "chi2: "<< frame->chiSquare() << endl; 
	 */

	
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


double G(double *x, double *par){
	double Sig=(par[0]/(1+par[3]+par[6]))*(TMath::Gaus(x[0],par[1],par[2],kTRUE)+par[3]*TMath::Gaus(x[0],par[4],par[5])+par[6]*TMath::Gaus(x[0],par[7],par[8]));
	double BKG=(1-par[0])/(par[9]+par[10]+par[11]+par[12]+par[13]);//*BKG_cheb(x,&par[9]); 
	return Sig+BKG; 
}

double pol(double *x, double *par){
	
	double y=par[0]*(1+par[1]*x[0])/(1.6+par[1]*2.56/2);
	if(x[0]>9.0 && x[0]<10.6){
		TF1::RejectPoint(); 
		return 0; 
	}
	return y; 
}



void get_LS(int iy, int ipt){
	cout << "get lineshape: " << endl; 
	string LS_name1; 
	string LS_name2;
	string LS_name3; 
	

		LS_name1=LS_name(1, iy,ipt); 
		LS_name2=LS_name(2, iy,ipt);
		LS_name3=LS_name(3, iy,ipt);

	
	
	cout <<"LS Name: " << LS_name1 << endl; 
	
	//LS1_shape=(TF2*)lineshape_file->FindObjectAny(LS_name1.c_str());
	lineshape_file = new TFile(Form("../Generate_LS/output_file_y%d_pt%d.root",iy,ipt),"READ"); 
	
	if(iy==0 && ipt==0) {
	
		TCanvas *C1 = (TCanvas*)lineshape_file->Get("M_mumu_fit_1S"); 
		CName[C1->GetName()]=C1;
		TCanvas *C2 = (TCanvas*)lineshape_file->Get("M_mumu_fit_2S"); 
		CName[C2->GetName()]=C2; 
		TCanvas *C3 = (TCanvas*)lineshape_file->Get("M_mumu_fit_3S"); 
		CName[C3->GetName()]=C3; 

	}
	
	LS1_shape=(TF2*)lineshape_file->FindObjectAny(LS_name1.c_str()); 
	LS2_shape=(TF2*)lineshape_file->FindObjectAny(LS_name2.c_str());
	LS3_shape=(TF2*)lineshape_file->FindObjectAny(LS_name3.c_str());
}

double signal_shape1(double *x, double *par){
	double m=x[0]; 
	double deltaM=par[0];
	double cw=par[1];
	return LS1_shape->Eval(m-deltaM, cw); 
}

double signal_shape2(double *x, double *par){
	double m=x[0]; 
	double deltaM=par[0];
	double cw=par[1];
	return LS2_shape->Eval(m-deltaM, cw); 
}						 

double signal_shape3(double *x, double *par){
	double m=x[0]; 
	double deltaM=par[0];
	double cw=par[1];
	return LS3_shape->Eval(m-deltaM, cw); 
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
	
	h_Name[name] = h;
}

void CreateHistogram(const char* name,   const char* title,
					 const char* xTitle, const char* yTitle,
					 Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
	TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);
	
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	
	h->Sumw2();
	
	h_Name[name] = h;
}

