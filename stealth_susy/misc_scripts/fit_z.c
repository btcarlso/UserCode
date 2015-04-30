/*
 *  fit_z.c
 *  
 *
 *  Created by Benjamin Carlson on 2/25/14.
 *  Copyright 2014 Carnegie Mellon University. All rights reserved.
 *
 */

#include "fit_z.h"

void fit_z(){
    gROOT->SetBatch();
    
	open_files();
 
    fit_bins();
    
    TFile *output_file = new TFile("fit_zpeak_out.root","RECREATE");
    
    for (std::map<TString,TCanvas*>::iterator it=CName.begin(); it!=CName.end(); it++) {
		CName[it->first]->Write();
        //TString name="/uscms/home/btcarlso/AN-13-083/tdr2/notes/AN-13-083/trunk/";
         TString name="/uscms_data/d3/btcarlso/figtmp/";
        name=name+it->first+".pdf";
        CName[it->first]->Print(name);
    }
    
    output_file->Close();
    
	
}

void fit_bins(){
    ///uscms/home/btcarlso/AN-13-083/tdr2/notes/AN-13-083/trunk/
    outputR.open("table_zfitR.tex");
    output.open("table_zfit_data.tex");
    output_injection.open("table_zfit_signal_injection.tex");
    output_mc.open("table_zfit_mc.tex");

    
    outputR << "\\begin{tabular}{";
	int N=1;
	for(int i=0; i<=N; i++) outputR << "c";
	outputR << "}" << endl;
    
    outputR << "n-jets & R \\\\ \\hline" << endl;
    
    output << "\\begin{tabular}{";
    N=7;
	for(int i=0; i<=N; i++) output << "c";
	output << "}" << endl;

    output << "n-jets & $N_{dy}^{MC}$ & $N_{top}^{MC}$ & $N_{sig}^{MC}$ & $N_{dy}^{fit}$ & $N_{diboson}$ & $N_{pol}$ & p1 \\\\ \\hline" << endl;
    
    
    output_injection << "\\begin{tabular}{";
	 N=7;
	for(int i=0; i<=N; i++) output_injection << "c";
	output_injection << "}" << endl;
    
    output_injection << "n-jets & $N_{dy}^{MC}$ & $N_{top}^{MC}$ & $N_{sig}^{MC}$ & $N_{dy}^{fit}$ & $N_{diboson}$ & $N_{pol}$ & p1 \\\\ \\hline" << endl;
    
    output_mc << "\\begin{tabular}{";
	 N=7;
	for(int i=0; i<=N; i++) output_mc << "c";
	output_mc << "}" << endl;
    
    output_mc << "n-jets & $N_{dy}^{MC}$ & $N_{top}^{MC}$ & $N_{sig}^{MC}$ & $N_{dy}^{fit}$ & $N_{diboson}$ & $N_{pol}$ & p1 \\\\ \\hline" << endl;
    
    
    double unc;
    for(int nJ=2; nJ<=6; nJ++){
        double R = fit_fraction(nJ,unc);
        //output << nJ << " & " ;
        //output << R << " $\\pm$ " << R*unc << "\\\\" << endl;
        
    }
    outputR << "\\end{tabular}" << endl;
    output << "\\end{tabular}" << endl;
    output_injection << "\\end{tabular}" << endl;
    output_mc <<"\\end{tabular}" << endl;
}

bool reject=true;
double pol_exc(double *x, double *par){
	if(reject && x[0]>70 && x[0]<105){
		TF1::RejectPoint();
		return 0;
	}
	return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

double integralE(TH1F *h,double x1,double x2){
	double tot=0;
	for(int i=1; i<=h->GetNbinsX();i++){
		double x=h->GetBinCenter(i);
		if(x>=x1 && x<=x2)tot+=h->GetBinError(i);
	}
	return tot; 
}

double integral(TH1F *h,double x1,double x2){
	double tot=0;
	for(int i=1; i<=h->GetNbinsX();i++){
		double x=h->GetBinCenter(i);
		if(x>=x1 && x<=x2)tot+=h->GetBinContent(i);
	}
	return tot; 
}

double fit_fraction(int nJet, double &unc){
	double hmin = 50;//hName["h_mumu_nJet5_dy"]->GetXaxis()->GetXmin();
	double hmax = 130;//hName["h_mumu_nJet5_dy"]->GetXaxis()->GetXmax();
		
	TString variable=Form("h_dilepton_mass_2Mu_nJets%d_stg300_0btag_SF",nJet);

	TString file="_dy";
	double ndy=integral(hName[variable+file],hmin,hmax); 
	file="_singleTop"; 
	double nt=integral(hName[variable+file],hmin,hmax); 
	file="_ttbar";
	double ntt=integral(hName[variable+file],hmin,hmax); 
	file="_diboson";
	double ndiboson=integral(hName[variable+file],hmin,hmax); 
	file="_stealth_300_200";
	double nsig=integral(hName[variable+file],hmin,hmax);
    file="_singleMu";
    double ndata=integral(hName[variable+file],hmin,hmax);
    
    cout << "Ndata: " << ndata << endl;
	cout << "Number of dy: " << ndy << " single Top: " << nt << " ttbar: " << ntt << " diboson: " << ndiboson << "nsig: " << nsig<< endl; 
	
	// Declare observable x
	RooRealVar mass("mass","mass",hmin,hmax) ;
	mass.setRange("RangeZ",81,101); 
	
	RooRealVar NDY("NDY","NDY",ndy,0,2500000);
	RooRealVar NT("NT","NT",nt,0,1000000);
	RooRealVar NTT("NTT","NTT",ntt,0,1000000);
	RooRealVar NDiboson("NDiboson","NDiboson",ndiboson,0,100000);
	RooRealVar NSig("NSig","NSig",nsig,0,1000000);
	RooRealVar Npol("Npol","Npol",ntt,0,1000000); 
	
	NT.setConstant(kTRUE);
	NDiboson.setConstant(kTRUE);
	//NTT.setConstant(kTRUE);

    
	file="_sigbkg";
	RooDataHist dh_sigbkg("dh_sigbkg","dh",mass,Import(*hName[variable+file])) ;
	
	file="_singleMu";
	RooDataHist dh_data("dh_data","dh",mass,Import(*hName[variable+file])) ;
	
	file="_allMC";
	RooDataHist dh_all("dh_all","dh",mass,Import(*hName[variable+file])) ;

	file="_ttbar";
	RooDataHist dh_ttbar("dh_ttbar","ttbar",mass,Import(*hName[variable+file])) ;
	file="_singleTop";
	RooDataHist dh_singleTop("dh_singleTop","singleTop",mass,Import(*hName[variable+file])) ;
	file="_diboson";
	RooDataHist dh_diboson("dh_diboson","diboson",mass,Import(*hName[variable+file])) ;
	file="_dy";
	RooDataHist dh_z("dh_z","dh_z",mass,Import(*hName[variable+file]));
	file="_stealth_300_200";
	RooDataHist dh_sig("dh_sig","dh_sig",mass,Import(*hName[variable+file]));

	//produce HistPDF
	
	//RooHistPdf ttbar_shape("TT_shape","TT_shape",mass,dh_ttbar,0);
	RooHistPdf singleTop_shape("T_shape","T_shape",mass,dh_singleTop,0);
	RooHistPdf diboson_shape("diboson_shape","diboson_shape",mass,dh_diboson,0); 
	RooHistPdf DY_shape("DY_shape","DY_shape",mass,dh_z,0); 
	//RooHistPdf sig_shape("sig_shape","sig_shape",mass,dh_sig,0);
	
	RooRealVar a1("a1","a1",0.5,-1,1); 
	
	RooChebychev pol("pol","pol",mass,RooArgList(a1)); 
	
	
	//RooAddPdf sum("sum","sum",RooArgList(DY_shape,singleTop_shape,diboson_shape,ttbar_shape),RooArgList(NDY,NT,NDiboson,NTT));
	//RooAddPdf sumSig("sumSig","sumSig",RooArgList(DY_shape,singleTop_shape,diboson_shape,ttbar_shape,sig_shape),RooArgList(NDY,NT,NDiboson,NTT,NSig));
	RooAddPdf sumPol("sumPol","sumPol",RooArgList(DY_shape,singleTop_shape,diboson_shape,pol), RooArgList(NDY,NT,NDiboson,Npol) );
    //RooAddPdf sumPol("sumPol","sumPol",RooArgList(DY_shape,diboson_shape,pol), RooArgList(NDY,NDiboson,Npol) );
    
	cout << "polynomial fit: " << endl; 
/*
	sumPol.fitTo(dh_all,Extended(kTRUE)); 
	cout << "Nz_fit: " << NDY.getVal() << endl; 
	cout << "polynomial fit sig+bkg: " << endl; 
	sumPol.fitTo(dh_sigbkg,Extended(kTRUE)); 
*/
//fit MC to itself
	sumPol.fitTo(dh_all,Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE));
//NDY,NT,NDiboson,Npol
	double Nz_fitMC=NDY.getVal();
	double Nz_fitMCE=NDY.getError();
    
    double NDB_fitMC=NDiboson.getVal();
    double NDB_fitMCE=NDiboson.getError();

    double Npol_fitMC=Npol.getVal();
    double Npol_fitMCE=Npol.getError();
    
	double p1MC=a1.getVal();
    double p1MCE=a1.getError();
    
//fit to signal + bkg
	sumPol.fitTo(dh_sigbkg,Extended(kTRUE),PrintEvalErrors(-1),PrintLevel(-1),Verbose(kFALSE));
	double Nz_fitSig=NDY.getVal();
	double Nz_fitSigE=NDY.getError();
    
    double NDB_fitSig=NDiboson.getVal();
    double NDB_fitSigE=NDiboson.getError();
    
    double Npol_fitSig=Npol.getVal();
    double Npol_fitSigE=Npol.getError();
	
    double p1Sig=a1.getVal();
    double p1SigE=a1.getError();
    
//fit to data
	sumPol.fitTo(dh_data,Extended(kTRUE));
	double Nz_fitData=NDY.getVal();
	double Nz_fitDataE=NDY.getError();
    
    double NDB_fitData=NDiboson.getVal();
    double NDB_fitDataE=NDiboson.getError();
    
    double Npol_fitData=Npol.getVal();
    double Npol_fitDataE=Npol.getError();
	
    double p1Data=a1.getVal();
    double p1DataE=a1.getError();
	
    CreateCanvas(Form("FitResult_nJets%d",nJet),"",600,600);
    CName[Form("FitResult_nJets%d",nJet)]->cd();
	
    gPad->SetLogy();
    RooPlot *plot_sig = mass.frame();
    dh_data.plotOn(plot_sig);
    sumPol.plotOn(plot_sig,LineColor(kBlack));
    sumPol.plotOn(plot_sig,Components("DY_shape"),LineColor(kBlue));
//    sumPol.plotOn(plot_sig,Components("T_shape"),LineColor(kGreen));
    sumPol.plotOn(plot_sig,Components("diboson_shape"),LineColor(kRed));
    sumPol.plotOn(plot_sig,Components("pol"),LineColor(kGreen),LineStyle(kDotted));
    //sumSig.plotOn(plot_sig,Components("sig_shape"),LineColor(kBlack),LineStyle(kDashed));
    plot_sig->DrawClone();

    TLegend *L = new TLegend(0.7,0.7,0.88,0.88);
	L->SetFillColor(10);
	L->SetLineColor(10);
	L->SetLineWidth(0);
    
    TF1 *ff = new TF1("ff","pol1",0,1);
    ff->SetLineColor(kBlue);
    TF1 *gg = new TF1("gg","pol1",0,1);
    gg->SetLineColor(kRed);
    TF1 *hh = new TF1("hh","pol1",0,1);
    hh->SetLineColor(kGreen);
    hh->SetLineStyle(kDotted);
    TF1 *ii = new TF1("ii","pol1",0,1);
    ii->SetLineColor(kBlack);
    
    file="_singleMu";
    hName[variable+file]->SetMarkerStyle(20);
    hName[variable+file]->SetMarkerSize(1);
    L->AddEntry(hName[variable+file],"Data", "EP");
    L->AddEntry(ff,"Drell-Yan","L");
    L->AddEntry(gg,"Diboson","L");
    L->AddEntry(hh,"Polynomial Background","L");
    L->AddEntry(ii,"sum bkg","L");
   // L->AddEntry(plot_sig->findObject("sig_shape"),"Stealth: M_{#tilde{q}}=300, M_{#chi^{#pm}}=200","L");

    L->DrawClone();
    draw_header();
    
//	cout << "Nz actual: " << ndy << " Nz fit MC: " << Nz_fitMC << " +/- " << Nz_fitMCE << "  Nz fit + signal: " << Nz_fit << " +/- " << Nz_fitE << " N sig: " << nsig << endl;
//	cout << "Nz data: " << Nz_fitData << " +/- " << Nz_fitDataE << endl;
    
    string pm=" $\\pm$ ";
    
    unc=Nz_fitDataE/Nz_fitData;
    
    outputR << nJet << " & " << Nz_fitData/ndy << " $\\pm$ " << (Nz_fitDataE/Nz_fitData)*(Nz_fitData/ndy) << "\\\\" << endl;
    
    output << nJet << setprecision(4) << " & " << ndy << " & " << ntt+nt << " & " <<"0" << " & "  << Nz_fitData << pm << Nz_fitDataE << " & ";
    output << NDB_fitData << pm << NDB_fitDataE << " & ";
    output << Npol_fitData << pm << Npol_fitDataE << " & ";
    output << p1Data << pm << p1DataE << "\\\\ " << endl;
    
    output_mc << nJet << setprecision(4) << " & " << ndy << " & " << ntt+nt << " & " <<"0" << " & " << Nz_fitMC << pm << Nz_fitMCE << " & ";
    output_mc << NDB_fitMC << pm << NDB_fitMCE << " & ";
    output_mc << Npol_fitMC << pm << Npol_fitMCE << " & ";
    output_mc << p1MC << pm << p1MCE << "\\\\ " << endl;
    
    output_injection << nJet << setprecision(4) << " & " << ndy << " & " << ntt+nt << " & " <<nsig << " & "  << Nz_fitSig << pm << Nz_fitSigE << " & ";
    output_injection << NDB_fitSig << pm << NDB_fitSigE << " & ";
    output_injection << Npol_fitSig << pm << Npol_fitSigE << " & ";
    output_injection << p1Sig << pm << p1SigE << "\\\\ " << endl;
    
    //    output_injection << nJet << " & " << ndy + ntt + nt + ndiboson << " & " <<  nsig << " & " << Nz_fit/ndy << " $\\pm$ " << (Nz_fitE/Nz_fit)*(Nz_fit/ndy) << "\\\\" << endl;
    
    delete ff;
    delete gg;
    delete hh;
    delete ii;
    
	return Nz_fitData/ndy;
	
}

void open_files(){
	//Open MC files
	
//	open_file("QCD");
	for(int irpv=300; irpv<=1000; irpv+=100){
		//	open_file(Form("RPV%d",irpv));
	}
	//open_file("UDD300"); 
/*
	open_file("stealth_500_400");
	open_file("stealth_300_200");
	open_file("wJets");
	open_file("ttbar");
	open_file("singleTop");
	open_file("diboson"); 
	open_file("singleMu");
 */
	open_file("diboson"); 
	open_file("singleTop");
	open_file("ttbar");
	open_file("sigbkg");
	open_file("allMC"); 
	open_file("dy");
	open_file("stealth_300_200");
	open_file("singleMu");
	
	
}


void open_file(TString name){
	cout << "open file:" << name << endl; 
	//TString file_name="/eos/uscms/store/user/btcarlso/histogram_files/output_file_"+name+".root"; 
	TString file_name="/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_"+name+".root";
	//if(name=="ttbar")file_name="/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_ttbar2.root";
	//if(name=="singleMu") file_name="/uscms_data/d3/btcarlso/TEST_DIR/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/output_file_allMC.root";
	TFile *f = new TFile(file_name,"READ");
	if(f->IsOpen()!=1) {
		cout << "File: " << file_name << " Failed to open!" << endl; 
		return; 
	}
	name="_"+name;
	fName[name]=f; 
	GetHistograms(name);
	
}

void GetHistograms(TString file){
	TH1F *TH1F_names=(TH1F*)fName[file]->FindObjectAny("TH1F_names"); 
	
	
	for (int i=1; i<=TH1F_names->GetNbinsX(); i++) {
		
		TString name=TH1F_names->GetXaxis()->GetBinLabel(i);
		if(name=="")continue; 
		if(name=="TH1F_names" || name=="TH2F_names" || name=="TProfile_names") continue; 
		//cout <<" filename: "<< file <<  " hist_name: "<< name << endl; 
		if(!(name.Contains("h_dilepton_mass")==1 || name.Contains("reconstructedZ")==1)) continue;
		TH1F *h=(TH1F*)fName[file]->FindObjectAny(name); 
		//if(file=="_singleMu" && name.Contains("st") && name.Contains("_GeV")==0 )setHistogramErrors(h); 
		h->SetStats(kFALSE); 
		h->SetName(h->GetName()+file); 
		hName[h->GetName()]=h; 
	}
	delete TH1F_names; 
	//cout << "end of GetHistograms() " << endl; 
}

void draw_header(){
	
	L1.SetNDC(kTRUE);
	L2.SetNDC(kTRUE);
	L1.SetTextSize(0.03);
	L2.SetTextSize(0.03);
	L1.DrawLatex(0.15,0.92, cms_pre);
	L2.DrawLatex(0.43,0.92, lumi);
	
}

void CreateCanvas(TString Name,TString Title, int x, int y ){
	TCanvas *createC = new TCanvas(Name, Title,x,y);
	CName[Name]=createC;
}
